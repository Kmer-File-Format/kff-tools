#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <queue>
#include <cassert>
#include <unordered_set>
#include <atomic>
#include <memory>
#include <thread>

#include "compact.hpp"
#include "encoding.hpp"
#include "sequences.hpp"
#include "merge.hpp"
#include "colinear.hpp"
#include "skmers.hpp"

using namespace std;


Compact::Compact() {
	input_filename = "";
	output_filename = "";
	sorted = false;
	nb_threads = 1;

	this->buffer_size = 1 << 10;
	this->next_free = 0;
	this->kmer_buffer = (uint8_t *)malloc(this->buffer_size);
	memset(this->kmer_buffer, 0, this->buffer_size);

	this->k = 0;
	this->m = 0;
	this->bytes_compacted = 0;
	this->offset_idx = 0;
}

Compact::Compact(std::string input_filename, std::string output_filename, bool sorted) {
    this->input_filename = input_filename;
    this->output_filename = output_filename;
    this->sorted = sorted;
    this->nb_threads = 1;

    this->buffer_size = 1 << 10;
    this->next_free = 0;
    this->kmer_buffer = (uint8_t *)malloc(this->buffer_size);
    memset(this->kmer_buffer, 0, this->buffer_size);

    this->k = 0;
    this->m = 0;
    this->bytes_compacted = 0;
    this->offset_idx = 0;
}


Compact::~Compact() {
	free(this->kmer_buffer);
}


void Compact::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("compact", "Read a kff file and try to compact the kmers from minimizer sections. The available ram must be sufficent to load a complete minimizer section into memory.");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "Input kff file to compact.");
	input_option->required();
	input_option->check(CLI::ExistingFile);
	CLI::Option * out_option = subapp->add_option("-o, --outfile", output_filename, "Kff to write (must be different from the input)");
	out_option->required();
	subapp->add_flag("-s, --sorted", sorted, "The output compacted superkmers will be sorted to allow binary search. Sorted superkmer have a lower compaction ratio (ie will be less compacted).");
	subapp->add_option("-t, --threads", nb_threads, "Number of threads: minimizer sections are compacted in parallel (the output file does not depend on this number). Default 1.");
}


void Compact::exec() {
	Kff_file infile(input_filename, "r");
	Kff_file outfile(output_filename, "w");

	outfile.write_encoding(infile.encoding);
	
	outfile.set_uniqueness(infile.uniqueness);
	outfile.set_canonicity(infile.canonicity);
	
	// Metadata transfer
	uint8_t * metadata = new uint8_t[infile.metadata_size];
	infile.read_metadata(metadata);
	outfile.write_metadata(infile.metadata_size, metadata);
	delete[] metadata;

	bool first_warning = true;

	// Minimizer sections are loaded by batches, compacted in parallel, then written in order.
	// A batch is written before any variable section, so the output does not depend on the
	// number of threads.
	if (this->nb_threads == 0)
		this->nb_threads = 1;
	vector<CompactSection> batch;
	size_t batch_bytes = 0;
	const size_t max_batch_sections = this->nb_threads == 1 ? 1 : 32 * this->nb_threads;
	const size_t max_batch_bytes = 256 << 20;
	vector<unique_ptr<Compact> > workers;
	for (uint t=0 ; t<this->nb_threads and this->nb_threads > 1 ; t++) {
		workers.emplace_back(new Compact());
		workers.back()->sorted = this->sorted;
	}

	auto flush = [&]() {
		if (batch.empty())
			return;
		if (workers.empty() or batch.size() == 1) {
			for (CompactSection & section : batch)
				this->compact_section(section);
		} else {
			atomic<size_t> next_section(0);
			auto work = [&](Compact * worker) {
				for (size_t i = next_section++ ; i < batch.size() ; i = next_section++)
					worker->compact_section(batch[i]);
			};
			vector<thread> threads;
			for (size_t t=0 ; t<workers.size() and t<batch.size() ; t++)
				threads.emplace_back(work, workers[t].get());
			for (thread & th : threads)
				th.join();
		}
		for (CompactSection & section : batch) {
			outfile.register_position('m');
			outfile.write(section.output.data(), section.output.size());
		}
		batch.clear();
		batch_bytes = 0;
	};

	while (infile.tellp() != infile.end_position) {
		char section_type = infile.read_section_type();

		if (section_type == 'v') {
			flush();
			Section_GV isgv(&infile);
			isgv.close();

			unordered_map<string, uint64_t> to_copy;
			for (auto & p : isgv.vars) {
				if (p.first != "first_index" and p.first != "footer_size") {
					to_copy[p.first] = p.second;
				}
			}

			if (to_copy.size() > 0) {
				Section_GV osgv(&outfile);
				for (auto & p : isgv.vars)
					osgv.write_var(p.first, p.second);
				osgv.close();
			}
		}
		else if (section_type == 'i') {
			Section_Index si(&infile);
			si.close();
		}
		else if (section_type == 'r') {
			if (first_warning) {
				first_warning = false;
				cerr << "WARNING: kff-tools has detected R sections inside of the file. The compact tool is only compacting kmers inside of M sections. The R sections are omitted." << endl;
			}

			Section_Raw sr(&infile);
			sr.close();
		}
		else if (section_type == 'm') {
			uint k = outfile.global_vars["k"];
			uint m = outfile.global_vars["m"];

			// Rewrite a value section if max is not sufficently large
			if (outfile.global_vars["max"] < k - m + 1) {
				flush();
				unordered_map<string, uint64_t> values(outfile.global_vars);
				Section_GV sgv(&outfile);

				for (auto & p : values)
					if (p.first != "max")
						sgv.write_var(p.first, p.second);
				sgv.write_var("max", k - m + 1);

				sgv.close();
			}

			// Load the section (compacted and saved with its batch)
			Section_Minimizer sm(&infile);
			batch.emplace_back();
			this->load_section(sm, outfile, batch.back());
			sm.close();
			batch_bytes += batch.back().bytes.size();
			if (batch.size() >= max_batch_sections or batch_bytes >= max_batch_bytes)
				flush();
		}
	}
	flush();

	infile.close();
	outfile.close();
}

void Compact::compact_section(Section_Minimizer & ism, Kff_file & outfile) {
	CompactSection section;
	this->load_section(ism, outfile, section);
	this->compact_section(section);
	outfile.register_position('m');
	outfile.write(section.output.data(), section.output.size());
}


void Compact::load_section(Section_Minimizer & ism, Kff_file & outfile, CompactSection & section) const {
	// Variables of the input section
	section.k = ism.k;
	section.m = ism.m;
	section.max = ism.max;
	section.data_size = ism.data_size;
	// Variables of the compacted section
	section.out_k = outfile.global_vars["k"];
	section.out_m = outfile.global_vars["m"];
	section.out_max = outfile.global_vars["max"];
	section.out_data_size = outfile.global_vars["data_size"];

	section.minimizer.assign(ism.minimizer, ism.minimizer + (2 * ism.m + 7) / 8);

	// Blocks: sequence without minimizer then data
	uint64_t max_nucl = ism.k + ism.max - 1 - ism.m;
	uint64_t max_seq_bytes = (max_nucl + 3) / 4;
	vector<uint8_t> seq_buffer(max_seq_bytes + 1);
	vector<uint8_t> data_buffer(ism.data_size * ism.max + 1);
	section.blocks.clear();
	section.bytes.clear();
	section.blocks.reserve(ism.nb_blocks);
	for (uint64_t n=0 ; n<ism.nb_blocks ; n++) {
		uint64_t mini_pos = 0xFFFFFFFFFFFFFFFF;
		uint nb_kmers = ism.read_compacted_sequence_without_mini(seq_buffer.data(), data_buffer.data(), mini_pos);
		uint64_t seq_bytes = (ism.k - ism.m + nb_kmers - 1 + 3) / 4;
		section.blocks.push_back({nb_kmers, mini_pos, section.bytes.size()});
		section.bytes.insert(section.bytes.end(), seq_buffer.data(), seq_buffer.data() + seq_bytes);
		section.bytes.insert(section.bytes.end(), data_buffer.data(), data_buffer.data() + nb_kmers * ism.data_size);
	}
}


void Compact::compact_section(CompactSection & section) {
	// General variables
	uint k = section.out_k;
	uint m = section.out_m;
	uint data_size = section.out_data_size;

	this->k = k;
	this->m = m;
	this->data_size = data_size;
	this->bytes_compacted = (k - m + 3) / 4;
    this->mini_pos_size = (static_cast<uint>(ceil(log2(section.k + section.max - section.m))) + 7) / 8;
	this->offset_idx = (4 - ((k - m) % 4)) % 4;


    // 1 - Load the input section (the kmers of the previous sections are not needed any more)
	this->next_free = 0;
	vector<vector<uint8_t *> > kmers_per_index = this->prepare_kmer_matrix(section);
	// write_paths() reads the minimizer positions with the width of the output section, which
	// can be larger than the width used in the buffer: the bytes after the last kmer are then
	// read. They are kept null, as they were when the buffer was never emptied (the loading
	// keeps 8 free bytes at the end of the buffer).
	memset(this->kmer_buffer + this->next_free, 0, 8);
	
	// 2 - Compact kmers
	vector<vector<uint8_t *> > paths;
	if (this->sorted) {
		paths = this->sorted_assembly(kmers_per_index);
	} else {
		vector<pair<uint8_t *, uint8_t *> > to_compact = this->greedy_assembly(kmers_per_index);
		paths = this->pairs_to_paths(to_compact);
	}

	// 3 - Bytes of the compacted section: type, minimizer, number of blocks, blocks
	MinimizerBlockEncoder encoder(k, m, section.out_max, data_size);
	section.output.clear();
	section.output.push_back('m');
	section.output.insert(section.output.end(), section.minimizer.begin(), section.minimizer.begin() + encoder.minimizer_bytes);
	const size_t nb_blocks_position = section.output.size();
	section.output.resize(nb_blocks_position + 8, 0);
	uint64_t nb_blocks = this->write_paths(paths, encoder, section.output);
	for (uint b=0 ; b<8 ; b++)
		section.output[nb_blocks_position + b] = static_cast<uint8_t>(nb_blocks >> (8 * (7 - b)));
}


long Compact::add_kmer_to_buffer(const uint8_t * seq, const uint8_t * data, uint64_t mini_pos) {
	// Realloc if needed
	if (this->buffer_size - this->next_free < this->bytes_compacted + this->data_size + this->mini_pos_size) {
		this->kmer_buffer = (uint8_t *) realloc((void *)this->kmer_buffer, this->buffer_size*2);
		memset(this->kmer_buffer + this->buffer_size, 0, this->buffer_size);
		this->buffer_size *= 2;
	}
	long position = this->next_free;
	uint8_t * buffer = this->kmer_buffer + this->next_free;

	// Copy kmer sequence
	memcpy(buffer, seq, this->bytes_compacted);
	buffer += this->bytes_compacted;
	// Copy data array
	memcpy(buffer, data, this->data_size);
	buffer += this->data_size;
	// Write mini position
	for (int b=mini_pos_size-1 ; b>=0 ; b--) {
		*(buffer + b) = mini_pos & 0xFF;
		mini_pos >>= 8;
	}
	this->next_free += this->bytes_compacted + this->data_size + this->mini_pos_size;

	return position;
}

vector<vector<uint8_t *> > Compact::prepare_kmer_matrix(Section_Minimizer & sm) {
	CompactSection section;
	section.k = sm.k;
	section.m = sm.m;
	section.max = sm.max;
	section.data_size = sm.data_size;
	uint64_t max_nucl = sm.k + sm.max - 1 - sm.m;
	vector<uint8_t> seq_buffer((max_nucl + 3) / 4 + 1);
	vector<uint8_t> data_buffer(sm.data_size * sm.max + 1);
	for (uint64_t n=0 ; n<sm.nb_blocks ; n++) {
		uint64_t mini_pos = 0xFFFFFFFFFFFFFFFF;
		uint nb_kmers = sm.read_compacted_sequence_without_mini(seq_buffer.data(), data_buffer.data(), mini_pos);
		uint64_t seq_bytes = (sm.k - sm.m + nb_kmers - 1 + 3) / 4;
		section.blocks.push_back({nb_kmers, mini_pos, section.bytes.size()});
		section.bytes.insert(section.bytes.end(), seq_buffer.data(), seq_buffer.data() + seq_bytes);
		section.bytes.insert(section.bytes.end(), data_buffer.data(), data_buffer.data() + nb_kmers * sm.data_size);
	}
	return this->prepare_kmer_matrix(section);
}


vector<vector<uint8_t *> > Compact::prepare_kmer_matrix(const CompactSection & sm) {
	vector<vector<long> > pos_matrix(sm.k - sm.m + 1, vector<long>());

	uint64_t kmer_bytes = (sm.k - sm.m + 3) / 4;
    uint64_t mini_pos_size = (static_cast<uint>(ceil(log2(sm.k + sm.max - sm.m))) + 7) / 8;

	// 1 - Load the input section
	for (const CompactSection::Block & block : sm.blocks) {
		uint64_t mini_pos = block.mini_pos;
		uint nb_kmers = block.nb_kmers;
		const uint8_t * seq_buffer = sm.bytes.data() + block.offset;
		const uint8_t * data_buffer = seq_buffer + (sm.k - sm.m + nb_kmers - 1 + 3) / 4;

		// Add kmer by index
		for (uint kmer_idx=0 ; kmer_idx<nb_kmers ; kmer_idx++) {
			uint kmer_pos = sm.k - (uint)sm.m - mini_pos + kmer_idx;

			// Realloc if needed (8 bytes are kept free after the last kmer, see compact_section)
			if (this->buffer_size - this->next_free < kmer_bytes + sm.data_size + mini_pos_size + 8) {
				this->kmer_buffer = (uint8_t *) realloc((void *)this->kmer_buffer, this->buffer_size*2);
				memset(this->kmer_buffer + this->buffer_size, 0, this->buffer_size);
				this->buffer_size *= 2;
			}

			// Copy kmer sequence
			subsequence(seq_buffer, sm.k - sm.m + nb_kmers - 1, this->kmer_buffer + next_free, kmer_idx, kmer_idx + sm.k - sm.m - 1);
			// Copy data array
			memcpy(this->kmer_buffer + next_free + kmer_bytes, data_buffer + kmer_idx * sm.data_size, sm.data_size);
			// Write mini position
			uint kmer_mini_pos = mini_pos - kmer_idx;

            for (int b=mini_pos_size-1 ; b>=0 ; b--) {
				*(this->kmer_buffer + next_free + kmer_bytes + sm.data_size + b) = kmer_mini_pos & 0xFF;
				kmer_mini_pos >>= 8;
			}
			// Update
			pos_matrix[kmer_pos].push_back(this->next_free);
			next_free += kmer_bytes + sm.data_size + mini_pos_size;
		}
	}

	// Transform the position matrix into the kmer matrix
	vector<vector<uint8_t *> > kmer_matrix(sm.k - sm.m + 1, vector<uint8_t *>());
	size_t idx = 0;

	for (vector<long> & positions : pos_matrix) {
		for (long pos : positions) {
			kmer_matrix[idx].push_back(this->kmer_buffer + pos);
		}
		idx += 1;
	}

	return kmer_matrix;
}


uint Compact::mini_pos_from_buffer(const uint8_t * kmer) const {
	uint mini_pos = 0;

	for (uint i=0 ; i<this->mini_pos_size ; i++) {
		mini_pos <<= 8;
		mini_pos += kmer[this->bytes_compacted + this->data_size + i];
	}

	return mini_pos;
}

int Compact::interleaved_compare_kmers(const uint8_t * kmer1, const uint8_t * kmer2) const {
	const uint mini_pos1 = this->mini_pos_from_buffer(kmer1);
	[[maybe_unused]] const uint mini_pos2 = this->mini_pos_from_buffer(kmer2);

	assert(mini_pos1 == mini_pos2);

	const uint used_nucl = this->k - this->m;
	const uint offset_nucl = (4 - (used_nucl % 4)) % 4;
	const uint pref_nucl = mini_pos1;
	const uint pref_bytes = (offset_nucl + pref_nucl + 3) / 4;
	const uint suff_nucl = used_nucl - pref_nucl;
	const uint suff_bytes = (suff_nucl + 3) / 4;
	const uint total_bytes = (used_nucl + 3) / 4;

	// --- Prefix ---
	int last_prefix_divergence = -1;
	// Prepare masks
	const uint pref_start_mask = (1u << (8 - 2 * offset_nucl)) - 1;
	const uint pref_stop_mask = ~((1u << (2 * (suff_nucl % 4))) - 1);
	// Iterate over all bytes
	for (uint pref_byte=0 ; pref_byte<pref_bytes ; pref_byte++) {
		uint8_t byte1 = kmer1[pref_byte];
		uint8_t byte2 = kmer2[pref_byte];

		bool end_byte = pref_byte == pref_bytes-1;

		// Mask useless bits
		if (pref_byte == 0) {
			byte1 &= pref_start_mask;
			byte2 &= pref_start_mask;
		}
		if (end_byte) {
			byte1 &= pref_stop_mask;
			byte2 &= pref_stop_mask;	
		}

		// Compare
		uint8_t result = byte1 xor byte2;
		// Get the rightmost bit set to 0: ie the last difference between sequences
		for (uint8_t i=0 ; i<4 ; i++) {
			if ((not end_byte) or (end_byte and (i >= (suff_nucl % 4)))) {
				if ((result & (0b11 << (2 * i))) != 0) {
					last_prefix_divergence = pref_byte * 4 + 3 - i - offset_nucl;	
					break;
				}
			}
		}
	}

	// --- Suffix ---
	int first_suffix_divergence = suff_nucl;
	int current_divergence_idx = 0;
	const uint suff_first_byte = total_bytes - suff_bytes;
	
	// Iterate over all the bytes from the suffix
	for (uint suff_byte=suff_first_byte ; suff_byte<total_bytes and first_suffix_divergence==(int)suff_nucl ; suff_byte++) {
		// Extract and compare bytes

		uint8_t byte1 = kmer1[suff_byte];
		uint8_t byte2 = kmer2[suff_byte];
		uint8_t result = byte1 xor byte2;

        for (uint8_t i=4 ; i>0 ; i--) {

			// Skip the first nucleotides of the first byte
            if (suff_byte == suff_first_byte and i == 4) {
                i = ((suff_nucl - 1) % 4) + 1;
			}

			// Check for divergeance
			if ((result & (0b11 << (2 * (i-1)))) != 0) {
				first_suffix_divergence = current_divergence_idx;
                break;
			} else {
				current_divergence_idx += 1;
			}
		}
	}

	// In case of sequence similarity
	if (last_prefix_divergence == -1 and first_suffix_divergence == (int)suff_nucl) {
		return 0;
	}
	// Check the first suffix divergence
	else {
		// Compute the divergence posi
		uint pref_div_distance = pref_nucl - last_prefix_divergence - 1;
		uint nucl_pos = offset_nucl;

		if (pref_div_distance == pref_nucl) {
			nucl_pos += pref_nucl + first_suffix_divergence;
		} else if (first_suffix_divergence == (int)suff_nucl) {
			nucl_pos += last_prefix_divergence;
		}
		// First interleaved divergence in the prefix
		else if ((int)pref_div_distance < first_suffix_divergence) {
			nucl_pos += last_prefix_divergence;
		}
		// First interleaved divergence in the suffix
		else {
			nucl_pos += pref_nucl + first_suffix_divergence;
		}


		// Extract the divergent nucleotides
		uint byte_pos = nucl_pos / 4;
		uint nucl_shift = 2 * (3 - (nucl_pos % 4));
		uint nucl1 = (kmer1[byte_pos] >> nucl_shift) & 0b11;
		uint nucl2 = (kmer2[byte_pos] >> nucl_shift) & 0b11;

		// Compare
		if (nucl1 < nucl2)
			return -1;
		else
			return +1;
	}
}


void Compact::interleaved_key(const uint8_t * kmer, uint64_t & high, uint64_t & low) const {
	const uint nb_nucl = this->k - this->m;
	const uint offset = (4 - (nb_nucl % 4)) % 4;
	interleaved_value(kmer, offset, nb_nucl, this->mini_pos_from_buffer(kmer), high, low);
}


void Compact::sort_matrix(vector<vector<uint8_t *> > & kmer_matrix) {
	// Sort by column

	// Keys of the comparison: same comparisons, so std::sort gives the same order
	struct Keyed {
		uint64_t high;
		uint64_t low;
		uint8_t * kmer;
	};
	vector<Keyed> keyed;

	for (uint i=0 ; i<kmer_matrix.size() ; i++) {
		if (this->k - this->m > 64) {
			// Comparison function (depends on minimizer position)
			auto comp_function = [this](const uint8_t * kmer1, const uint8_t * kmer2) {
				return this->interleaved_compare_kmers(kmer1, kmer2) < 0;
			};
			sort(kmer_matrix[i].begin(), kmer_matrix[i].end(), comp_function);
			continue;
		}

		keyed.resize(kmer_matrix[i].size());
		for (size_t j=0 ; j<keyed.size() ; j++) {
			keyed[j].kmer = kmer_matrix[i][j];
			this->interleaved_key(keyed[j].kmer, keyed[j].high, keyed[j].low);
		}
		sort(keyed.begin(), keyed.end(), [](const Keyed & a, const Keyed & b) {
			return a.high != b.high ? a.high < b.high : a.low < b.low;
		});
		for (size_t j=0 ; j<keyed.size() ; j++)
			kmer_matrix[i][j] = keyed[j].kmer;
	}
}

#include <bitset>
/* The nb_nucl <= 64 nucleotides of a kmer of the buffer (aligned on the right of its bytes) as a
 * 128-bit integer (high, low), the last nucleotide in the 2 lowest bits, without the padding bits.
 */
static inline void kmer_value(const uint8_t * kmer, uint nb_nucl, uint64_t & high, uint64_t & low) {
	high = low = 0;
	for (uint b=0 ; b<(nb_nucl + 3) / 4 ; b++) {
		high = (high << 8) | (low >> 56);
		low = (low << 8) | kmer[b];
	}
	if (nb_nucl < 32)
		low &= (uint64_t(1) << (2 * nb_nucl)) - 1;
	else if (nb_nucl < 64)
		high &= (uint64_t(1) << (2 * (nb_nucl - 32))) - 1;
}


vector<pair<uint64_t, uint64_t> > Compact::pair_kmers(const vector<uint8_t *> & column1, const vector<uint8_t *> & column2) const {
	const uint nb_nucl = k - m;
    vector<pair<uint64_t, uint64_t> > pairs;
	pairs.reserve(max(column1.size(), column2.size()));

	// Hashes of the k-m-1 prefixes (second column) and suffixes (first column): their value
	// (subseq_to_uint() truncates them to their 32 last nucleotides). From 2 to 64 nucleotides,
	// they are extracted from the value of the kmer. When the prefixes have at most 32
	// nucleotides, the hash is the whole prefix: kmers with the same hash overlap.
	const bool from_value = nb_nucl >= 2 and nb_nucl <= 64;
	const bool exact_hash = from_value and nb_nucl - 1 <= 32;
	const uint64_t sub_mask = nb_nucl - 1 >= 32 ? ~uint64_t(0) : (uint64_t(1) << (2 * (nb_nucl - 1))) - 1;
	auto prefix_hash = [&](const uint8_t * kmer) -> uint64_t {
		if (not from_value)
			return subseq_to_uint(kmer, nb_nucl, 0, nb_nucl-2);
		uint64_t high, low;
		kmer_value(kmer, nb_nucl, high, low);
		return ((low >> 2) | (high << 62)) & sub_mask;
	};
	auto suffix_hash = [&](const uint8_t * kmer) -> uint64_t {
		if (not from_value)
			return subseq_to_uint(kmer, nb_nucl, 1, nb_nucl-1);
		uint64_t high, low;
		kmer_value(kmer, nb_nucl, high, low);
		return low & sub_mask;
	};

	// Index the second column by their prefix hash: (hash, kmer position) sorted by hash, then
	// by position (the candidates of a hash are tested in the column order)
	vector<pair<uint64_t, uint64_t> > index;
	index.reserve(column2.size());
	for (uint64_t idx=0 ; idx<column2.size() ; idx++)
		index.emplace_back(prefix_hash(column2[idx]), idx);
	// By hash, then by position
	auto by_hash = [](const pair<uint64_t, uint64_t> & a, const pair<uint64_t, uint64_t> & b) {
		return a.first != b.first ? a.first < b.first : a.second < b.second;
	};
	sort(index.begin(), index.end(), by_hash);


	// Looks for suffix matches of the first column
	for (uint64_t idx=0 ; idx<column1.size() ; idx++) {
		uint8_t * kmer = column1[idx];
		uint64_t hash = suffix_hash(kmer);

		// Test each candidate with the same hash
		auto candidate = lower_bound(index.begin(), index.end(), make_pair(hash, uint64_t(0)), by_hash);
		for ( ; candidate != index.end() and candidate->first == hash ; candidate++) {
			uint64_t candidate_vect_idx = candidate->second;
            if (exact_hash or sequence_compare(
						column2[candidate_vect_idx], nb_nucl, 0, nb_nucl-2,
						kmer, nb_nucl, 1, nb_nucl-1
					) == 0) {
				pairs.emplace_back(idx, candidate_vect_idx);
			}
		}
	}

	return pairs;
}





ostream& operator<<(ostream& os, const PairInt& p)
{
    os << p.first << ", " << p.second;
    return os;
}



vector<pair<uint64_t, uint64_t> > Compact::colinear_chaining(vector<pair<uint64_t, uint64_t> > & candidates) const {

    if (candidates.empty()) {
        return {};
    }

	// Candidates increasing on both sides (frequent): they are the only longest chain
	bool increasing = true;
	for (size_t i=1 ; i<candidates.size() and increasing ; i++)
		increasing = candidates[i].first > candidates[i-1].first and candidates[i].second > candidates[i-1].second;
	if (increasing)
		return candidates;

	return colinear_longest_chain(candidates);
}


vector<vector<uint8_t *> > Compact::polish_sort(const vector<vector<uint8_t *> > & matrix , const vector<vector<pair<uint64_t, uint64_t> > > & pairs) const {
	const size_t nb_cols = matrix.size();
	const uint64_t none = numeric_limits<uint64_t>::max();

	// Kmer ids: the kmer of the row r of the column c is col_start[c] + r
	vector<uint64_t> col_start(nb_cols + 1, 0);
	for (size_t col=0 ; col<nb_cols ; col++)
		col_start[col + 1] = col_start[col] + matrix[col].size();
	const uint64_t nb_kmers = col_start[nb_cols];

	// Unify kmers into skmers. The links between 2 columns are a colinear chain (increasing on
	// both sides): a kmer is in at most one link on each side and the skmers are the paths of
	// linked kmers, from column to column (lonely kmers are skmers of 1 kmer).
	vector<uint64_t> next_kmer(nb_kmers, none);
	vector<bool> has_prev(nb_kmers, false);
	for (size_t col=0 ; col<pairs.size() ; col++)
		for (const pair<uint64_t, uint64_t> & link : pairs[col]) {
			next_kmer[col_start[col] + link.first] = col_start[col + 1] + link.second;
			has_prev[col_start[col + 1] + link.second] = true;
		}
	// Skmer of each kmer, first kmer, its column and size of each skmer
	vector<uint64_t> sk_of(nb_kmers);
	vector<uint64_t> sk_first, sk_col, sk_size;
	for (size_t col=0 ; col<nb_cols ; col++)
		for (uint64_t id=col_start[col] ; id<col_start[col + 1] ; id++) {
			if (has_prev[id])
				continue;
			const uint64_t sk = sk_first.size();
			uint64_t size = 0;
			for (uint64_t kmer=id ; kmer!=none ; kmer=next_kmer[kmer]) {
				sk_of[kmer] = sk;
				size += 1;
			}
			sk_first.push_back(id);
			sk_col.push_back(col);
			sk_size.push_back(size);
		}
	vector<bool>().swap(has_prev);

	// Number of kmers of each skmer that are the first remaining kmer of their column (0 for the
	// skmers already selected). A skmer can be selected when all its kmers are first.
	vector<uint64_t> sk_counts(sk_first.size(), 0);
	// Row of the first remaining kmer of each column
	vector<size_t> current_kmers(nb_cols, 0);

	// Interleaved values of the first remaining kmers: the kmers of the column c have their
	// minimizer after num_kmer_nucl - c nucleotides (the bytes are read from their first bit, as
	// interleaved() does). Up to 64 nucleotides, they are 128-bit integers aligned on the left.
	const size_t num_kmer_nucl = this->k - this->m;
	const bool use_values = num_kmer_nucl <= 64;
	vector<pair<uint64_t, uint64_t> > current_values(use_values ? nb_cols : 0);
	vector<uint8_t *> memory;
	vector<interleved_t> current_interleaves;
	if (not use_values)
		for (size_t col=0 ; col<nb_cols ; col++) {
			memory.push_back(new uint8_t[(num_kmer_nucl + 3) / 4]);
			current_interleaves.push_back((interleved_t){0,0,0});
		}
	auto set_current = [&](size_t col) {
		const uint64_t id = col_start[col] + current_kmers[col];
		sk_counts[sk_of[id]] += 1;
		const uint8_t * kmer = matrix[col][current_kmers[col]];
		if (not use_values) {
			current_interleaves[col] = interleaved(kmer, memory[col], num_kmer_nucl, num_kmer_nucl - col);
			return;
		}
		uint64_t high, low;
		interleaved_value(kmer, 0, num_kmer_nucl, num_kmer_nucl - col, high, low);
		const uint shift = 128 - 2 * num_kmer_nucl;
		if (shift >= 64) {
			high = low << (shift - 64);
			low = 0;
		} else if (shift > 0) {
			high = (high << shift) | (low >> (64 - shift));
			low <<= shift;
		}
		current_values[col] = make_pair(high, low);
	};
	// inf_interleaved(): the values are compared on the interleaved part of the 2 kmers (the
	// shortest one), then on the whole kmers if they have the same prefix size.
	auto inter_size = [num_kmer_nucl](size_t pref) -> size_t {
		const size_t suf = num_kmer_nucl - pref;
		return pref < suf ? 2 * pref + 1 : 2 * suf;
	};
	auto inf_current = [&](size_t col1, size_t col2) -> bool {
		if (not use_values)
			return inf_interleaved(current_interleaves[col1], current_interleaves[col2]);
		uint64_t high1 = current_values[col1].first, low1 = current_values[col1].second;
		uint64_t high2 = current_values[col2].first, low2 = current_values[col2].second;
		if (col1 != col2) {
			const size_t bits = 2 * min(inter_size(num_kmer_nucl - col1), inter_size(num_kmer_nucl - col2));
			if (bits <= 64) {
				const uint64_t mask = bits == 0 ? 0 : ~uint64_t(0) << (64 - bits);
				high1 &= mask; high2 &= mask;
				low1 = low2 = 0;
			} else {
				const uint64_t mask = ~uint64_t(0) << (128 - bits);
				low1 &= mask; low2 &= mask;
			}
		}
		return high1 != high2 ? high1 < high2 : low1 < low2;
	};

	for (size_t col=0 ; col<nb_cols ; col++)
		if (matrix[col].size() > 0)
			set_current(col);

	// Columns in the order of the candidates: the one in the middle first, then alternate
	// left-right
	vector<size_t> candidate_order(nb_cols);
	for (size_t i=0 ; i<nb_cols ; i++)
		candidate_order[i] = i % 2 == 0 ? nb_cols / 2 + i / 2 : nb_cols / 2 - 1 - i / 2;

	vector<vector<uint8_t *> > result;
	result.reserve(sk_first.size());

	// Add iteratively the skmers into the sorted list
	uint64_t remaining = nb_kmers;
	while (remaining > 0) {
		// Select the min first kmer whose skmer has all its kmers first (min_interleaved())
		size_t selected = nb_cols;
		for (size_t position : candidate_order) {
			if (current_kmers[position] >= matrix[position].size())
				continue;
			const uint64_t sk = sk_of[col_start[position] + current_kmers[position]];
			if (sk_counts[sk] == sk_size[sk] and (selected == nb_cols or inf_current(position, selected)))
				selected = position;
		}
		if (selected == nb_cols)
			throw runtime_error("compact: no superkmer can be selected (the links are not colinear)");

		// Add the related skmer to the output
		const uint64_t sk = sk_of[col_start[selected] + current_kmers[selected]];
		sk_counts[sk] = 0;
		result.emplace_back();
		vector<uint8_t *> & skmer = result.back();
		skmer.reserve(sk_size[sk]);
		// Jump over its kmers in the matrix
		size_t col = sk_col[sk];
		for (uint64_t id=sk_first[sk] ; id!=none ; id=next_kmer[id], col++) {
			skmer.push_back(matrix[col][current_kmers[col]]);
			current_kmers[col] += 1;
			remaining -= 1;
			if (current_kmers[col] < matrix[col].size())
				set_current(col);
		}
	}

	// free memory
	for (uint8_t * mem : memory)
		delete[](mem);

	return result;
}


vector<vector<uint8_t *> > Compact::sorted_assembly(vector<vector<uint8_t *> > & kmers) {
	// 1 - Sort Matrix per column
	this->sort_matrix(kmers);

    vector<vector<pair<uint64_t, uint64_t> > > kmer_pairs;

	// Pair columns
	for (uint i=0 ; i<this->k-this->m ; i++) {
		// 2 - Find all the possible overlaps of kmers
        vector<pair<uint64_t, uint64_t> > candidate_links = this->pair_kmers(kmers[i], kmers[i+1]);

		// 3 - Filter out kmer pairs that are not in optimal colinear chainings
		const vector<pair<uint64_t, uint64_t> > colinear_links = this->colinear_chaining(candidate_links);

		kmer_pairs.push_back(colinear_links);
	}

	// // 4 - Finish the ordering by sorting skmers that could have been interchanged
	const vector<vector<uint8_t *> > skmers = polish_sort(kmers, kmer_pairs);

	// return skmers;
	return skmers;
}


vector<pair<uint8_t *, uint8_t *> > Compact::greedy_assembly(vector<vector<uint8_t *> > & kmers) {
	uint nb_nucl = k - m;
	vector<pair<uint8_t *, uint8_t *> > assembly;

	// Index kmers from the 0th set
	for (uint8_t * kmer : kmers[0]) {
		assembly.emplace_back(nullptr, kmer);
	}

	// Kmers of the ith set sorted by suffix value. The sort is stable: kmers with the same
	// value stay in their order of the set (the candidates are tested in this order).
	struct Candidate {
		uint64_t value;
		uint8_t * kmer;
		bool used;
	};
	vector<Candidate> index;

	for (uint i=0 ; i<nb_nucl ; i++) {
		// Index kmers in ith set
		index.clear();
		for (uint8_t * kmer : kmers[i]) {
			// Get the suffix
			uint64_t val = subseq_to_uint(kmer, nb_nucl, 1, nb_nucl-1);
			index.push_back({val, kmer, false});
		}
		stable_sort(index.begin(), index.end(), [](const Candidate & a, const Candidate & b) {
			return a.value < b.value;
		});

		// link kmers from (i+1)th set to ith kmers.
		for (uint8_t * kmer : kmers[i+1]) {
			uint64_t val = subseq_to_uint(kmer, nb_nucl, 0, nb_nucl-2);
			auto candidate = lower_bound(index.begin(), index.end(), val, [](const Candidate & c, uint64_t v) {
				return c.value < v;
			});

			bool chaining_found = false;
			// verify complete matching for the remaining candidates kmers
			for ( ; candidate != index.end() and candidate->value == val ; candidate++) {
				if (candidate->used)
					continue;
				// If the kmers can be assembled
				if (sequence_compare(
							kmer, nb_nucl, 0, nb_nucl-2,
							candidate->kmer, nb_nucl, 1, nb_nucl-1
						) == 0) {
					// Update status
					chaining_found = true;
					assembly.emplace_back(candidate->kmer, kmer);
					// remove candidate from list
					candidate->used = true;
					// Quit candidate searching
					break;
				}
			}
			// If no assembly possible, create a new superkmer
			if (not chaining_found) {
				assembly.emplace_back(nullptr, kmer);
			}
		}
	}

	// Index last kmers without compaction
	int assembly_idx = assembly.size()-1;
	for (auto it=kmers[nb_nucl].end() ; it>kmers[nb_nucl].begin() ; it--) {
		uint8_t * kmer = *(it-1);

		if (assembly_idx < 0 or kmer != assembly[assembly_idx].second) {
			assembly.emplace_back(nullptr, kmer);
		} else {
			assembly_idx--;
		}
	}

	return assembly;
}

vector<vector<uint8_t *> > Compact::pairs_to_paths(const vector<pair<uint8_t *, uint8_t *> > & to_compact) {
	vector<vector<uint8_t *> > paths;
	// Path of the last kmer of each path, indexed by the position of the kmer in the kmer
	// buffer. 0 for the other kmers (the default value of the previous std::unordered_map).
	const uint64_t record_size = this->bytes_compacted + this->data_size + this->mini_pos_size;
	vector<uint> path_registry(this->next_free / record_size + 1, 0);
	auto registry = [&](uint8_t * kmer) -> uint & {
		return path_registry[(kmer - this->kmer_buffer) / record_size];
	};

	for (const pair<uint8_t *, uint8_t *> & p : to_compact) {
		// First element of a compaction path
		if (p.first == nullptr) {
			registry(p.second) = paths.size();
			paths.emplace_back(vector<uint8_t *>());

			uint vec_idx = registry(p.second);
			paths[vec_idx].push_back(p.second);
		}
		// Extending existing path
		else {
			uint vec_idx = registry(p.first);
			paths[vec_idx].push_back(p.second);
			registry(p.first) = 0;
			registry(p.second) = vec_idx;
		}
	}

	return paths;
}

uint64_t Compact::write_paths(const vector<vector<uint8_t *> > & paths, const MinimizerBlockEncoder & encoder, vector<uint8_t> & output) {
	const uint data_size = encoder.data_size;
	uint64_t nb_blocks = 0;
	uint kmer_bytes = (k - m + 3) / 4;
	uint kmer_offset = (4 - ((k - m) % 4)) % 4;
	// Width of the minimizer positions in the kmer buffer (computed from the input section)
	uint mini_pos_size = this->mini_pos_size;

	uint max_skmer_bytes = (2 * (k - m) + 3) / 4;
	uint8_t * skmer_buffer = new uint8_t[max_skmer_bytes + 1];
	uint data_bytes = (k - m + 1) * data_size;
	uint8_t * data_buffer = new uint8_t[data_bytes];

	// uint8_t encoding[] = {0, 1, 3, 2};
	// Stringifyer strif(encoding);

	// Write skmer per skmer
	for (const vector<uint8_t *> & path : paths) {
		// cout << "path " << ++idx << "/" << paths.size() << endl;
		// Cleaning previous skmers/data
		memset(skmer_buffer, 0, max_skmer_bytes + 1);
		memset(data_buffer, 0, data_bytes);

		// Get the skmer minimizer position
		uint mini_pos = 0;
		uint8_t * mini_pos_pointer = path[0] + kmer_bytes + data_size;
		for (uint b=0 ; b<mini_pos_size ; b++) {
			mini_pos <<= 8;
			mini_pos += *mini_pos_pointer;
			mini_pos_pointer += 1;
		}

		// Usefull variables
		uint skmer_size = k - m - 1 + path.size();
		// uint skmer_bytes = (skmer_size + 3) / 4;
		uint skmer_offset = (4 - (skmer_size % 4)) % 4;

		// cout << "first kmer" << endl;
		// Save the first kmer
		memcpy(skmer_buffer, path[0], kmer_bytes);
		leftshift8(skmer_buffer, kmer_bytes, 2 * kmer_offset);
		rightshift8(skmer_buffer, kmer_bytes + 1, 2 * skmer_offset);
		// Save the first data
		memcpy(data_buffer, path[0] + kmer_bytes, data_size);

		// Compact kmer+data one by one
		for (uint kmer_idx = 1 ; kmer_idx<path.size() ; kmer_idx++) {
			// cout << "kmer " << kmer_idx << endl;
			uint8_t * kmer = path[kmer_idx];

			// Compute compaction position
			uint compact_nucl_pos = skmer_offset + k - m - 1 + kmer_idx;
			uint compact_byte = compact_nucl_pos / 4;
			uint compact_shift = 3 - (compact_nucl_pos % 4);
			// Compact the nucleotide
			uint8_t last_nucl = kmer[kmer_bytes - 1] & 0b11;
			skmer_buffer[compact_byte] |= last_nucl << (2 * compact_shift);
			// Copy data
			memcpy(data_buffer + kmer_idx * data_size, path[kmer_idx] + kmer_bytes, data_size);
		}

		// cout << "write_compacted_sequence_without_mini" << endl;
		// Write everything in the section
		encoder.append_without_mini(output, skmer_buffer, skmer_size, mini_pos, data_buffer);
		nb_blocks += 1;
	}

	// cout << "delete" << endl;

	delete[] skmer_buffer;
	delete[] data_buffer;
	return nb_blocks;
}


