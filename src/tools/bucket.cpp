#include <utility>
#include <vector>
#include <limits>
#include <string>
#include <cstring>
#include <utility>
#include <unordered_map>
#include <algorithm>
#include <queue>

#include "bucket.hpp"
#include "sequences.hpp"
#include "merge.hpp"

#include "encoding.hpp"

using namespace std;


Bucket::Bucket(uint8_t m, bool revcomp) {
	input_filename = "";
	output_filename = "";

	this->m = m;
	this->singleside = !revcomp;

	kmer_buffer = new uint8_t[1];
	data_buffer = new uint8_t[1];

	copy_buffer = new uint8_t[1];
}

Bucket::Bucket(std::string input_filename, std::string outpur_filename, uint8_t m, bool revcomp) {
    this->input_filename = std::move(input_filename);
    this->output_filename = std::move(outpur_filename);

    this->m = m;
    this->singleside = !revcomp;

    kmer_buffer = new uint8_t[1];
    data_buffer = new uint8_t[1];

    copy_buffer = new uint8_t[1];
}


Bucket::~Bucket() {
	delete[] kmer_buffer;
	delete[] data_buffer;
	delete[] copy_buffer;
}


void Bucket::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("bucket", "Read a kff file and split the kmers into buckets. Each bucket corresponds to all the kmers sharing the same minimizer. The minimizer of size m is the one that minimize the alphabetic order regarding the encoding. WARNING: If the minimzer is on the reverse strand of a kmer, the kmer will be reverse complemented in the output.");//" To avoid such thing, you can use the single-side flag.");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "Input kff file to bucket.");
	input_option->required();
	input_option->check(CLI::ExistingFile);
	CLI::Option * out_option = subapp->add_option("-o, --outfile", output_filename, "Kff to write (must be different from the input)");
	out_option->required();
	CLI::Option * mini_size = subapp->add_option("-m, --minimizer-size", m, "Minimizer size [Max 31].");
	mini_size->required();
	// subapp->add_flag("-s, --single-side", singleside, "Look for the minimizer only on the forward strand.");
}


namespace {

/* Asks the processor to load the cache line of address (no effect without GCC/Clang builtins) */
static inline void prefetch(const void * address) {
#if defined(__GNUC__) || defined(__clang__)
	__builtin_prefetch(address);
#else
	(void)address;
#endif
}

/* Minimizer -> index of its bucket, open addressing (linear probing), doubled when half full */
class MinimizerIndex {
public:
	static constexpr uint32_t none = numeric_limits<uint32_t>::max();

	MinimizerIndex() : keys(1 << 10), values(1 << 10, none), bits(10), count(0) {}

	uint32_t find(uint64_t key) const {
		for (uint64_t slot=this->slot_of(key) ; ; slot=(slot + 1) & (this->keys.size() - 1)) {
			if (this->values[slot] == none or this->keys[slot] == key)
				return this->values[slot];
		}
	}

	/** Loads the first slot of the key in the cache */
	void prefetch_slot(uint64_t key) const {
		const uint64_t slot = this->slot_of(key);
		prefetch(&this->keys[slot]);
		prefetch(&this->values[slot]);
	}

	void insert(uint64_t key, uint32_t value) {
		if (2 * (this->count + 1) > this->keys.size())
			this->grow();
		this->place(key, value);
		this->count += 1;
	}

private:
	vector<uint64_t> keys;
	vector<uint32_t> values;
	uint bits;
	size_t count;

	uint64_t slot_of(uint64_t key) const {
		return (key * 0x9E3779B97F4A7C15ULL) >> (64 - this->bits);
	}

	void place(uint64_t key, uint32_t value) {
		uint64_t slot = this->slot_of(key);
		while (this->values[slot] != none)
			slot = (slot + 1) & (this->keys.size() - 1);
		this->keys[slot] = key;
		this->values[slot] = value;
	}

	void grow() {
		vector<uint64_t> old_keys(this->keys.size() * 2);
		vector<uint32_t> old_values(this->values.size() * 2, none);
		old_keys.swap(this->keys);
		old_values.swap(this->values);
		this->bits += 1;
		for (size_t slot=0 ; slot<old_keys.size() ; slot++)
			if (old_values[slot] != none)
				this->place(old_keys[slot], old_values[slot]);
	}
};

/** Variables of the input file when a bucket is created: its section is written with them */
struct BucketVariables {
	uint64_t k;
	uint64_t max;
	uint64_t data_size;
};

/** A bucket, kept in memory until the end, then written as a minimizer section through
 * Section_Minimizer::write_compacted_sequence() */
struct BucketSection {
	// Index of the variables of the input file when the bucket was created
	uint32_t variables;
	// The minimizer (m nucleotides, right aligned)
	uint8_t minimizer[8];
	// The superkmers one after the other: number of nucleotides and minimizer position (varints),
	// nucleotides (2 bits each, right aligned), data of the kmers
	vector<uint8_t> superkmers;
};

/* Writes value as a varint at bytes (7 bits per byte, high bit set when bytes follow), returns the
 * address after it */
static inline uint8_t * put_varint(uint8_t * bytes, uint64_t value) {
	while (value >= 0x80) {
		*(bytes++) = static_cast<uint8_t>(value) | 0x80;
		value >>= 7;
	}
	*(bytes++) = static_cast<uint8_t>(value);
	return bytes;
}

/* Reads a varint written by put_varint and moves bytes after it */
static inline uint64_t get_varint(uint8_t * & bytes) {
	uint64_t value = 0;
	for (uint shift=0 ; ; shift += 7) {
		const uint8_t byte = *(bytes++);
		value |= static_cast<uint64_t>(byte & 0x7F) << shift;
		if (byte < 0x80)
			return value;
	}
}

}


void Bucket::exec() {
	uint nb_mutex = 64;

	// Open the sequence stream
	KffSeqStream stream(this->input_filename);
	// minimizer -> index of its bucket. The iteration order of these maps gives the order
	// of the sections in the output file.
	vector<unordered_map<uint64_t, uint32_t> > buckets;
	buckets.resize(nb_mutex);
	// Same index for the lookups (the maps are only used for the insertions, in the same order)
	MinimizerIndex bucket_index;
	vector<BucketSection> sections;
	// The configurations of variables of the input (usually only one)
	vector<BucketVariables> variable_sets;

	RevComp rc(stream.reader.get_encoding());

	{
	// Variables init by thread
	MinimizerSearcher * searcher = new MinimizerSearcher(0, m, stream.reader.get_encoding());
	vector<uint8_t> subseq(1);
	// Bucket of each skmer of the current sequence
	vector<uint32_t> skmer_buckets;
	uint max_seq = 1;
	uint8_t * seq = new uint8_t[(max_seq + 3) / 4];
	uint max_data = 1;
	uint8_t * data = new uint8_t[max_data];
	uint prev_k = 0;
	int64_t variables_idx = -1;
	uint64_t max = 0;
	uint cached_k = 0;
	uint cached_data_size = 0;
	uint64_t cached_reader_max = 0;

	// Read the stream line by line
	while(true) {
		int nb_kmers = 0;
		nb_kmers = stream.next_sequence(seq, max_seq, data, max_data);

		if (nb_kmers == 0)
			break;
		else if (nb_kmers < 0) {
			// Seq buffer update
			uint new_seq_size = stream.reader.k + stream.reader.max - 1;
			if (new_seq_size > max_seq) {
				delete[] seq;
				max_seq = new_seq_size;
				seq = new uint8_t[(max_seq + 3) / 4];
			}

			// Data buffer update
			uint new_data_size = stream.reader.max * stream.reader.data_size;
			if (new_data_size > max_data) {
				delete[] data;
				max_data = new_data_size;
				data = new uint8_t[max_data];
			}
			continue;
		}

		uint k = stream.reader.k;
		uint data_size = stream.reader.data_size;
		uint seq_size = k - 1 + nb_kmers;
		if (k != prev_k) {
			prev_k = k;
			delete searcher;
			searcher = new MinimizerSearcher(k, m, stream.reader.get_encoding());
		}
		// Room for the longest super-kmer (+1 byte needed by subsequence())
		if (subseq.size() < (seq_size + 3) / 4 + 1)
			subseq.resize((seq_size + 3) / 4 + 1);

		// Current variables (the variable max is read again only when the variables of the
		// reader change)
		if (k != cached_k or data_size != cached_data_size or stream.reader.max != cached_reader_max) {
			cached_k = k;
			cached_data_size = data_size;
			cached_reader_max = stream.reader.max;
			max = stream.reader.get_var("max");
		}
		if (variables_idx < 0 or variable_sets[variables_idx].k != k or variable_sets[variables_idx].max != max
				or variable_sets[variables_idx].data_size != data_size) {
			variables_idx = -1;
			for (size_t v=0 ; v<variable_sets.size() ; v++)
				if (variable_sets[v].k == k and variable_sets[v].max == max and variable_sets[v].data_size == data_size)
					variables_idx = v;
			if (variables_idx < 0) {
				variables_idx = variable_sets.size();
				variable_sets.push_back({k, max, data_size});
			}
		}

        // Skmer deduction
        vector<skmer> skmers = searcher->get_skmers(seq, seq_size);
		// 1 - Bucket of each skmer (new buckets created in the order of the skmers). With several
		// skmers, the slots of their minimizers in the index are loaded first.
		const bool prefetching = skmers.size() >= 4;
		if (prefetching)
			for (const skmer & sk : skmers)
				bucket_index.prefetch_slot(sk.minimizer);
		skmer_buckets.resize(skmers.size());
		for (size_t sk_idx=0 ; sk_idx<skmers.size() ; sk_idx++) {
			const skmer & sk = skmers[sk_idx];
			uint mutex_idx = sk.minimizer % nb_mutex;
			uint32_t bucket_idx = bucket_index.find(sk.minimizer);
			// New bucket
			if (bucket_idx == MinimizerIndex::none) {
				bucket_idx = sections.size();
				buckets[mutex_idx][sk.minimizer] = bucket_idx;
				bucket_index.insert(sk.minimizer, bucket_idx);
				sections.emplace_back();
				BucketSection & section = sections.back();
				section.variables = variables_idx;
				// The minimizer
				if (sk.minimizer_position < 0) {
					int mini_pos = - sk.minimizer_position - 1;
					subsequence(seq, seq_size, subseq.data(), mini_pos, mini_pos + m - 1);
					rc.rev_comp(subseq.data(), m);
				} else {
					subsequence(seq, seq_size, subseq.data(), sk.minimizer_position, sk.minimizer_position + m - 1);
				}
				memcpy(section.minimizer, subseq.data(), (m + 3) / 4);
			}
			skmer_buckets[sk_idx] = bucket_idx;
		}
		// The sections and the ends of their superkmers are loaded before the appends
		if (prefetching)
			for (uint32_t bucket_idx : skmer_buckets) {
				const BucketSection & section = sections[bucket_idx];
				prefetch(&section);
				prefetch(section.superkmers.data() + section.superkmers.size());
			}

		// 2 - Append the skmers to their buckets
		for (size_t sk_idx=0 ; sk_idx<skmers.size() ; sk_idx++) {
			const skmer & sk = skmers[sk_idx];
			const uint32_t bucket_idx = skmer_buckets[sk_idx];

			uint mini_pos = k + 2;
			// Get the subsequence
			subsequence(seq, seq_size, subseq.data(), sk.start_position, sk.stop_position);
			uint subseq_size = sk.stop_position - sk.start_position + 1;
			if (sk.minimizer_position >= 0) {
				mini_pos = sk.minimizer_position - sk.start_position;
			}
			// Get the rev subsequence
			else {
				rc.rev_comp(subseq.data(), subseq_size);
				mini_pos = sk.stop_position + sk.minimizer_position - m + 2;
				// Reverse data
				rc.rev_data(data + sk.start_position * data_size, data_size, subseq_size - k + 1);
			}

			// Save the skmer and the data of its kmers (as many bytes as its section will write)
			BucketSection & section = sections[bucket_idx];
			const BucketVariables & vars = variable_sets[section.variables];
			const uint64_t seq_bytes = (subseq_size + 3) / 4;
			const uint64_t data_bytes = (subseq_size - vars.k + 1) * vars.data_size;
			vector<uint8_t> & superkmers = section.superkmers;
			const size_t start = superkmers.size();
			superkmers.resize(start + 2 * 10 + seq_bytes + data_bytes);
			uint8_t * superkmer = put_varint(superkmers.data() + start, subseq_size);
			superkmer = put_varint(superkmer, mini_pos);
			memcpy(superkmer, subseq.data(), seq_bytes);
			superkmer += seq_bytes;
			if (data_bytes > 0)
				memcpy(superkmer, data + sk.start_position * data_size, data_bytes);
			superkmers.resize(superkmer + data_bytes - superkmers.data());
		}
	}

	delete searcher;
	delete[] seq;
	delete[] data;
	}

	// Write the buckets in one file. The file is the same as the merge of one file per
	// bucket: a variable section each time the variables change, then the minimizer sections.
	Kff_file outfile(output_filename, "w");
	outfile.set_indexation(true);
	outfile.write_encoding(stream.reader.get_encoding());
	std::string meta = "Merged file";
	outfile.write_metadata(meta.length(), (uint8_t *)meta.c_str());

	for (uint mutex_idx=0 ; mutex_idx<nb_mutex ; mutex_idx++)
		for (auto& it: buckets[mutex_idx]) {
			BucketSection & section = sections[it.second];
			const BucketVariables & vars = variable_sets[section.variables];

			// Variables of the section
			const pair<const char *, uint64_t> variables[] = {
				{"k", vars.k}, {"m", m}, {"max", vars.max}, {"data_size", vars.data_size}};
			bool v_section_needed = false;
			for (auto & p : variables)
				if (outfile.global_vars.find(p.first) == outfile.global_vars.end()
						or outfile.global_vars[p.first] != p.second)
					v_section_needed = true;
			if (v_section_needed) {
				Section_GV sgv(&outfile);
				for (auto & p : variables)
					sgv.write_var(p.first, p.second);
				sgv.close();
			}

			// Minimizer section: one block per superkmer
			Section_Minimizer sm(&outfile);
			sm.write_minimizer(section.minimizer);
			uint8_t * superkmer = section.superkmers.data();
			uint8_t * const end = superkmer + section.superkmers.size();
			while (superkmer < end) {
				const uint64_t seq_size = get_varint(superkmer);
				const uint64_t mini_pos = get_varint(superkmer);
				uint8_t * nucleotides = superkmer;
				superkmer += (seq_size + 3) / 4;
				uint8_t * kmer_data = superkmer;
				superkmer += (seq_size - vars.k + 1) * vars.data_size;
				sm.write_compacted_sequence(nucleotides, seq_size, mini_pos, kmer_data);
			}
			sm.close();
			// Free the memory as soon as the section is written
			vector<uint8_t>().swap(section.superkmers);
		}

	outfile.close();
}
