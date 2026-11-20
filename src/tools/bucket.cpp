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
#include "blocks.hpp"

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

/** The minimizer section of a bucket, kept in memory until the end */
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

struct BucketSection {
	// Encoder of the variables declared when the bucket was created (k, m, max, data_size)
	uint32_t encoder;
	// Bytes of the minimizer as written in the section header
	uint8_t minimizer[8];
	uint64_t nb_blocks;
	// Bytes of the blocks of the section
	vector<uint8_t> blocks;
};

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
	// One encoder per configuration of variables (usually only one)
	vector<MinimizerBlockEncoder> encoders;

	RevComp rc(stream.reader.get_encoding());

	{
	// Variables init by thread
	MinimizerSearcher * searcher = new MinimizerSearcher(0, m, stream.reader.get_encoding());
	vector<uint8_t> subseq(1);
	uint max_seq = 1;
	uint8_t * seq = new uint8_t[(max_seq + 3) / 4];
	uint max_data = 1;
	uint8_t * data = new uint8_t[max_data];
	uint prev_k = 0;
	int64_t encoder_idx = -1;
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

		// Encoder of the current variables (the variable max is read again only when the
		// variables of the reader change)
		if (k != cached_k or data_size != cached_data_size or stream.reader.max != cached_reader_max) {
			cached_k = k;
			cached_data_size = data_size;
			cached_reader_max = stream.reader.max;
			max = stream.reader.get_var("max");
		}
		if (encoder_idx < 0 or encoders[encoder_idx].k != k or encoders[encoder_idx].max != max
				or encoders[encoder_idx].data_size != data_size) {
			encoder_idx = -1;
			for (size_t e=0 ; e<encoders.size() ; e++)
				if (encoders[e].k == k and encoders[e].max == max and encoders[e].data_size == data_size)
					encoder_idx = e;
			if (encoder_idx < 0) {
				encoder_idx = encoders.size();
				encoders.emplace_back(k, m, max, data_size);
			}
		}

        // Skmer deduction
        vector<skmer> skmers = searcher->get_skmers(seq, seq_size);
        for (skmer sk : skmers) {
			uint mutex_idx = sk.minimizer % nb_mutex;
			uint32_t bucket_idx = bucket_index.find(sk.minimizer);
			// New bucket
			if (bucket_idx == MinimizerIndex::none) {
				bucket_idx = sections.size();
				buckets[mutex_idx][sk.minimizer] = bucket_idx;
				bucket_index.insert(sk.minimizer, bucket_idx);
				sections.emplace_back();
				BucketSection & section = sections.back();
				section.encoder = encoder_idx;
				section.nb_blocks = 0;
				// The minimizer
				if (sk.minimizer_position < 0) {
					int mini_pos = - sk.minimizer_position - 1;
					subsequence(seq, seq_size, subseq.data(), mini_pos, mini_pos + m - 1);
					rc.rev_comp(subseq.data(), m);
				} else {
					subsequence(seq, seq_size, subseq.data(), sk.minimizer_position, sk.minimizer_position + m - 1);
				}
				memcpy(section.minimizer, subseq.data(), encoders[encoder_idx].minimizer_bytes);
			}

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

			// Save the skmer and its related data (encoded with the variables of its section)
			BucketSection & section = sections[bucket_idx];
			encoders[section.encoder].append(section.blocks, subseq.data(), subseq_size, mini_pos,
			                                 data + sk.start_position * data_size);
			section.nb_blocks += 1;
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
			const MinimizerBlockEncoder & encoder = encoders[section.encoder];

			// Variables of the section
			const pair<const char *, uint64_t> variables[] = {
				{"k", encoder.k}, {"m", encoder.m}, {"max", encoder.max}, {"data_size", encoder.data_size}};
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

			// Minimizer section: type, minimizer, number of blocks, blocks
			outfile.register_position('m');
			uint8_t header[1 + 8 + 8];
			header[0] = 'm';
			memcpy(header + 1, section.minimizer, encoder.minimizer_bytes);
			for (uint b=0 ; b<8 ; b++)
				header[1 + encoder.minimizer_bytes + b] = static_cast<uint8_t>(section.nb_blocks >> (8 * (7 - b)));
			outfile.write(header, 1 + encoder.minimizer_bytes + 8);
			outfile.write(section.blocks.data(), section.blocks.size());
			// Free the memory as soon as the section is written
			vector<uint8_t>().swap(section.blocks);
		}

	outfile.close();
}
