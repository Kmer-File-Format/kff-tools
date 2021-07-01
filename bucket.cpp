#include <vector>
#include <string>
#include <cstring>
#include <utility>
#include <unordered_map>
#include <algorithm>

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


Bucket::~Bucket() {
	delete[] kmer_buffer;
	delete[] data_buffer;
	delete[] copy_buffer;
}


void Bucket::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("bucket", "Read a kff file and split the kmers into buckets. Each bucket corresponds to all the kmers sharing the same minimizer. The minimizer of size m is the one that minimize the alphabetic order regarding the encoding. WARNING: If the minimzer is on the reverse strand of a kmer, the kmer will be reverse complemented in the output. To avoid such thing, you can use the single-side flag.");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "Input kff file to bucket.");
	input_option->required();
	CLI::Option * out_option = subapp->add_option("-o, --outfile", output_filename, "Kff to write (must be different from the input)");
	out_option->required();
	CLI::Option * mini_size = subapp->add_option("-m, --minimizer-size", m, "Minimizer size [Max 31].");
	mini_size->required();
	subapp->add_flag("-s, --single-side", singleside, "Look for the minimizer only on the forward strand.");
}

void Bucket::exec() {
	// Open the sequence stream
	KffSeqStream stream(this->input_filename);
	unordered_map<uint64_t, Section_Minimizer *> buckets;
	RevComp rc(stream.reader.get_encoding());
	Stringifyer strif(stream.reader.get_encoding());

	uint prev_k = 0;
	uint8_t * subseq = nullptr;

	uint8_t * seq;
	uint8_t * data;
	uint nb_kmers;
	while((nb_kmers = stream.next_sequence(seq, data)) != 0) {
		uint k = stream.reader.get_var("k");
		uint data_size = stream.reader.get_var("data_size");
		uint seq_size = k - 1 + nb_kmers;
		if (k != prev_k) {
			prev_k = k;
			delete[] subseq;
			subseq = new uint8_t[(k * 2 + 3) / 4];
			memset(subseq, 0, (k * 2 + 3) / 4);
		}

		// Minimizers finding
		vector<pair<int, uint64_t> > minimizers = compute_minizers(seq, nb_kmers + k - 1, k, m, rc, singleside);

		// Skmer deduction
		// vector<pair<uint, uint> > skmers;
		vector<pair<int, int> > skmers = compute_skmers(nb_kmers + k - 1, k, m, minimizers);
		for (uint i=0 ; i<minimizers.size() ; i++) {
			pair<int, int> & skmer_boundaries = skmers[i];
			pair<int, uint64_t> & minimizer = minimizers[i];

			uint64_t mini_val = minimizer.second;
			// New bucket
			if (buckets.find(mini_val) == buckets.end()) {
				// Create the file
				Kff_file * outfile = new Kff_file(output_filename + "_" + to_string(mini_val) + ".kff", "w");
				outfile->write_encoding(stream.reader.get_encoding());
				// Usefull values
				Section_GV sgv(outfile);
			  sgv.write_var("k", k);
			  sgv.write_var("m", m);
			  sgv.write_var("max", stream.reader.get_var("max"));
			  sgv.write_var("data_size", data_size);
			  sgv.close();
			  // Create the bucket by itself
			  Section_Minimizer * sm = new Section_Minimizer(outfile);
			  buckets[mini_val] = sm;
			  // Write the minimizer
			  if (minimizer.first < 0) {
			  	int mini_pos = seq_size + minimizer.first - m + 1;
			  	subsequence(seq, seq_size, subseq, mini_pos, mini_pos + m - 1);
			  	rc.rev_comp(subseq, m);
			  } else {
			  	subsequence(seq, seq_size, subseq, minimizer.first, minimizer.first + m - 1);
			  }
			  sm->write_minimizer(subseq);
			}

			uint seq_size = k - 1 + nb_kmers;
			uint subseq_size;
			uint mini_pos = k + 2;
			// Get the fwd subsequence
			if (skmer_boundaries.first >= 0) {
				subseq_size = skmer_boundaries.second - skmer_boundaries.first + 1;
				subsequence(seq, seq_size, subseq, skmer_boundaries.first, skmer_boundaries.second);
				mini_pos = minimizer.first - skmer_boundaries.first;
			}
			// Get the rev subsequence
			else {
				subseq_size = skmer_boundaries.first - skmer_boundaries.second + 1;
				subsequence(seq, seq_size, subseq, seq_size + skmer_boundaries.second, seq_size + skmer_boundaries.first);
				rc.rev_comp(subseq, subseq_size);
				mini_pos = - (minimizer.first - skmer_boundaries.first);
			}

			// Save the skmer and its related data
			buckets[mini_val]->write_compacted_sequence(
					subseq, subseq_size, mini_pos,
					data + skmer_boundaries.first * data_size
			);
		}
	}

	// Close all the buckets
	for (auto& it: buckets) {
		Section_Minimizer * sm = it.second;
		Kff_file * outfile = sm->file;
		sm->close();
		delete sm;

		outfile->close();
		delete outfile;
	}

	delete[] subseq;

	// Prepare bucket merging
	vector<string> bucket_names;
	for (const pair<uint64_t, Section_Minimizer *> & bucket : buckets)
		bucket_names.push_back(output_filename + "_" + to_string(bucket.first) + ".kff");
	sort(bucket_names.begin(), bucket_names.end());

	// Merge all the buckets in one file
	Merge mg;
	mg.merge(bucket_names, output_filename);

	for (string & filename : bucket_names) {
		remove(filename.c_str());
	}
}
