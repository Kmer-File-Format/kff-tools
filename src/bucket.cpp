#include <vector>
#include <string>
#include <cstring>
#include <utility>
#include <unordered_map>
#include <algorithm>
#include <queue>

#include <sys/resource.h>

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


	unordered_map<uint64_t, Kff_file *> opened_files;
	unordered_map<uint64_t, uint64_t> opened_file_counts;
	queue<uint64_t> mini_fifo;
	struct rlimit nb_file_descriptors;
	getrlimit(RLIMIT_NOFILE, &nb_file_descriptors);

	const uint opened_limit = nb_file_descriptors.rlim_cur - 50;

	RevComp rc(stream.reader.get_encoding());
	uint max_seq_size = stream.reader.k - 1 +  stream.reader.get_var("max");
	MinimizerSearcher ms(stream.reader.k, this->m, max_seq_size, this->singleside, stream.reader.get_encoding());

	uint prev_k = 0;
	uint8_t * subseq = nullptr;

	uint8_t * seq;
	uint8_t * data;

	uint8_t * rev_comp = new uint8_t[(max_seq_size + 3) / 4];
	
	uint nb_kmers;
	while((nb_kmers = stream.next_sequence(seq, data)) != 0) {
		uint k = stream.reader.k;
		uint data_size = stream.reader.data_size;
		uint seq_size = k - 1 + nb_kmers;
		if (k != prev_k) {
			prev_k = k;
			delete[] subseq;
			subseq = new uint8_t[(k * 2 + 3) / 4];
			memset(subseq, 0, (k * 2 + 3) / 4);
		}

		// Reverse complement
		if (not this->singleside) {
			memcpy(rev_comp, seq, (seq_size + 3) / 4);
			rc.rev_comp(rev_comp, seq_size);
		}

		// TODO : What to do if k of max change ?

		// Get superkmers
		const vector<skmer> & skmers = ms.get_skmers(seq, seq_size);
		// cout << "NEW SEQUENCE" << endl;

		// Skmer deduction
		for (const skmer & skmer : skmers) {
			uint64_t mini_val = skmer.minimizer;

			// Limit the number of simultaneously opened bucket files
			while (opened_file_counts.size() == opened_limit) {
				uint64_t front_mini = mini_fifo.front();
				mini_fifo.pop();

				if (opened_file_counts[front_mini] == 1) {
					opened_file_counts.erase(front_mini);
					opened_files[front_mini]->tmp_close();
				} else {
					opened_file_counts[front_mini] -= 1;
				}
			}
			
			// Save the minimizer usage
			if (mini_fifo.size() == 0 or mini_fifo.back() != mini_val) {
				mini_fifo.push(mini_val);
				// Register file usage
				if (opened_file_counts.find(mini_val) == opened_file_counts.end()) {
					opened_file_counts[mini_val] = 0;
				}
				opened_file_counts[mini_val] += 1;
			}

			// New bucket
			if (buckets.find(mini_val) == buckets.end()) {
				// Create the file
				// cout << output_filename + "_" + to_string(mini_val) + ".kff" << endl;
				Kff_file * outfile = new Kff_file(output_filename + "_" + to_string(mini_val) + ".kff", "w");
				opened_files[mini_val] = outfile;
				outfile->write_encoding(stream.reader.get_encoding());
				outfile->set_uniqueness(stream.reader.file->uniqueness);
				outfile->set_canonicity(stream.reader.file->canonicity);

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
			  if (skmer.minimizer_position < 0) {
			  	int abs_mini_pos = -skmer.minimizer_position - 1;
			  	uint start = RevComp::rev_position(abs_mini_pos + m - 1, seq_size);
			  	uint stop = RevComp::rev_position(abs_mini_pos, seq_size);
			  	subsequence(rev_comp, seq_size, subseq, start, stop);
			  } else {
			  	subsequence(seq, seq_size, subseq, skmer.minimizer_position, skmer.minimizer_position + m - 1);
			  }
			  sm->write_minimizer(subseq);
			}

			uint seq_size = k - 1 + nb_kmers;
			uint mini_pos = k + 2;
			// Get the fwd subsequence
			if (skmer.minimizer_position >= 0) {
				subsequence(seq, seq_size, subseq, skmer.start_position, skmer.stop_position);
				mini_pos = skmer.minimizer_position - skmer.start_position;
			}
			// Get the rev subsequence
			else {
				uint start = RevComp::rev_position(skmer.stop_position, seq_size);
			  uint stop = RevComp::rev_position(skmer.start_position, seq_size);
				subsequence(rev_comp, seq_size, subseq, start, stop);
				// cout << "reverse " << skmer.minimizer_position << " " << start << " " << stop << endl;
				mini_pos = RevComp::rev_position(-skmer.minimizer_position - 1 + m - 1, seq_size) - start;
				// cout << mini_pos << endl;
				uint sub_nb_kmer = skmer.stop_position - skmer.start_position + 1 - k;
				rc.rev_data(data + skmer.start_position * data_size, data_size, sub_nb_kmer);
			}

			buckets[mini_val]->write_compacted_sequence(
					subseq, skmer.stop_position - skmer.start_position + 1, mini_pos,
					data + skmer.start_position * data_size
			);
		}
	}

	// Close all the buckets
	// uint i=1;
	for (auto& it: buckets) {
		Section_Minimizer * sm = it.second;
		sm->close();
		delete sm;
	}

	for (auto& it: opened_files) {
		Kff_file * outfile = it.second;
		outfile->close();
		delete outfile;
	}

	delete[] subseq;
	delete[] rev_comp;

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
