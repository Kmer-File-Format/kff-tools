#include <vector>
#include <string>
#include <cstring>
#include <utility>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include "omp.h"

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

	RevComp rc(stream.reader.get_encoding());
	Stringifyer strif(stream.reader.get_encoding());

	// Prepare synchronisation locks
	vector<omp_lock_t> bucket_mutexes;
	bucket_mutexes.resize(1024);
	for (uint64_t i(0); i < 1024; ++i) {
		omp_init_lock(&bucket_mutexes[i]);
	}

	#pragma omp parallel num_threads(2)
	{
	// Variables init by thread
	MinimizerSearcher * searcher = new MinimizerSearcher(0, m, stream.reader.get_encoding());
	uint8_t * subseq = new uint8_t[1];
	uint max_seq = 1;
	uint8_t * seq = new uint8_t[(max_seq + 3) / 4];
	uint max_data = 1;
	uint8_t * data = new uint8_t[max_data];
	uint prev_k = 0;
	cout << "data " << (int *)data << " " << omp_get_thread_num() << endl;

	// Read the stream line by line
	while(true) {
		uint nb_kmers;
		#pragma omp critical
		{
			cout << "Start " << omp_get_thread_num() << endl;
			while ((nb_kmers = stream.next_sequence(seq, max_seq, data, max_data)) < 0) {
				cout << "Realloc " << omp_get_thread_num() << endl;
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
			}
			cout << nb_kmers << endl;
			cout << "Stop critic " << omp_get_thread_num() << endl;
		}
		if (nb_kmers == 0)
			break;

		uint k = stream.reader.k;
		uint data_size = stream.reader.data_size;
		uint seq_size = k - 1 + nb_kmers;
		if (k != prev_k) {
			prev_k = k;
			delete[] subseq;
			subseq = new uint8_t[(k * 2 + 3) / 4];
			memset(subseq, 0, (k * 2 + 3) / 4);
			delete searcher;
			searcher = new MinimizerSearcher(k, m, stream.reader.get_encoding());
		}

		// Skmer deduction
		vector<skmer> skmers = searcher->get_skmers(seq, seq_size);
		for (skmer sk : skmers) {
			omp_set_lock(&bucket_mutexes[sk.minimizer % 1024]);
			cout << "mutex " << sk.minimizer % 1024 << endl;
			// New bucket
			if (buckets.find(sk.minimizer) == buckets.end()) {
				cout << "New bucket " << sk.minimizer << endl;
				// Create the file
				Kff_file * outfile = new Kff_file(output_filename + "_" + to_string(sk.minimizer) + ".kff", "w");
				opened_files[sk.minimizer] = outfile;
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
			  buckets[sk.minimizer] = sm;
			  // Write the minimizer
			  if (sk.minimizer_position < 0) {
			  	int mini_pos = - sk.minimizer_position - 1;
			  	subsequence(seq, seq_size, subseq, mini_pos, mini_pos + m - 1);
			  	rc.rev_comp(subseq, m);
			  } else {
			  	subsequence(seq, seq_size, subseq, sk.minimizer_position, sk.minimizer_position + m - 1);
			  }
			  sm->write_minimizer(subseq);
			} else {
				cout << "No new bucket" << endl;
			}

			uint seq_size = k - 1 + nb_kmers;
			uint mini_pos = k + 2;
			// Get the subsequence
			cout << "Subsequence" << endl;
			subsequence(seq, seq_size, subseq, sk.start_position, sk.stop_position);
			uint subseq_size = sk.stop_position - sk.start_position + 1;
			if (sk.minimizer_position >= 0) {
				mini_pos = sk.minimizer_position - sk.start_position;
			}
			// Get the rev subsequence
			else {
				rc.rev_comp(subseq, subseq_size);
				mini_pos = sk.stop_position + sk.minimizer_position - m + 2;
				// Reverse data
				rc.rev_data(data + sk.start_position * data_size, data_size, subseq_size - k + 1);
			}

			cout << "File write" << endl;
			// Save the skmer and its related data
			buckets[sk.minimizer]->write_compacted_sequence(
					subseq, subseq_size, mini_pos,
					data + sk.start_position * data_size
			);

			cout << "Xmutex " << sk.minimizer % 1024 << endl;
			omp_unset_lock(&bucket_mutexes[sk.minimizer % 1024]);
		}
	}

	delete searcher;
	delete[] seq;
	delete[] data;
	delete[] subseq;
	}

	// Close all the buckets
	for (auto& it: buckets) {
		Section_Minimizer * sm = it.second;
		Kff_file * outfile = sm->file;
		sm->close();
		delete sm;

		outfile->close(false);
	}


	// Prepare bucket merging
	vector<Kff_file *> files;
	files.reserve(opened_files.size());
	for(auto & kv : opened_files) {
		// cout << kv.second->filename << endl;
		kv.second->open("r");
    files.push_back(kv.second);
	}

	// Merge all the buckets in one file
	Merge mg;
	mg.merge(files, output_filename);

	for (Kff_file * file : files)
		delete file;
}