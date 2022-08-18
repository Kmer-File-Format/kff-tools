#include <utility>
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
	this->subapp = app->add_subcommand("bucket", "Read a kff file and split the kmers into buckets. Each bucket corresponds to all the kmers sharing the same minimizer. The minimizer of size m is the one that minimize the alphabetic order regarding the encoding. WARNING: If the minimzer is on the reverse strand of a kmer, the kmer will be reverse complemented in the output. To avoid such thing, you can use the single-side flag.");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "Input kff file to bucket.");
	input_option->required();
	input_option->check(CLI::ExistingFile);
	CLI::Option * out_option = subapp->add_option("-o, --outfile", output_filename, "Kff to write (must be different from the input)");
	out_option->required();
	CLI::Option * mini_size = subapp->add_option("-m, --minimizer-size", m, "Minimizer size [Max 31].");
	mini_size->required();
	subapp->add_flag("-s, --single-side", singleside, "Look for the minimizer only on the forward strand.");
}


void Bucket::exec() {
	uint nb_mutex = 64;

	// Open the sequence stream
	KffSeqStream stream(this->input_filename);
	vector<unordered_map<uint64_t, Section_Minimizer *> > buckets;
	buckets.resize(nb_mutex);

	RevComp rc(stream.reader.get_encoding());
	Stringifyer strif(stream.reader.get_encoding());

	// Prepare synchronisation locks
	vector<omp_lock_t> bucket_mutexes;
	bucket_mutexes.resize(nb_mutex);
	for (uint64_t i(0); i < nb_mutex; ++i) {
		omp_init_lock(&bucket_mutexes[i]);
	}

	// #pragma omp parallel num_threads(8)
	{
	// Variables init by thread
	MinimizerSearcher * searcher = new MinimizerSearcher(0, m, stream.reader.get_encoding());
	uint8_t * subseq = new uint8_t[1];
	uint max_seq = 1;
	uint8_t * seq = new uint8_t[(max_seq + 3) / 4];
	uint max_data = 1;
	uint8_t * data = new uint8_t[max_data];
	uint prev_k = 0;

	// Read the stream line by line
	while(true) {
		int nb_kmers = 0;
		#pragma omp critical
		{
			nb_kmers = stream.next_sequence(seq, max_seq, data, max_data);
		}

		if (nb_kmers == 0)
			break;
		else if (nb_kmers < 0) {
			// cout << "Realloc " << omp_get_thread_num() << endl;
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
			delete[] subseq;
			subseq = new uint8_t[(k * 2 + 3) / 4];
			memset(subseq, 0, (k * 2 + 3) / 4);
			delete searcher;
			searcher = new MinimizerSearcher(k, m, stream.reader.get_encoding());
		}

        // Skmer deduction
        vector<skmer> skmers = searcher->get_skmers(seq, seq_size);
        for (skmer sk : skmers) {
			uint mutex_idx = sk.minimizer % nb_mutex;
			omp_set_lock(&bucket_mutexes[mutex_idx]);
			// cout << "mutex " << mutex_id << endl;
			// New bucket
			if (buckets[mutex_idx].find(sk.minimizer) == buckets[mutex_idx].end()) {
				// cout << "New bucket " << sk.minimizer << endl;
				// Create the file
				Kff_file * outfile = new Kff_file(output_filename + "_" + to_string(sk.minimizer) + ".kff", "w");
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
			  buckets[mutex_idx][sk.minimizer] = sm;
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
				// cout << "No new bucket" << endl;
			}
			omp_unset_lock(&bucket_mutexes[mutex_idx]);

			uint seq_size = k - 1 + nb_kmers;
			uint mini_pos = k + 2;
			// Get the subsequence
			// cout << "Subsequence" << endl;
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

			// cout << "File write " << sk.minimizer << " " << sk.minimizer % 1024 << endl;
			// Save the skmer and its related data
			omp_set_lock(&bucket_mutexes[mutex_idx]);
			buckets[mutex_idx][sk.minimizer]->write_compacted_sequence(
					subseq, subseq_size, mini_pos,
					data + sk.start_position * data_size
			);

			// cout << "Xmutex " << sk.minimizer % 1024 << endl;
			omp_unset_lock(&bucket_mutexes[mutex_idx]);
		}
	}

	delete searcher;
	delete[] seq;
	delete[] data;
	delete[] subseq;
	}

	// Close all the buckets
	vector<Kff_file *> files;
	for (uint mutex_idx=0 ; mutex_idx<nb_mutex ; mutex_idx++)
		for (auto& it: buckets[mutex_idx]) {
			Section_Minimizer * sm = it.second;
			Kff_file * outfile = sm->file;
			sm->close();
			delete sm;

			files.push_back(outfile);
			outfile->close(false);
			outfile->open("r");
		}

	// Merge all the buckets in one file
	Merge mg;
	mg.merge(files, output_filename);

	for (Kff_file * file : files)
		delete file;
}

