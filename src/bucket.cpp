#include <vector>
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

	uint prev_k = 0;
	uint8_t * subseq = nullptr;

	uint8_t * seq;
	uint8_t * data;
	uint nb_kmers;


	// Minimizers finding
	// Minimizer_Creator * mc = new Minimizer_Creator(0, 0, rc);
	MinimizerSearcher * ms = new MinimizerSearcher(0, m, stream.reader.get_encoding());

	while((nb_kmers = stream.next_sequence(seq, data)) != 0) {
		uint k = stream.reader.k;
		uint data_size = stream.reader.data_size;
		uint seq_size = k - 1 + nb_kmers;
		if (k != prev_k) {
			prev_k = k;
			delete[] subseq;
			subseq = new uint8_t[(k * 2 + 3) / 4];
			memset(subseq, 0, (k * 2 + 3) / 4);
			delete ms;
			ms = new MinimizerSearcher(k, m, stream.reader.get_encoding());
		}

		// Skmer deduction
		vector<skmer> skmers = ms->get_skmers(seq, seq_size);
		for (skmer sk : skmers) {
			// New bucket
			if (buckets.find(sk.minimizer) == buckets.end()) {
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
			}

			uint seq_size = k - 1 + nb_kmers;
			uint mini_pos = k + 2;
			// Get the subsequence
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

			// Save the skmer and its related data
			buckets[sk.minimizer]->write_compacted_sequence(
					subseq, subseq_size, mini_pos,
					data + sk.start_position * data_size
			);
		}
	}

	delete ms;

	// Close all the buckets
	for (auto& it: buckets) {
		Section_Minimizer * sm = it.second;
		Kff_file * outfile = sm->file;
		sm->close();
		delete sm;

		outfile->close(false);
	}

	delete[] subseq;

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