#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <queue>
#include <limits>

#include "bucket.hpp"
#include "sequences.hpp"

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

	// TODO: merge
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
	Stringifyer strif(stream.reader.get_encoding());

	uint8_t * seq;
	uint8_t * data;
	uint nb_kmers;
	while((nb_kmers = stream.next_sequence(seq, data)) != 0) {
		uint k = stream.reader.get_var("k");
		cout << strif.translate(seq, nb_kmers + k - 1) << endl;
		vector<uint> candidates = compute_mini_candidates(seq, nb_kmers + k - 1, k, m);
		for (uint candidate : candidates) {
			uint8_t val = candidate;
			cout << "\t" << strif.translate(&val, 4) << endl;
		}
	}
}
