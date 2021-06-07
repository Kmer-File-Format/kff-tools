#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <queue>
#include <limits>

#include "bucket.hpp"
#include "encoding.hpp"
#include "sequences.hpp"


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
	// Read the encoding and prepare the translator
	Kff_reader reader = Kff_reader(input_filename);

	// Prepare complement values
	uint8_t * encoding = reader.get_encoding();
	for (uint i=0 ; i<4 ; i++)
		this->complement[encoding[i]] = encoding[3-i];

	// Prepare sequence and data buffers
	uint8_t * nucleotides = nullptr;
	uint8_t * data = nullptr;

	Stringifyer strif(reader.get_encoding());

	// Read the kmers one at a time
	while (reader.has_next()) {
		// Get the kmer
		reader.next_kmer(nucleotides, data);

		// Controle the variables
		uint64_t k = reader.get_var("k");
		if (m > k) {
			cerr << "m larger than k, kmer ignored " << strif.translate(nucleotides, k);
			continue;
		}

		// Bucketize
		this->bucketize(nucleotides, data, k);
		// return;
		// Merge is done on the object destruction
	}
}

void Bucket::bucketize(uint8_t * kmer, uint8_t * data, uint64_t k) {
	// Compute the minimizer
	int mini_pos;
	uint mini = get_minimizer(kmer, k, mini_pos);
	cout << mini_pos << " " << mini << endl;

	// Save to the right file
}

/** Compute the encoding dependant minimizer. Will return the minimizer regarding the encoding.
  *
  * @param kmer The 2-bits encoded byte sequence for the kmer.
  * @param k kmer length
  * @param position The variable where the minimizer position will be registered.
  *
  * @return The minimizer value.
  **/
uint Bucket::get_minimizer(const uint8_t * kmer, const uint64_t k, int & position) const {
	uint minimizer_value = std::numeric_limits<uint>::max();

	uint current_minimizer = 0;
	uint current_rev_minimzer = 0;
	uint kmer_nb_bytes = (k+3)/4;

	// Prepare for the first minimizer
	for (uint r_idx=0; r_idx<m-1 ; r_idx++) {
		// Get the nucleotide
		uint byte_idx = (kmer_nb_bytes - 1) - (r_idx / 4);
		uint nucl_offset = r_idx % 4;
		uint nucleotide = (kmer[byte_idx] >> (2 * nucl_offset)) & 0b11;

		// Set the nucleotide to the minimizer candidate
		current_minimizer >>= 2;
		current_minimizer += nucleotide << (2 * (m-1));

		// Set the nucleotide to the rev minizer candidate
		current_rev_minimzer <<= 2;
		current_rev_minimzer += this->complement[nucleotide];
	}

	// Compare all the minimizers
	uint rev_mask = (1 << (2 * m)) - 1;
	for (uint r_idx=m-1 ; r_idx<k ; r_idx++) {
		// Get the nucleotide
		uint byte_idx = (kmer_nb_bytes - 1) - (r_idx / 4);
		uint nucl_offset = r_idx % 4;
		uint nucleotide = (kmer[byte_idx] >> (2 * nucl_offset)) & 0b11;

		// Set the nucleotide to the minimizer candidate
		current_minimizer >>= 2;
		current_minimizer += nucleotide << (2 * (m-1));

		// Set the nucleotide to the rev minizer candidate
		current_rev_minimzer <<= 2;
		current_rev_minimzer += this->complement[nucleotide];
		current_rev_minimzer &= rev_mask;

		if (current_minimizer <= minimizer_value) {
			minimizer_value = current_minimizer;
			position = k - m - (r_idx - m + 1);
			// cout << r_idx << " " << position << "    ";
		} else if (!singleside and current_rev_minimzer < minimizer_value) {
			minimizer_value = current_rev_minimzer;
			position = - (r_idx - m + 1);
			// cout << "HOLA" << endl;
		}
	}
	cout << endl;

	return minimizer_value;
}
