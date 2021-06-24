#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <queue>
#include <cstring>
#include <cstdio>

#include "instr.hpp"
#include "encoding.hpp"
#include "merge.hpp"
#include "sequences.hpp"


using namespace std;



Instr::Instr() {
	input_filename = "";
	output_filename = "";
	data_size = 0;
	is_counts = false;
	k = 0;
	max_kmerseq = 255;
}

void Instr::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("instr", "Convert a text kmer file into one or multuple kff file(s). Instr suppose that data are signed integers (max 64 bits).");
	CLI::Option * k_opt = subapp->add_option("-k, --kmer-size", k, "Mandatory kmer size");
	k_opt->required();
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "A text file with one sequence per line (sequence omitted if its size < k). Empty data is added (size defined by -d option).");
	input_option->required();
	subapp->add_flag("-c, --counts", is_counts, "Tell instr to consider the input file as count list. One kmer and one count per line are expected (separeted by any char delimiter).");
	CLI::Option * output_option = subapp->add_option("-o, --outprefix", output_filename, "The kff output file prefix. For a single file the name will be <prefix>.kff, otherwise, one file per minimizer is created (<prefix>_<minimizer>.kff).");
	output_option->required();
	subapp->add_option("-d, --data-size", data_size, "Data size in Bytes (Default 0, max 8).");
	subapp->add_option("-m, --max-kmer-seq", max_kmerseq, "The maximum number of kmer that can be inside of sequence in the output (default 255).");
}


/** Read txt sequence file (1 sequence per line)
 */
class TxtSeqStream : public SequenceStream {
private:
  std::fstream fs;
  Binarizer bz;

  uint buffer_size;
  uint8_t * buffer;

  uint k;
  uint data_size;

public:
  TxtSeqStream(const std::string filename, const uint8_t encoding[4]) 
      : fs(filename, std::fstream::in)
      , bz(encoding)
      , buffer_size(1024)
      , buffer(new uint8_t[1024])
      , k(0)
      , data_size(0)
  {};
  ~TxtSeqStream() {
    this->fs.close();
    delete[] buffer;
  }

  void set_counts(uint k, uint data_size) {
  	this->k = k;
  	this->data_size = data_size;
  }
  
  uint next_sequence(uint8_t * & seq, uint8_t * & data) {
		// Verify stream integrity
		if (not this->fs)
			return 0;

		// read next sequence
		string line;
		getline(this->fs, line);
		cout << line << endl;
		// Update buffer
		uint seq_size = line.size();
		if (seq_size > this->buffer_size * 4) {
			delete[] this->buffer;
			this->buffer_size = (seq_size + 3) / 4;
			this->buffer = new uint8_t[this->buffer_size];
		}
		// convert/copy the sequence
		this->bz.translate(line, this->buffer);
		seq = this->buffer;

		// If counts
		if (this->k > 0) {
			string str_count = line.substr(k+1);
			unsigned long count = stoul(str_count);
			for (uint i=0 ; i<this->data_size ; i++) {
				data[data_size-1-i] = (uint8_t)count & 0xFF;
				count >>= 8;
			}
		}

		cout << "Seq size " << seq_size << endl << endl;
		return seq_size;
	}
};


void Instr::exec() {
	// reset data size to 0 if data are not counts
	if (this->is_counts) {
		this->max_kmerseq = 1;
	} else
		this->data_size = 0;

	// Open a KFF for output
	Kff_file outfile(this->output_filename, "w");
	// Write needed variables
	Section_GV sgv(&outfile);
	sgv.write_var("k", this->k);
	sgv.write_var("data_size", this->data_size);
	sgv.write_var("ordered", 0);
	sgv.write_var("max", this->max_kmerseq);
	sgv.close();

	cout << "k " << k << endl;
	cout << "max " << max_kmerseq << endl;

	// Open the input seq stream
	const uint8_t encoding[4] = {0, 1, 3, 2};
	TxtSeqStream stream(input_filename, encoding);
	if (this->is_counts)
		stream.set_counts(this->k, this->data_size);

	// Write the sequences inside of a raw section
	Section_Raw sr(&outfile);

	uint8_t * seq;
	uint8_t * sub_seq = new uint8_t[(max_kmerseq + 3) / 4];
	uint8_t * data = new uint8_t[data_size];
	uint seq_size = 0;
	// read the next line from the txt file
	while ((seq_size = stream.next_sequence(seq, data)) > 0) {
		// Sequence too small
		if (seq_size < k)
			continue;

		uint nb_kmers = seq_size - this->k + 1;
		uint first_nucl = 0;
		while (nb_kmers > 0) {
			uint8_t * to_copy;
			uint nb_kmer_copied = min(nb_kmers, this->max_kmerseq);
			uint copy_size = nb_kmer_copied + (k - 1);
			cout << "copy size " << copy_size << endl;
			// If sequence < max
			if (nb_kmers < this->max_kmerseq)
				to_copy = seq;
			// Get subsequence if needed
			else {
				uint last_nucl = first_nucl + (copy_size - 1);
				subsequence(seq, seq_size, sub_seq, first_nucl, last_nucl);
				first_nucl = last_nucl + 1;
				to_copy = sub_seq;
			}

			// Write the sequence
			sr.write_compacted_sequence(to_copy, copy_size, data);

			// reduce the number of remaining kmers
			nb_kmers -= nb_kmer_copied;
		}
	}
	sr.close();

	delete[] sub_seq;
	delete[] data;
	outfile.close();
}

