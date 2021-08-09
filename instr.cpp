#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <queue>
#include <cstring>
#include <cstdio>
#include <cstdlib>

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
	delimiter = ' ';
	data_delimiter = ',';
}

void Instr::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("instr", "Convert a text kmer file or a text sequence file into a kff file. Kmers or sequences must be 1 per line. If data size is more than 0, then the delimiters are used to split each line.");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "A text file with one sequence per line (sequence omitted if its size < k). Empty data is added (size defined by -d option).");
	input_option->required();
	CLI::Option * output_option = subapp->add_option("-o, --outfile", output_filename, "The kff output file name.");
	output_option->required();
	CLI::Option * k_opt = subapp->add_option("-k, --kmer-size", k, "Mandatory kmer size");
	k_opt->required();
	subapp->add_option("-d, --data-size", data_size, "Data size in Bytes (Default 0, max 8).");
	subapp->add_flag("-c, --counts", is_counts, "Tell instr to consider the input file as count list. One kmer and one count per line are expected (separeted by any char delimiter).");
	subapp->add_option("--delimiter", delimiter, "Character used as a delimiter between the sequence and the data (default ' ').");
	subapp->add_option("--data-delimiter", data_delimiter, "Character used as a delimiter between two kmers data from the same sequence (default ',').");
	subapp->add_option("-m, --max-kmer-seq", max_kmerseq, "The maximum number of kmer that can be inside of sequence in the output (default 255).");
}


/** Read txt sequence file (1 sequence per line)
 */
class TxtSeqStream : public SequenceStream {
private:
  std::fstream fs;
  Binarizer bz;

  uint buffer_size;
  uint8_t * seq_buffer;
  uint8_t * data_buffer;

  uint k;
  uint data_size;

  char delimiter;
  char data_delimiter;

public:
  TxtSeqStream(const std::string filename, const uint8_t encoding[4], uint k, uint data_size, char delim, char data_delim) 
      : fs(filename, std::fstream::in)
      , bz(encoding)
      , buffer_size(1024)
      , seq_buffer(new uint8_t[1024])
      , data_buffer(new uint8_t[(1024 * 4 - k + 1	) * data_size])
      , k(k)
      , data_size(data_size)
      , delimiter(delim)
      , data_delimiter(data_delim)
  {};
  ~TxtSeqStream() {
    this->fs.close();
    delete[] this->seq_buffer;
    delete[] this->data_buffer;
  }
  
  uint next_sequence(uint8_t * & seq, uint8_t * & data) {
		// Verify stream integrity
		if (not this->fs)
			return 0;

		// read next sequence
		string line;
		getline(this->fs, line);
		if (line.size() == 0)
			return 0;

		// Get the split limit
		size_t seq_size = line.find(this->delimiter);
		if (seq_size == string::npos)
			seq_size = line.size();

		// Update buffers
		if (seq_size > this->buffer_size * 4) {
			delete[] this->seq_buffer;
			this->buffer_size = (seq_size + 3) / 4;
			this->seq_buffer = new uint8_t[this->buffer_size];

			delete[] this->data_buffer;
			this->data_buffer = new uint8_t[(seq_size - k + 1) * this->data_size];
		}
		// convert/copy the sequence
		this->bz.translate(line, seq_size, this->seq_buffer);
		seq = this->seq_buffer;

		// If counts
		if (this->data_size > 0) {
			char * str_data = (char *)line.c_str() + seq_size + 1;

			unsigned long count = strtoul(str_data, &str_data, 10);
			for (uint i=0 ; i<this->data_size ; i++) {
				data[data_size-1-i] = (uint8_t)count & 0xFF;
				count >>= 8;
			}
		}

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

	// Open the input seq stream
	const uint8_t encoding[4] = {0, 1, 3, 2};
	TxtSeqStream stream(input_filename, encoding, this->k, this->data_size, this->delimiter, this->data_delimiter);

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
		// Full sequence copy
		if (nb_kmers <= this->max_kmerseq) {
			sr.write_compacted_sequence(seq, seq_size, data);
			continue;
		}

		// Sequence saved slice per slice
		uint first_nucl = 0;
		while (nb_kmers > 0) {
			uint nb_kmer_copied = min(nb_kmers, this->max_kmerseq);
			uint copy_size = nb_kmer_copied + (k - 1);
			
			uint last_nucl = first_nucl + (copy_size - 1);
			subsequence(seq, seq_size, sub_seq, first_nucl, last_nucl);
			first_nucl = last_nucl + 1 - (k - 1);

			// Write the sequence
			sr.write_compacted_sequence(sub_seq, copy_size, data);

			// reduce the number of remaining kmers
			nb_kmers -= nb_kmer_copied;
		}
	}
	sr.close();

	delete[] sub_seq;
	delete[] data;
	outfile.close();
}

