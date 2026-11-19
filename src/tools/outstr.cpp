#include <vector>
#include <string>
#include <charconv>
#include <cstdio>
#include <cstring>

#include "outstr.hpp"
#include "encoding.hpp"


using namespace std;


Outstr::Outstr() {
	input_filename = "";
	revcomp = false;
}

void Outstr::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("outstr", "Output kmer sequences as string on stdout. One kmer per line is printed");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "The file to print");
	input_option->required();
	input_option->check(CLI::ExistingFile);
	subapp->add_flag("-c, --reverse-complement", revcomp, "Print the minimal value between a kmer and its reverse complement");
}


char * outstr_write_data(const uint8_t * data, size_t data_size, char * out) {
	if (data_size == 0)
		return out;
	else if (data_size < 8) {
		// Value accumulated on an unsigned int (bytes beyond 32 bits are lost)
		uint val = data[0];
		for (uint i=1 ; i<data_size ; i++) {
			val <<= 8;
			val += data[i];
		}
		return std::to_chars(out, out + 10, val).ptr;
	} else {
		for (uint i=0 ; i<data_size ; i++) {
			*(out++) = '[';
			out = std::to_chars(out, out + 3, (uint)data[i]).ptr;
			*(out++) = ']';
		}
		return out;
	}
}


size_t outstr_data_chars(size_t data_size) {
	if (data_size == 0)
		return 0;
	else if (data_size < 8)
		return 10;
	else
		return 5 * data_size;
}


/** Characters written to stdout by large blocks */
class OutputBuffer {
private:
	std::vector<char> buffer;
	size_t used;

public:
	OutputBuffer() : buffer(1 << 20), used(0) {}
	~OutputBuffer() { this->flush(); }

	/** Space for at least size characters */
	char * reserve(size_t size) {
		if (this->used + size > this->buffer.size()) {
			this->flush();
			if (size > this->buffer.size())
				this->buffer.resize(size);
		}
		return this->buffer.data() + this->used;
	}
	/** Validate the characters written up to end */
	void commit(char * end) { this->used = end - this->buffer.data(); }
	void flush() {
		if (this->used > 0)
			fwrite(this->buffer.data(), 1, this->used, stdout);
		this->used = 0;
	}
};


bool inf_eq(uint8_t * seq1, uint8_t * seq2, uint64_t size) {
	uint64_t nb_bytes = (size + 3) / 4;

	// Test first byte
	uint16_t mask = (1 << ((size % 4) * 2)) - 1;
	uint8_t byte_seq1 = seq1[0] & mask;
	uint8_t byte_seq2 = seq2[0] & mask;
	if (byte_seq1 != byte_seq2)
		return byte_seq1 < byte_seq2;

	// Test all remaining bytes
	for (uint idx=1 ; idx<nb_bytes ; idx++)
		if (seq1[idx] != seq2[idx])
			return seq1[idx] < seq2[idx];

	// Equals
	return true;
}


void Outstr::exec() {
	// Read the encoding and prepare the translator
	Kff_reader reader = Kff_reader(input_filename);
	Stringifyer strif(reader.get_encoding());
	
	// Prepare revcomp
	RevComp rc(reader.get_encoding());
	vector<uint8_t> rc_copy;
	uint64_t k = 0;

	// Prepare sequence and data buffers
	uint8_t * nucleotides = nullptr;
	uint8_t * data = nullptr;

	// One line per kmer: "<kmer> <data>\n"
	OutputBuffer output;
	while (reader.next_kmer(nucleotides, data)) {
		const uint8_t * printed = nucleotides;

		if (revcomp) {
			// Change the size of rev comp datastruct if k changes
			if (reader.k != k) {
				k = reader.k;
				rc_copy.resize((k+3) / 4);
			}

			// Get the reverse complement
			memcpy(rc_copy.data(), nucleotides, (k+3)/4);
			rc.rev_comp(rc_copy.data(), k);

			if (not inf_eq(nucleotides, rc_copy.data(), k))
				printed = rc_copy.data();
		}

		char * line = output.reserve(reader.k + 2 + outstr_data_chars(reader.data_size));
		line = strif.write(printed, reader.k, line);
		*(line++) = ' ';
		line = outstr_write_data(data, reader.data_size, line);
		*(line++) = '\n';
		output.commit(line);
	}
}
