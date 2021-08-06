#include <vector>
#include <string>
#include <cstring>

#include "outstr.hpp"
#include "encoding.hpp"


using namespace std;


Outstr::Outstr() {
	input_filename = "";
}

void Outstr::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("outstr", "Output kmer sequences as string on stdout. One kmer per line is printed");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "The file to print");
	input_option->required();
	subapp->add_flag("-c, --reverse-complement", revcomp, "Print the minimal value between a kmer and its reverse complement");
}


string format_data(uint8_t * data, size_t data_size) {
	if (data_size == 0)
		return "";
	else if (data_size < 8) {
		uint val = data[0];
		for (uint i=1 ; i<data_size ; i++) {
			val <<= 8;
			val += data[i];
		}
		return std::to_string(val);
	} else {
		string val = "";
		for (uint i=0 ; i<data_size ; i++) {
			val += "[" + std::to_string((uint)data[i]) + "]";
		}
		return val;
	}
}



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
	uint8_t * rc_copy = new uint8_t[1];
	uint64_t k = 0;

	// Prepare sequence and data buffers
	uint8_t * nucleotides = nullptr;
	uint8_t * data = nullptr;

	while (reader.has_next()) {
		reader.next_kmer(nucleotides, data);

		if (not revcomp) {
			cout << strif.translate(nucleotides, reader.get_var("k")) << " ";
			cout << format_data(data, reader.get_var("data_size")) << endl;
		} else {
			// Change the size of rev comp datastruct if k changes
			if (reader.get_var("k") != k) {
				k = reader.get_var("k");
				delete[] rc_copy;
				rc_copy = new uint8_t[(k+3) / 4];
			}

			// Get the reverse complement
			memcpy(rc_copy, nucleotides, (k+3)/4);
			rc.rev_comp(rc_copy, k);

			if (inf_eq(nucleotides, rc_copy, k)) {
				cout << strif.translate(nucleotides, k) << " ";
				cout << format_data(data, reader.get_var("data_size")) << endl;
			} else {
				cout << strif.translate(rc_copy, k) << " ";
				cout << format_data(data, reader.get_var("data_size")) << endl;
			}
		}
	}
}
