#include <vector>
#include <string>

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
}


string format_data(uint8_t * data, size_t data_size) {
	if (data_size == 0)
		return "";
	else if (data_size < 8) {
		return to_string((uint)*data);
	} else {
		string val = "";
		for (uint i=0 ; i<data_size ; i++) {
			val += "[" + to_string(data[i]) + "]";
		}
		return val;
	}
}


void Outstr::exec() {
	// Read the encoding and prepare the translator
	Kff_reader reader = Kff_reader(input_filename);
	Stringifyer strif(reader.get_encoding());

	// Prepare sequence and data buffers
	uint8_t * nucleotides;
	uint8_t * data;

	while (reader.has_next()) {
		reader.next_kmer(&nucleotides, &data);
		cout << strif.translate(nucleotides, reader.get_var("k")) << " ";
		cout << format_data(data, reader.get_var("data_size")) << endl;
	}
}
