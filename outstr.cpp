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

void Outstr::exec() {
	// Read the encoding and prepare the translator
	Kff_reader reader = Kff_reader(input_filename);
	Stringifyer strif(reader.get_encoding());

	// Prepare sequence and data buffers
	uint8_t * nucleotides;
	uint8_t * data;

	while (reader.has_next()) {
		reader.next_kmer(&nucleotides, &data);
		cout << strif.translate(nucleotides, reader.get_var("k")) << " " << (uint)*data << endl;
	}
}
