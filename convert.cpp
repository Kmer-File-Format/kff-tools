#include <string>
#include <fstream>
#include <unordered_map>

#include "convert.hpp"
#include "encoding.hpp"


using namespace std;



Convert::Convert() {
	input_filename = "";
	output_filename = "";
	verbose = false;
}

void Convert::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("convert", "Convert a kmer text file to a kff file.");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "A list of kmers in text format (1 per line)");
	input_option->required();
	CLI::Option * output_option = subapp->add_option("-o, --outfile", output_filename, "The kff output file");
	output_option->required();
}

void Convert::exec() {
		// prepare outfile
		Kff_file outfile(output_filename, "w");
		outfile.write_encoding(0, 1, 3, 2);

		ifstream txt_file(input_filename);
		string line;
		getline(txt_file, line);
		txt_file.close();

		uint64_t k = line.substr(0, line.find(" ")).length();

		Section_GV sgv(&outfile);
		sgv.write_var("k", k);
		sgv.write_var("max", 1);
		sgv.write_var("data_size", 0);
		sgv.close();
		
		Section_Raw sr(&outfile);

		txt_file = ifstream(input_filename);
		Binarizer brz(outfile.encoding);
		uint8_t * bin = new uint8_t[k / 4 + 1];
		while (getline(txt_file, line)) {
			string kmer = line.substr(0, line.find(" "));

			// Change the size of k
			if (kmer.length() != k) {
				sr.close();
				k = kmer.length();
				Section_GV sgv(&outfile);
				sgv.write_var("k", k);
				sgv.close();
				sr = Section_Raw(&outfile);
			}

			brz.translate(kmer, bin);
			sr.write_compacted_sequence(bin, k, nullptr);
		}

		delete[] bin;
		sr.close();
		outfile.close();
}
