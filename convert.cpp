#include <string>
#include <fstream>
#include <unordered_map>

#include "convert.hpp"
#include "encoding.hpp"


using namespace std;



Convert::Convert() {
	input_filename = "";
	output_prefix = "";
	m = 0;
	verbose = false;
}

void Convert::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("convert", "Convert a kmer text file to a kff file.");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "A list of kmers in text format (1 per line)");
	input_option->required();
	CLI::Option * output_option = subapp->add_option("-o, --outprefix", output_prefix, "The kff output file prefix. If m is not set, one file <prefix>.kff is output, otherwise, one file per minimizer is created (<prefix>_<minimizer>.kff)");
	output_option->required();
	subapp->add_option("-m, --mini-size", m, "Set a minimizer size. If set, all the kmers are separated into bins, one bin per minimizer.");
}

void Convert::exec() {
	if (m == 0) {
		this->monofile();
	} else {
		this->multifile();
	}
}

void Convert::multifile() {
	map<uint, Kff_file *> outfiles;
	map<uint, Section_Minimizer *> outsections;

	// Prepare first k value
	ifstream txt_file(input_filename);
	string line;
	getline(txt_file, line);
	txt_file.close();

	uint64_t k = line.substr(0, line.find(" ")).length();

	// Prepare the binarized datastructures
	uint8_t encoding[4] = {0, 1, 3, 2};
	Binarizer brz(encoding);
	uint8_t * bin = new uint8_t[k / 4 + 1];
	// Read the input kmer stream
	txt_file = ifstream(input_filename);
	while (getline(txt_file, line)) {
		string kmer = line.substr(0, line.find(" "));

		// Change the size of k
		if (kmer.length() != k) {
			k = kmer.length();
			delete[] bin;
			bin = new uint8_t[k / 4 + 1];
		}

		// Binarize the kmer
		brz.translate(kmer, bin);
		// Search for minimizer
		uint minimizer = 0;
		uint minimizer_position = 0;
		bool reversed = false;
		// Reverse complement if needed
		if (reversed) {
			cerr << "TODO: reverse complement" << endl;
		}
		// Create file for minimizer if missing
		if (outfiles.find(minimizer) == outfiles.end()) {
			// Create the file and write its header
			outfiles[minimizer] = new Kff_file(output_prefix + "_" + to_string(minimizer) + ".kff", "w");
			outfiles[minimizer]->write_encoding(0, 1, 3, 2);
			// Write the global variable needed
			Section_GV sgv(outfiles[minimizer]);
			sgv.write_var("k", k);
			sgv.write_var("m", m);
			sgv.write_var("max", 1);
			sgv.write_var("data_size", 0);
			sgv.close();
			// Open the minimizer block section
			outsections[minimizer] = new Section_Minimizer(outfiles[minimizer]);
		}
		// Verify k size in the output file
		if (outfiles[minimizer]->global_vars["k"] != k) {
			// Close current block section
			outsections[minimizer]->close();
			delete outsections[minimizer];
			// Update k
			Section_GV sgv(outfiles[minimizer]);
			sgv.write_var("k", k);
			sgv.close();
			// Open a new block section
			outsections[minimizer] = new Section_Minimizer(outfiles[minimizer]);
		}
		// Write the minimizer in the right file
		outsections[minimizer]->write_compacted_sequence (bin, k, minimizer_position, nullptr);
	}

	delete[] bin;
	// close all the outfiles
	for (auto & it : outsections) {
		// Close the section
		it.second->close();
		delete it.second;
		// Close the file
		outfiles[it.first]->close();
		delete outfiles[it.first];
	}
}


void Convert::monofile() {
	// prepare outfile
	Kff_file outfile(output_prefix + ".kff", "w");
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
			delete[] bin;
			bin = new uint8_t[k / 4 + 1];
		}

		brz.translate(kmer, bin);
		sr.write_compacted_sequence(bin, k, nullptr);
	}

	delete[] bin;
	sr.close();
	outfile.close();
}
