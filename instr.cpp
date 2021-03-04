#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <cstring>
#include <cstdio>

#include "instr.hpp"
#include "encoding.hpp"
#include "merge.hpp"
#include "sequences.hpp"


using namespace std;



Instr::Instr() {
	input_filename = "";
	output_prefix = "";
	m = 0;
	split = false;
	data_size = 0;
}

void Instr::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("instr", "Convert a text kmer file into one or multuple kff file(s). Instr suppose that data are signed integers (max 64 bits).");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "A text kmer file in tsv format. One kmer per line and data as the second field (can be omitted for no data)");
	input_option->required();
	CLI::Option * output_option = subapp->add_option("-o, --outprefix", output_prefix, "The kff output file prefix. For a single file the name will be <prefix>.kff, otherwise, one file per minimizer is created (<prefix>_<minimizer>.kff).");
	output_option->required();
	subapp->add_option("-m, --mini-size", m, "Set a minimizer size. If set, all the kmers are separated into one section per lexicagraphic minimizer. (Max 31)");
	subapp->add_flag("-s, --split", split, "The output is splited into one file per minimizer.");
	subapp->add_option("-d, --data-size", data_size, "Data size in Bytes (Default 0, max 8).");
}

void Instr::exec() {
	if (m == 0) {
		this->monofile();
	} else {
		this->multifile();
	}
}

void Instr::multifile() {
	vector<string> filenames;
	map<uint, Kff_file *> outfiles;
	map<uint, Section_Minimizer *> outsections;

	// Prepare first k value
	ifstream txt_file(input_filename);
	string line;
	getline(txt_file, line);
	txt_file.close();

	uint64_t k = line.substr(0, line.find("\t")).length();

	// Prepare the binarized datastructures
	uint8_t encoding[4] = {0, 1, 3, 2};
	Binarizer brz(encoding);
	uint8_t * bin = new uint8_t[k / 4 + 1];
	uint8_t * bin_mini = new uint8_t[m / 4 + 1];
	uint8_t * data = new uint8_t[data_size];
	// Read the input kmer stream
	txt_file = ifstream(input_filename);
	while (getline(txt_file, line)) {
		string kmer = line.substr(0, line.find("\t"));
		uint64_t data_int = stoi(line.substr(line.find("\t")+1));

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
		// bool reversed = false;

		search_mini(bin, k, m, minimizer, minimizer_position);

		// Create file for minimizer if missing
		if (outfiles.find(minimizer) == outfiles.end()) {
			// Create the file and write its header
			string filename = output_prefix + "_" + to_string(minimizer) + ".kff";
			filenames.push_back(filename);
			outfiles[minimizer] = new Kff_file(filename, "w");
			outfiles[minimizer]->write_encoding(0, 1, 3, 2);
			// Write the global variable needed
			Section_GV sgv(outfiles[minimizer]);
			sgv.write_var("k", k);
			sgv.write_var("m", m);
			sgv.write_var("max", 1);
			sgv.write_var("data_size", data_size);
			sgv.close();
			// Open the minimizer block section
			outsections[minimizer] = new Section_Minimizer(outfiles[minimizer]);
			uint_to_seq(minimizer, bin_mini, m);
			outsections[minimizer]->write_minimizer(bin_mini);
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
			uint_to_seq(minimizer, bin_mini, m);
			outsections[minimizer]->write_minimizer(bin_mini);
		}
		// Translate int values
		for (int i=data_size-1 ; i>=0 ; i--) {
			data[i] = data_int & 0xFF;
			data_int >>= 8;
		}
		// Write the minimizer in the right file
		outsections[minimizer]->write_compacted_sequence (bin, k, minimizer_position, data);
	}

	delete[] bin;
	delete[] bin_mini;
	delete[] data;
	// close all the outfiles
	for (auto & it : outsections) {
		// Close the section
		it.second->close();
		delete it.second;
		// Close the file
		outfiles[it.first]->close();
		delete outfiles[it.first];
	}

	if (not split) {
		Merge merger;
		merger.merge(filenames, output_prefix + ".kff");

		for (string filename : filenames)
			remove(filename.c_str());
	}
}


void Instr::monofile() {
	// prepare outfile
	Kff_file outfile(output_prefix + ".kff", "w");
	outfile.write_encoding(0, 1, 3, 2);

	ifstream txt_file(input_filename);
	string line;
	getline(txt_file, line);
	txt_file.close();

	uint64_t k = line.substr(0, line.find("\t")).length();

	Section_GV sgv(&outfile);
	sgv.write_var("k", k);
	sgv.write_var("max", 1);
	sgv.write_var("data_size", data_size);
	sgv.close();
	
	Section_Raw sr(&outfile);

	txt_file = ifstream(input_filename);
	Binarizer brz(outfile.encoding);
	uint8_t * bin = new uint8_t[k / 4 + 1];
	uint8_t * data = new uint8_t[data_size];
	while (getline(txt_file, line)) {
		string kmer = line.substr(0, line.find("\t"));
		uint64_t data_int = stoi(line.substr(line.find("\t")+1));

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

		// Translate nucleotides
		brz.translate(kmer, bin);
		// Translate int values
		for (int i=data_size-1 ; i>=0 ; i--) {
			data[i] = data_int & 0xFF;
			data_int >>= 8;
		}
		// Write into the kff file
		sr.write_compacted_sequence(bin, k, data);
	}

	delete[] bin;
	delete[] data;
	sr.close();
	outfile.close();
}
