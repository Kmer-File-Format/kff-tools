#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <queue>

#include "encoding.hpp"
#include "sequences.hpp"
#include "compact.hpp"
#include "merge.hpp"


using namespace std;


Compact::Compact() {
	input_filename = "";
	output_filename = "";
	split = false;
	m = 0;

	load_mem_size = 1;
	loading_memory = new uint8_t[load_mem_size];
	kmer_buffer = new uint8_t[1];
	skmer_buffer = new uint8_t[1];
	data_buffer = new uint8_t[1];
}


Compact::~Compact() {
	delete[] loading_memory;
	delete[] kmer_buffer;
	delete[] skmer_buffer;
	delete[] data_buffer;
}


void Compact::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("compact", "Read a kff file and try to compact the kmers from minimizer sections. The available ram must be sufficent to load a complete minimizer section into memory.");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "Input kff file to compact.");
	input_option->required();
	CLI::Option * out_option = subapp->add_option("-o, --outfile", output_filename, "Kff to write (must be different from the input)");
	out_option->required();
}


void Compact::compact() {
	Kff_file infile(input_filename, "r");
	Kff_file outfile(output_filename, "w");

	// TODO: Write encoding
	// TODO: Set the flags
	// TODO: Write header

	while (input_file.fs.tellp() != size - 3) {
		char section_type = infile.read_section_type();

		if (section_type == 'v') {
			Section_GV sgv(&infile);
			sgv.close();
		}
		else if (section_type == 'i') {

		}
		else if (section_type == 'r') {

		}
		else if (section_type == 'm') {
			uint k = outfile.global_vars["k"];
			uint m = outfile.global_vars["m"];
			uint data_size = outfile.global_vars["data_size"];

			// Rewrite a value section if max is not sufficently large
			if (outfile.global_vars["max"] < (k-m) * 2) {
				map<string, uint64_t> values(outfile.global_vars);
				Section_GV sgv(&outfile);

				for (pair<string, uint64_t> & p : values)
					if (p.first != "max")
						sgv.write_var(p.first, p.second);
				sgv.write_var("max", (k-m)*2);

				sgv.close();
			}

			// Compact and save the kmers
			Section_Minimizer sm(&infile);
			this->compact_section(sm, outfile);
			sm.close();
		}
	}

	infile.close();
	outfile.close();
}

void Compact::compact_section(const Section_Minimizer & ism, Kff_file & outfile) {
	// General variables
	uint k = outfile.global_vars["k"];
	uint m = outfile.global_vars["m"];
	uint data_size = outfile.global_vars["data_size"];

	// Buffers
	uint8_t * seq_buffer = new uint8_t[(2 * k + 3) / 4];
	uint8_t * data_buffer = new uint8_t[2 * k * data_size];
	uint64_t mini_pos;

	// Kmer storing space
	uint kmer_buffer_size = 10000;
	uint next_free = 0;
	uint8_t * kmers = malloc(10000);
	vector<vector<uint8_t *> > kmers_per_index;

	// 1 - Load the input section
	for (uint n=0 ; n<ism.nb_blocks ; n++) {
		uint nb_kmers = read_compacted_sequence_without_mini(
			seq_buffer, data_buffer, mini_pos);
	}

	free(kmers);

	// 2 - Compact kmers


	// Cleaning
	delete[] seq_buffer;
	delete[] data_buffer;
}

