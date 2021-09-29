#include <vector>
#include <string>

#include "datarm.hpp"


using namespace std;



DataRm::DataRm() {
	input_filename = "";
	output_filename = "";
}

void DataRm::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("data-rm", "Read the input file and rewrite it in the output, removing all the data associated with the kmers.");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "The file to copy");
	input_option->required();
	input_option->check(CLI::ExistingFile);

	CLI::Option * out_option = subapp->add_option("-o, --outfile", output_filename, "Outfile without data.");
	out_option->required();
}

void DataRm::exec() {
	// Rewrite the encoding
	Kff_file infile(input_filename, "r");
	Kff_file outfile(output_filename, "w");
	outfile.set_indexation(true);
	outfile.write_encoding(infile.encoding);

	// Rewrite metadata
	uint8_t * metadata = new uint8_t[infile.metadata_size];
	infile.read_metadata(metadata);
	outfile.write_metadata(infile.metadata_size, metadata);
	delete[] metadata;

	// Prepare sequence buffer
	uint8_t * nucleotides = new uint8_t[1];
	uint8_t * data = new uint8_t[1];
	uint64_t real_data_size = 1;

	// Read and write section per section
	char section_type = infile.read_section_type();
	while (infile.tellp() != infile.end_position) {
		// Read variables
		if (section_type == 'v') {
			// Load variables
			Section_GV isgv(&infile);

			// Remove old footer
			if (isgv.vars.find("footer_size") != isgv.vars.end())
				continue;

			Section_GV osgv(&outfile);

			bool nucl_buffer_changed = false;
			bool data_buffer_changed = false;
			for (auto var_tuple : isgv.vars) {
				if (var_tuple.first == "data_size") {
					osgv.write_var("data_size", 0);
					real_data_size = var_tuple.second;
				} else {
					osgv.write_var(var_tuple.first, var_tuple.second);

					if (var_tuple.first == "k" or var_tuple.first == "max")
						nucl_buffer_changed = true;
					if (var_tuple.first == "data_size" or var_tuple.first == "max")
						data_buffer_changed = true;
				}
			}
			isgv.close();
			osgv.close();

			// Buffer update
			if (nucl_buffer_changed) {
				delete[] nucleotides;
				uint max_nucl = outfile.global_vars["k"] + outfile.global_vars["max"] - 1;
				nucleotides = new uint8_t[max_nucl / 4 + 1];
			}
			if (data_buffer_changed) {
				delete[] data;
				uint nb_data = outfile.global_vars["max"];
				data = new uint8_t[nb_data * real_data_size];
			}
		}
		else if (section_type == 'i') {
			Section_Index si(&infile);
			si.close();
		}
		// rewrite a raw block
		else if (section_type == 'r') {
			// Open sections
			Section_Raw in_section(&infile);
			Section_Raw out_section(&outfile);

			// Rewrite block per block
			uint64_t k = outfile.global_vars["k"];
			for (uint i=0 ; i<in_section.nb_blocks ; i++) {
				uint nb_kmers = in_section.read_compacted_sequence(nucleotides, data);
				out_section.write_compacted_sequence(nucleotides, k + nb_kmers - 1, nullptr);
			}

			in_section.close();
			out_section.close();
		}
		// Revwrite a minimizer block
		else if (section_type == 'm') {
			// Open sections
			Section_Minimizer in_section(&infile);
			Section_Minimizer out_section(&outfile);

			// Rewrite the minimizer
			uint k = outfile.global_vars["k"];
			uint m = outfile.global_vars["m"];
			out_section.write_minimizer(in_section.minimizer);

			// Rewrite block per block
			for (uint i=0 ; i<in_section.nb_blocks ; i++) {
				// Read
				uint64_t mini_pos;
				uint nb_kmers = in_section.read_compacted_sequence_without_mini(nucleotides, data, mini_pos);
				// Write
				uint nucl_size = k - m + nb_kmers - 1;
				out_section.write_compacted_sequence_without_mini(nucleotides, nucl_size, mini_pos, nullptr);
			}

			in_section.close();
			out_section.close();
		}

		section_type = infile.read_section_type();
	}

	delete[] nucleotides;
	delete[] data;
	infile.close();
	outfile.close();
}
