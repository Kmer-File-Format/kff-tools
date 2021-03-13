#include "split.hpp"


using namespace std;


Split::Split() {
	// Paths
	this->input_filename = "";
	this->output_dirname = "./";
}

void Split::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("split", "Split a kff file into one file per section.");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "Input file in kff format");
	input_option->required();
	subapp->add_option("-o, --outdir", output_dirname, "Output directory where to put all the subkff files");
}

void Split::exec() {
	// IO Prepare
	if (output_dirname[output_dirname.length()-1] != '/')
		output_dirname += "/";

	Kff_file input_file(input_filename, "r");
	uint8_t * input_metadata = new uint8_t[input_file.metadata_size];
	input_file.read_metadata(input_metadata);

	long buffer_size = 1048576; // 1 MB
	char buffer[1048576];

	uint idx = 0;
	char section_type = input_file.read_section_type();
	while (not input_file.fs.eof()) {
		if (section_type == 'v') {
			Section_GV sgv(&input_file);
		} else {
			// Beginning of a section with kmers. Open a new kff file
			stringstream ss;
			ss << output_dirname << section_type << "_" << idx << ".kff";
			Kff_file output_file(ss.str(), "w");

			// Write header
			output_file.write_encoding(
				input_file.encoding[0],
				input_file.encoding[1],
				input_file.encoding[2],
				input_file.encoding[3]
			);
			output_file.write_metadata(input_file.metadata_size, input_metadata);

			// Write needed variables
			Section_GV sgv(&output_file);
			sgv.write_var("k", input_file.global_vars["k"]);
			sgv.write_var("max", input_file.global_vars["max"]);
			sgv.write_var("data_size", input_file.global_vars["data_size"]);
			if (section_type == 'm')
				sgv.write_var("m", input_file.global_vars["m"]);
			sgv.close();

			// Analyse the section size
			auto begin_byte = input_file.fs.tellp();
			if (not input_file.jump_next_section()) {
				cerr << "Error inside of the input file." << endl;
				cerr << "Impossible to jump over the section" << endl;
				exit(1);
			}
			auto end_byte = input_file.fs.tellp();
			long size = end_byte - begin_byte;
			input_file.fs.seekp(begin_byte);

			// Read from input and write into output
			while (size > 0) {
				size_t to_copy = size > buffer_size ? buffer_size : size;

				input_file.fs.read(buffer, to_copy);
				output_file.fs.write(buffer, to_copy);

				size -= to_copy;
			}

			output_file.close();
			idx += 1;
		}
		section_type = input_file.read_section_type();
	}

	delete[] input_metadata;
}
