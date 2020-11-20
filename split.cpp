#include "split.hpp"


using namespace std;


Split::Split() {
	// Paths
	this->input_filename = "";
	this->output_dirname = "./";
}

void Split::cli_prepare(CLI::App * subapp) {
	CLI::Option * input_option = subapp->add_option("-i, --input", input_filename, "Input file in kff format");
	input_option->required();
	subapp->add_option("-o, --outdir", output_dirname, "Output directory where to put all the subkff files");
}

void Split::exec() {
	// IO Prepare
	if (output_dirname[output_dirname.length()-1] != '/')
		output_dirname += "/";

	Kff_file input_file(input_filename, "r");
	input_file.read_encoding();
	uint8_t input_metadata[1024];
	uint32_t metadata_size = input_file.size_metadata();
	input_file.read_metadata(metadata_size, input_metadata);

	long buffer_size = 1048576; // 1 MB
	char buffer[1048576];

	uint idx = 0;
	char section_type = input_file.read_section_type();
	while (not input_file.fs.eof()) {
		if (section_type == 'v') {
			input_file.open_section_GV();
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
			output_file.write_metadata(metadata_size, input_metadata);

			// Write needed variables
			Section_GV sgv = output_file.open_section_GV();
			sgv.write_var("k", input_file.global_vars["k"]);
			sgv.write_var("max", input_file.global_vars["max"]);
			sgv.write_var("data_size", input_file.global_vars["data_size"]);
			if (section_type == 'm')
				sgv.write_var("m", input_file.global_vars["m"]);
			sgv.close();

			// Analyse the section size
			auto begin_byte = input_file.fs.tellp();
			Block_section_reader * section = Block_section_reader::construct_section(section_type, &input_file);
			section->jump_section();
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
}
