#include <vector>
#include <string>

#include "merge.hpp"


using namespace std;


Merge::Merge() {
	output_filename = "";
}

void Merge::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("merge", "Merge a list of kff files into one. All the files must have the same encoding.");
	CLI::Option * input_option = subapp->add_option("-i, --inputs", input_filenames, "A list of input file names. The order of the file list will be preserved in the output file.");
	input_option->required();
	input_option->expected(2, -1);

	CLI::Option * out_option = subapp->add_option("-o, --outfile", output_filename, "Kff file where all the input will be merged");
	out_option->required();
}

void Merge::exec() {
	// Useful variables
	long buffer_size = 1048576; // 1 MB
	char buffer[1048576];
	uint8_t global_encoding[4];

	// Read the encoding of the first file and push it as outcoding
	Kff_file infile(input_filenames[0], "r");
	for (uint i=0 ; i<4 ; i++)
		global_encoding[i] = infile.encoding[i];
	infile.close();

	// Write header of the output
	Kff_file outfile(output_filename, "w");
	outfile.write_encoding(
		global_encoding[0],
		global_encoding[1],
		global_encoding[2],
		global_encoding[3]
	);
	// Set metadata
	std::string meta = "Merged file";
	outfile.write_metadata(meta.length(), (uint8_t *)meta.c_str());

	// Append each file one by one
	for (string in_filename : input_filenames) {
		// Open the file
		Kff_file infile(in_filename, "r");

		// Encoding verification
		for (uint i=0 ; i<4 ; i++) {
			if (infile.encoding[i] != global_encoding[i]) {
				cerr << "Wrong encoding for file " << in_filename << endl;
				cerr << "Its nucleotide encoding is different from previous kff files." << endl;
				cerr << "Please first use 'kff-tools translate' to have the same encoding" << endl;
				exit(1);
			}
		}

		// NB: Automatic jump over metadata due to API

		// Read section by section
		char section_type = infile.read_section_type();
		while(not infile.fs.eof()) {
			vector<string> to_copy;
			long size, end_byte, begin_byte;

			switch (section_type) {
				// Write the variables that change from previous sections (possibly sections from other input files)
				case 'v':
				{
					// Read the values
					Section_GV in_sgv(&infile);

					// Verify the presence and value of each variable in output
					for (auto& tuple : in_sgv.vars) {
						if (outfile.global_vars.find(tuple.first) == outfile.global_vars.end()
								or outfile.global_vars[tuple.first] != tuple.second)
							to_copy.push_back(tuple.first);
					}
				}
				// Create a global variable section if needed
				if (to_copy.size() > 0) {
					Section_GV sgv(&outfile);
					// Write variables
					for (string s : to_copy)
						sgv.write_var(s, infile.global_vars[s]);
					sgv.close();
				}
				break;

				// copy the sequence section from input to output
				case 'r':
				case 'm':
				// Analyse the section size
				begin_byte = infile.fs.tellp();
				infile.jump_next_section();
				end_byte = infile.fs.tellp();
				size = end_byte - begin_byte;
				infile.fs.seekp(begin_byte);

				// Read from input and write into output
				while (size > 0) {
					size_t size_to_copy = size > buffer_size ? buffer_size : size;

					infile.fs.read(buffer, size_to_copy);
					outfile.fs.write(buffer, size_to_copy);

					size -= size_to_copy;
				}
				break;

				default:
					cerr << "Unknown section type " << section_type << " in file " << in_filename << endl;
					exit(2);
			}

			// Prepare next section
			section_type = infile.read_section_type();
		}

		infile.close();
	}

	outfile.close();
}
