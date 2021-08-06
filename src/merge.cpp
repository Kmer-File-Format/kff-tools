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

void Merge::merge(vector<string> inputs, string output) {
	// Useful variables
	const long buffer_size = 1048576; // 1 MB
	char buffer[1048576];
	uint8_t global_encoding[4];

	// Read the encoding of the first file and push it as outcoding
	Kff_file infile(inputs[0], "r");
	for (uint i=0 ; i<4 ; i++)
		global_encoding[i] = infile.encoding[i];
	infile.close();

	// Write header of the output
	Kff_file outfile(output, "w");
	outfile.set_indexation(false);
	outfile.write_encoding(
		global_encoding[0],
		global_encoding[1],
		global_encoding[2],
		global_encoding[3]
	);
	// Set metadata
	std::string meta = "Merged file";
	outfile.write_metadata(meta.length(), (uint8_t *)meta.c_str());

	// remember index previous position for chaining
	long last_index = 0;
	// Footers
	map<string, uint64_t> footer_values;

	// Append each file one by one
	for (string in_filename : inputs) {
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

		long current_pos = infile.fs.tellp();
		infile.fs.seekg(0, infile.fs.end);
		long filesize = infile.fs.tellp();
		infile.fs.seekp(current_pos);

		// NB: Automatic jump over metadata due to API
		// Read section by section
		char section_type = infile.read_section_type();
		while(infile.fs.tellp() != filesize - 3) {
			vector<string> to_copy;
			long size, end_byte, begin_byte;

			switch (section_type) {
				// Write the variables that change from previous sections (possibly sections from other input files)
				case 'v':
				{
					unordered_map<string, uint64_t> variables;

					// Read variables
					while (section_type == 'v') {
						Section_GV sgv(&infile);

						// Is it a footer ?
						if (sgv.vars.find("footer_size") != sgv.vars.end()) {
							for (auto& tuple : sgv.vars)
								if (tuple.first != "footer_size" and tuple.first != "first_index") {
									if (footer_values.find(tuple.first) == footer_values.end())
										footer_values[tuple.first] = tuple.second;
									else
										// Sum up the common footer values
										footer_values[tuple.first] += tuple.second;
								}
							break;
						}
						// Not a footer
						for (auto & p : sgv.vars)
							variables[p.first] = p.second;
						sgv.close();

						// Update section_type
						section_type = infile.read_section_type();
					}

					// Does it need a rewrite ?
					bool v_section_needed = false;
					for (auto & p : variables) {
						if (outfile.global_vars.find(p.first) == outfile.global_vars.end()
								or outfile.global_vars[p.first] != p.second) {
							v_section_needed = true;
							break;
						}
					}

					// cout << "V needed ?"
					// Rewrite
					if (v_section_needed) {
						Section_GV sgv(&outfile);
						for (auto & p : variables) {
							sgv.write_var(p.first, p.second);
						}
						sgv.close();
					}
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
				case 'i': {
				// read section and compute its size
				Section_Index si(&infile);
				si.close();
				long file_size = infile.fs.tellp() - si.beginning - 8l;
				infile.fs.seekp(si.beginning);

				// Save the position in the file for later chaining
				long i_position = outfile.fs.tellp();
				size = file_size;
				// Copy section (except the chaining part)
				// Read from input and write into output
				while (size > 0) {
					size_t size_to_copy = size > buffer_size ? buffer_size : size;

					infile.fs.read(buffer, size_to_copy);
					outfile.fs.write(buffer, size_to_copy);

					size -= size_to_copy;
				}
				// Jump over the last value of infile
				infile.fs.seekp(infile.fs.tellp() + 8l);
				// Chain the section and save its position
				long i_relative = last_index - (i_position + file_size + 8l);
				for (uint i=0 ; i<8 ; i++) {
					char val = (char)(i_relative >> (56 - 8 * i));
					outfile.fs.write(&val, 1);
				}
				// write_value(last_index, outfile.fs);
				last_index = i_position;
				} break;

				default:
					cerr << "Unknown section type " << section_type << " in file " << in_filename << endl;
					exit(2);
			}

			// Prepare next section
			section_type = infile.read_section_type();
		}

		infile.close();
	}

	// Write footer
	if (last_index != 0) {
		Section_GV sgv(&outfile);
		long size = 9;

		for (const auto & tuple : footer_values) {
			sgv.write_var(tuple.first, tuple.second);
			size += tuple.first.length() + 1 + 8;
		}
		sgv.write_var("first_index", last_index);
		sgv.write_var("footer_size", size + 2 * (12 + 8));
		sgv.close();

	}

	outfile.close();
}

void Merge::exec() {
	this->merge(input_filenames, output_filename);
}
