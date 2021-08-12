#include <vector>
#include <string>
#include <algorithm>

#include "sort.hpp"


using namespace std;


Sort::Sort () {
	output_filename = "";
}

void Sort::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("sort", "Sort the order of the sequences within each section (m and r) of a kff files.");
	CLI::Option * input_option = subapp->add_option("-i, --input", input_filename, "Input KFF file");
	input_option->required();

	CLI::Option * out_option = subapp->add_option("-o, --outfile", output_filename, "Sortd output KFF file");
	out_option->required();
}

// the code of this function is largely inspired by merge.cpp
void Sort::sort(string input, string output) {
	// Useful variables
	const long buffer_size = 1048576; // 1 MB
	uint8_t buffer[1048576];
	uint8_t global_encoding[4];

	// Read the encoding of the input file and push it as output encoding 
	Kff_file infile_enc(input, "r");
	for (uint i=0 ; i<4 ; i++)
		global_encoding[i] = infile_enc.encoding[i];
    infile_enc.close();

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
	std::string meta = "Sortd file";
	outfile.write_metadata(meta.length(), (uint8_t *)meta.c_str());


	// remember index previous position for chaining
	long last_index = 0;
	// Footers
	map<string, uint64_t> footer_values;

	    Kff_file infile(input, "r");

		// NB: Automatic jump over metadata due to API
		// Read section by section
		char section_type = infile.read_section_type();
		while(infile.tellp() != infile.end_position) {
			vector<string> to_copy;
			long size;

			switch (section_type) {
				// Write the variables that change from previous sections (possibly sections from other input files)
				case 'v':
				{
					// Read the values
					Section_GV in_sgv(&infile);

					// Discard footers
					if (in_sgv.vars.find("footer_size") != in_sgv.vars.end()) {
						for (auto& tuple : in_sgv.vars) {
							if (tuple.first != "footer_size" and tuple.first != "first_index") {
								if (footer_values.find(tuple.first) == footer_values.end())
									footer_values[tuple.first] = tuple.second;
								else
									// Sum up the common footer values
									footer_values[tuple.first] += tuple.second;
							}
						}
						break;
					}

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

				// process a raw sequence section
				case 'r':
                {
                    // get some params, allocate seq/data buffers
                    // can't do that earlier, as we need to read global variable section first
                    uint k = infile.global_vars["k"];
                    uint max = infile.global_vars["max"];
                    uint data_size = infile.global_vars["data_size"];
                    uint max_nucl = k + max - 1;
                    uint8_t * seq_bytes = new uint8_t[max_nucl / 4 + 1];
                    uint8_t * data_bytes = new uint8_t[data_size * max];
                    cerr << "new 'r' section found" << endl;
                    cerr << "Max sequence size: " << max_nucl << " data size: " << data_size<< endl;

                    vector<pair<vector<uint8_t>,vector<uint8_t>>> everything;

                    // Open sections
                    Section_Raw in_section(&infile);
                    Section_Raw out_section(&outfile);

                    // Read whole block and store everything inside a vector
                    for (uint i=0 ; i<in_section.nb_blocks ; i++) {
                        // Read the block.
                        string seq;
                        int nb_kmers = in_section.read_compacted_sequence(seq_bytes, data_bytes);
                        int seq_len = k + nb_kmers - 1;
                        //cerr << "Read seq of length: " << seq_len << " data size: " << data_size*nb_kmers << endl;
                        vector<uint8_t> seq_vec(&seq_bytes[0],&seq_bytes[seq_len/4+1]);
                        vector<uint8_t> data_vec(&data_bytes[0],&data_bytes[data_size*nb_kmers]);
                        //cerr << "Converted to pair of lengths " << seq_vec.size() << " / " << data_vec.size() << endl;
                        everything.push_back(make_pair(seq_vec,data_vec));
                    }

                    // Sort the vector
                    std::sort(everything.begin(),everything.end(), [](auto &left, auto &right) {
                                return left.first < right.first;
                                });

                    // Write the vector to out block
                    for (pair<vector<uint8_t>,vector<uint8_t>> elt : everything)
                    {
                        int nb_kmers =  elt.second.size() / data_size;
                        //cerr << "Writing seq of length: " << k + nb_kmers - 1 << " data size: " << elt.second.size() << endl;
                        out_section.write_compacted_sequence(
                                elt.first.data(),
                                k + nb_kmers - 1,
                                elt.second.data());
                    }
                    in_section.close();
                    out_section.close();
                    break;
                }

                /*
                 // not implemented for now, but here's the code I'd start from.
                case 'm':
                {

                    // Open sections
                    Section_Minimizer in_section(&infile);
                    Section_Minimizer out_section(&outfile);

                    vector<uint8_t *> saved_kmers;
                    vector<uint8_t *> saved_data;

                    // Rewrite the minimizer
                    uint k = outfile.global_vars["k"];
                    uint m = outfile.global_vars["m"];
                    uint data_size = outfile.global_vars["data_size"];
                    out_section.write_minimizer(in_section.minimizer);

                    // Rewrite block per block
                    for (uint i=0 ; i<in_section.nb_blocks ; i++) {
                        // Read
                        uint64_t mini_pos;
                        uint64_t nb_kmers = in_section.read_compacted_sequence_without_mini(seq, data, mini_pos);

                        out_section.write_compacted_sequence_without_mini(
                                seq,
                                seq.length() - m,
                                mini_pos,
                                data
                                );
                    }

                    break;
                }
                */

                // nb: is probably overkill code here, taken from disjoin.cpp, could be simplified as we're processing just 1 file
                case 'i': {
				// read section and compute its size
				Section_Index si(&infile);
				si.close();
				long file_size = infile.tellp() - si.beginning - 8l;
				infile.jump(-file_size - 8);

				// Save the position in the file for later chaining
				long i_position = outfile.tellp();
				size = file_size;
				// Copy section (except the chaining part)
				// Read from input and write into output
				while (size > 0) {
					size_t size_to_copy = size > buffer_size ? buffer_size : size;

					infile.read(buffer, size_to_copy);
					outfile.write(buffer, size_to_copy);

					size -= size_to_copy;
				}
				// Jump over the last value of infile
				infile.jump(8);
				// Chain the section and save its position
				long i_relative = last_index - (i_position + file_size + 8l);
				for (uint i=0 ; i<8 ; i++) {
					uint8_t val = (uint8_t)(i_relative >> (56 - 8 * i));
					outfile.write(&val, 1);
				}
				// write_value(last_index, outfile.fs);
				last_index = i_position;
				} break;

				default:
					cerr << "Unsupported section type " << section_type << " in file " << input << endl;
					exit(2);
			}

			// Prepare next section
			section_type = infile.read_section_type();
		}

		infile.close();

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

void Sort::exec() {
	this->sort(input_filename, output_filename);
}
