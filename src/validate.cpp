#include <vector>
#include <string>

#include "validate.hpp"
#include "encoding.hpp"


using namespace std;



Validate::Validate() {
	input_filename = "";
	verbose = false;
}

void Validate::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("validate", "Read the input file and tell if some values are unexpected.");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "The file to copy");
	input_option->required();
	input_option->check(CLI::ExistingFile);

	subapp->add_flag("-v, --verbose", verbose, "Print all the validation process instead of only unexpected values.");
}

void Validate::exec() {
	try {
		Kff_file infile(input_filename, "r");
		Stringifyer strif(infile.encoding);

		if (verbose) {
			cout << "=== Header ===" << endl;
			// Important values
			cout << "-> KFF file version " << (uint)infile.major_version << "." << (uint)infile.minor_version << endl;
			cout << "-> Encoding: A=" << (uint)infile.encoding[0] << " C=" << (uint)infile.encoding[1] << " G=" << (uint)infile.encoding[2] << " T=" << (uint)infile.encoding[3] << endl;
			cout << "-> Uniqueness: " << (uint)infile.uniqueness << endl;


			// Metadata
			cout << "-> Metadata (" << ((uint)infile.metadata_size) << "B)" << endl;
			uint8_t * metadata = new uint8_t[infile.metadata_size];
			infile.read_metadata(metadata);
			for (uint i=0 ; i<infile.metadata_size ; i++) {
				cout << (uint)metadata[i] << "\t";
				if (i % 16 == 15)
					cout << endl;
			}
			if (infile.metadata_size % 16 != 0)
				cout << endl;
			delete[] metadata;

			cout << endl;
		}

		char section_type = infile.read_section_type();
		while(infile.tellp() < infile.end_position) {
			if (verbose)
				cout << "=== Section " << section_type << " ===" << endl;

			// Global variable section
			if (section_type == 'v') {
				Section_GV sgv(&infile);

				if (verbose) {
					cout << "Start Byte " << sgv.beginning << endl;
					cout << "-> " << sgv.nb_vars << " variables" << endl;
					for (auto tuple : sgv.vars) {
						cout << tuple.first << " = " << tuple.second << endl;
					}
				}

				sgv.close();
			}
			// Index section
			else if (section_type == 'i') {
				Section_Index si(&infile);
				long end_byte = si.beginning + 17 + 9 * si.index.size();

				if (this->verbose) {
					cout << "Start Byte " << si.beginning << endl;
					cout << "Section\trelative\tabsolute" << endl;
				}

				for (const auto & pair : si.index) {
					if (this->verbose)
						cout << pair.second << "\t" << pair.first << "\t" << (end_byte + pair.first) << endl;

					// Jump to the section
					long section_pos = end_byte + pair.first;
					long current_pos = infile.tellp();
					infile.jump_to(section_pos);
					// Read the byte at this position
					uint8_t type = 0;
					infile.read(&type, 1);
					if (type != pair.second) {
						cerr << "Wrong section at position " << section_pos << ". Found a section " << type << endl;
						exit(1);
					}
					// Go back to previous position
					infile.jump_to(current_pos);
				}

				if (this->verbose) {
					cout << "Next index position " << si.next_index << endl;
				}
				
				if (si.next_index != 0) {
					// Jump to the section
					long section_pos = end_byte + si.next_index;
					long current_pos = infile.tellp();
					infile.jump_to(section_pos);
					// Read the byte at this position
					uint8_t type = 0;
					infile.read(&type, 1);
					if (type != 'i') {
						cerr << end_byte << " " << si.next_index << endl;
						cerr << "No index found at position " << section_pos << "." << endl;
						exit(1);
					}
					// Go back to previous position
					infile.jump_to(current_pos);	
				}

				si.close();
			}
			// Raw sequence section
			else if (section_type == 'r') {
				Section_Raw sr(&infile);

				uint k = infile.global_vars["k"];
				uint max = infile.global_vars["max"];
				uint data_size = infile.global_vars["data_size"];

				if (verbose) {
					cout << "Start Byte " << sr.beginning << endl;
					cout << "-> Number of blocks: " << sr.nb_blocks << endl;
				}

				uint max_nucl = k + max - 1;
				uint8_t * seq_bytes = new uint8_t[max_nucl / 4 + 1];
				uint8_t * data_bytes = new uint8_t[data_size * max];

				for (uint i=0 ; i<sr.nb_blocks ; i++) {
					if (infile.tellp() >= infile.end_position) {
						cerr << "/!\\ End of the file reached before the end of the section." << endl;
						exit(1);
					}
					uint nb_kmers = sr.read_compacted_sequence(seq_bytes, data_bytes);

					if (nb_kmers == 0) {
						cerr << "Block containing 0 kmer detected." << endl;
						exit(1);
					}

					if (verbose) {
						cout << "* Number of kmers: " << nb_kmers << endl;
						cout << strif.translate(seq_bytes, k + nb_kmers - 1) << endl;

						if (data_size != 0) {
							cout << "data array: ";
							for (uint i_data=0 ; i_data<nb_kmers ; i_data++) {
								if (i_data > 0)
									cout << ",\t";

								cout << "[";
								for (uint b_idx=0 ; b_idx<data_size ; b_idx++) {
									if (b_idx > 0)
										cout << "\t";
									cout << (uint)data_bytes[data_size*i_data+b_idx];
								}
								cout << "]";
							}
							cout << endl;
						}
					}
				}

				delete[] seq_bytes;
				delete[] data_bytes;
			}
			// Minimizer sequence section
			else if (section_type == 'm') {
				Section_Minimizer sm(&infile);

				uint k = infile.global_vars["k"];
				uint m = infile.global_vars["m"];
				uint max = infile.global_vars["max"];
				uint data_size = infile.global_vars["data_size"];

				if (verbose) {
					cout << "Start Byte " << sm.beginning << endl;
					cout << "-> Minimizer: " << strif.translate(sm.minimizer, m) << endl;
					cout << "-> Number of blocks: " << sm.nb_blocks << endl;
					cout << endl;
				}

				uint max_nucl = k - m + max - 1;
				uint8_t * seq_bytes = new uint8_t[max_nucl / 4 + 1];
				uint8_t * data_bytes = new uint8_t[data_size * max];

				for (uint i=0 ; i<sm.nb_blocks ; i++) {
					if (infile.tellp() >= infile.end_position) {
						cerr << "/!\\ End of the file reached before the end of the section." << endl;
						exit(1);
					}
					uint64_t mini_pos;
					uint nb_kmers = sm.read_compacted_sequence_without_mini(seq_bytes, data_bytes, mini_pos);

					if (nb_kmers == 0) {
						cerr << "Block containing 0 kmer detected." << endl;
						exit(1);
					}

					if (mini_pos > k - m + nb_kmers) {
						cerr << "* minimizer position out of sequence. position=" << mini_pos << " , skmer_size=" <<  (k-m+nb_kmers) << endl;
					}

					if (verbose) {
						cout << "* minimizer position: " << mini_pos << "\tNumber of kmers: " << nb_kmers << endl;
						string seq = strif.translate(seq_bytes, k - m + nb_kmers - 1);
						cout << seq.substr(0, mini_pos) << "|" << seq.substr(mini_pos, k - m + nb_kmers - 1 - mini_pos) << endl;
						if (data_size != 0) {
							cout << "data array: ";
							for (uint i_data=0 ; i_data<nb_kmers ; i_data++) {
								cout << "[";
								for (uint b_idx=0 ; b_idx<data_size ; b_idx++) {
									if (b_idx > 0)
										cout << "\t";
									cout << (uint)data_bytes[data_size*i_data+b_idx];
								}
								cout << "],\t";
							}
							cout << endl;
						}
					}
				}

				delete[] seq_bytes;
				delete[] data_bytes;
			}
			// Unknown section
			else {
				cerr << "/!\\ Unknown section " << section_type << " (uint value " << (uint)section_type << ")" << endl;
				exit(1);
			}

			if (verbose)
				cout << endl;
			section_type = infile.read_section_type();
		}

		char kff[3];
		infile.read((uint8_t *)kff, 3);
		if (kff[0] != 'K' or kff[1] != 'F' or kff[2] !='F')
			cout << "No KFF signature found at the end of the file. The file must be corrupted." << endl;

		if (infile.tellp() < infile.end_position + 3) {
			cout << "/!\\ Remaining bytes at the end of the file !" << endl;
		}
		infile.close();
		if (verbose)
			cout << "=== End of the file ===" << endl << endl;

		cout << "No problem detected" << endl;
	} catch (const char* msg) {
		cerr << "/!\\ " << msg << endl;
		exit(1);
	}
}
