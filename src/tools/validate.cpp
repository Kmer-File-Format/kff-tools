#include <vector>
#include <string>
#include <cstring>

#include "validate.hpp"
#include "sequences.hpp"


using namespace std;



Validate::Validate() : input_filename(" ")
										 , index_only(false)
										 , verbose(false)
										 , sort_verif(false)
{};

void Validate::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("validate", "Read the input file and tell if some values are unexpected.");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "The file to copy");
	input_option->required();
	input_option->check(CLI::ExistingFile);

	subapp->add_flag("--index-only", index_only, "Restrict the verification to the index only.");
	subapp->add_flag("-v, --verbose", verbose, "Print all the validation process instead of only unexpected values.");
	subapp->add_flag("-s, --is-sorted", sort_verif, "Verify if the kmers in each seactions are sorted. The tool looks for encoding order. For R sections, blocks have to be composed of 1 kmer. For M sections, the verification is done by column. The comparison is done between kmers sharing the same minimizer position.");
}

void Validate::exec() {
	try {
		Kff_file infile(input_filename, "r");
		this->strif = Stringifyer(infile.encoding);

		if (verbose) {
			cout << "=== Header ===" << endl;
			// Important values
			cout << "-> KFF file version " << (uint)infile.major_version << "." << (uint)infile.minor_version << endl;
			cout << "-> Encoding: A=" << (uint)infile.encoding[0] << " C=" << (uint)infile.encoding[1] << " G=" << (uint)infile.encoding[2] << " T=" << (uint)infile.encoding[3] << endl;
			cout << "-> Uniqueness: " << (uint)infile.uniqueness << endl;
			cout << "-> Canonicity: " << (uint)infile.canonicity << endl;


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
				if (not this->is_valid_i_section(infile))
					exit(1);
			}
			// Raw sequence section
			else if (section_type == 'r') {
				if (not this->is_valid_r_section(infile))
					exit(1);
			}
			// Minimizer sequence section
			else if (section_type == 'm') {
				if (not this->is_valid_m_section(infile))
					exit(1);
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

bool Validate::is_valid_r_section(Kff_file & infile) {
	Section_Raw sr(&infile);

	uint k = infile.global_vars["k"];
	uint max = infile.global_vars["max"];
	uint data_size = infile.global_vars["data_size"];
	bool verif_sorting = this->sort_verif;

	if (verbose) {
		cout << "Start Byte " << sr.beginning << endl;
		cout << "-> Number of blocks: " << sr.nb_blocks << endl;
	}

	uint max_nucl = k + max - 1;
	uint64_t seq_max_bytes = (max_nucl + 3) / 4;
	uint8_t * seq_bytes = new uint8_t[seq_max_bytes];
	uint8_t * prev_bytes = new uint8_t[seq_max_bytes];
	memset(prev_bytes, 0, seq_max_bytes);
	uint8_t * data_bytes = new uint8_t[data_size * max];

	for (uint64_t i=0 ; i<sr.nb_blocks ; i++) {
		if (infile.tellp() >= infile.end_position) {
			cerr << "/!\\ End of the file reached before the end of the section." << endl;
			exit(1);
		}
		uint nb_kmers = sr.read_compacted_sequence(seq_bytes, data_bytes);

		if (nb_kmers == 0)
		{
			cerr << "Block containing 0 kmer detected." << endl;
			exit(1);
		}
		// Sort verif
		else if (verif_sorting)
		{
			if (nb_kmers > 1)
			{
				cerr << "WARNING: R block at position " << sr.beginning << " not sorted." << endl;
				cerr << "\tcause: compacted kmers cannot be sorted" << endl;
				verif_sorting = false;
			}
			else
			{
				int comp = sequence_compare( prev_bytes, k, 0, k-1,
                                             seq_bytes , k, 0, k-1 );
				if (comp > 0)
				{
					cerr << "WARNING: R block at position " << sr.beginning << " not sorted." << endl;
					cerr << "\tcause: kmers not in the encoging order:" << endl;
					cerr << strif.translate(prev_bytes, k) << endl;
					cerr << strif.translate(seq_bytes, k) << endl;
					verif_sorting = false;
				}
			}
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

		uint8_t * tmp = prev_bytes;
		prev_bytes = seq_bytes;
		seq_bytes = tmp;
	}

	delete[] seq_bytes;
	delete[] prev_bytes;
	delete[] data_bytes;

	return true;
}


bool Validate::is_valid_m_section(Kff_file & infile) {
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
			return false;
		}
		uint64_t mini_pos;
		uint nb_kmers = sm.read_compacted_sequence_without_mini(seq_bytes, data_bytes, mini_pos);

		if (nb_kmers == 0) {
			cerr << "Block containing 0 kmer detected." << endl;
			return false;
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

	return true;		
}


bool Validate::is_valid_i_section(Kff_file & infile) {
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
			return false;
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
			return false;
		}
		// Go back to previous position
		infile.jump_to(current_pos);
	}

	si.close();
	return true;
}

