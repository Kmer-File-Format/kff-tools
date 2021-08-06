#include <vector>
#include <string>
#include <cstring>

#include "disjoin.hpp"
#include "sequences.hpp"


using namespace std;



Disjoin::Disjoin() {
	input_filename = "";
	output_filename = "";
}

void Disjoin::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("disjoin", "Read the input file, and for each block in each section, rewrite n blocks, where n is the number of kmer in the block. The output file contains exactly one kmer per block.");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "The file to disjoin");
	input_option->required();

	CLI::Option * out_option = subapp->add_option("-o, --outfile", output_filename, "Disjoined file.");
	out_option->required();
}

void Disjoin::exec() {
	// Rewrite the encoding
	Kff_file infile(input_filename, "r");
	Kff_file outfile(output_filename, "w");

	outfile.write_encoding(infile.encoding);

	// Rewrite metadata
	uint8_t * metadata = new uint8_t[infile.metadata_size];
	infile.read_metadata(metadata);
	outfile.write_metadata(infile.metadata_size, metadata);
	delete[] metadata;

	// Compute file size
	long current_pos = infile.fs.tellp();
	infile.fs.seekg(0, infile.fs.end);
	long size = infile.fs.tellp();
	infile.fs.seekp(current_pos);

	// Prepare sequence buffer
	uint8_t * nucleotide_shifts[4];
	for (uint i=0 ; i<4 ; i++)
		nucleotide_shifts[i] = new uint8_t[1];
	uint8_t * data = new uint8_t[1];
	uint64_t real_max = 1;

	// Read and write section per section
	char section_type = infile.read_section_type();
	while (infile.fs.tellp() != size - 3) {
		// Read variables
		if (section_type == 'v') {
			// Load variables
			Section_GV isgv(&infile);
			if (isgv.vars.find("footer_size") != isgv.vars.end()) {
				continue;
			}

			Section_GV osgv(&outfile);

			bool nucl_buffer_changed = false;
			bool data_buffer_changed = false;
			for (auto var_tuple : isgv.vars) {
				if (var_tuple.first == "max") {
					osgv.write_var("max", 1);
					real_max = var_tuple.second;
					nucl_buffer_changed = true;
					data_buffer_changed = true;
				} else {
					osgv.write_var(var_tuple.first, var_tuple.second);
					if (var_tuple.first == "k")
						nucl_buffer_changed = true;
					if (var_tuple.first == "data_size")
						data_buffer_changed = true;
				}
			}

			isgv.close();
			osgv.close();

			// Buffer update
			if (nucl_buffer_changed) {
				uint max_nucl = outfile.global_vars["k"] + real_max - 1;
				for (uint i=0 ; i<4 ; i++) {
					delete[] nucleotide_shifts[i];
					nucleotide_shifts[i] = new uint8_t[max_nucl / 4 + 1];
				}
			}
			if (data_buffer_changed) {
				delete[] data;
				data = new uint8_t[outfile.global_vars["data_size"] * real_max];
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
				// Read the block.
				uint nb_kmers = in_section.read_compacted_sequence(nucleotide_shifts[0], data);
				uint seq_nucl = k + nb_kmers - 1;
				uint useful_bytes = seq_nucl%4==0 ? seq_nucl/4 : seq_nucl/4+1;
				// Generate all possible shifts
				for (uint shift=1 ; shift<4 ; shift++) {
					memcpy(nucleotide_shifts[shift], nucleotide_shifts[shift-1], useful_bytes);
					rightshift8(nucleotide_shifts[shift], useful_bytes, 2);
				}
				// Write kmer per kmer
				uint first_nucl = (4 - (seq_nucl % 4)) % 4;
				for (uint kmer_idx=0 ; kmer_idx<nb_kmers ; kmer_idx++) {
					uint shift_idx = (nb_kmers - kmer_idx - 1) % 4;
					uint first_byte = (first_nucl + kmer_idx + shift_idx) / 4;

					out_section.write_compacted_sequence(
							nucleotide_shifts[shift_idx] + first_byte,
							k,
							data + kmer_idx * outfile.global_vars["data_size"]
					);
				}
			}

			in_section.close();
			out_section.close();
		}
		// Revwrite a minimizer block
		else if (section_type == 'm') {
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
				uint64_t nb_kmers = in_section.read_compacted_sequence_without_mini(nucleotide_shifts[0], data, mini_pos);
				uint64_t seq_nucl = k - m + nb_kmers - 1;
				uint useful_bytes = seq_nucl%4==0 ? seq_nucl/4 : seq_nucl/4+1;
				// Generate all possible shifts
				for (uint shift=1 ; shift<4 ; shift++) {
					memcpy(nucleotide_shifts[shift], nucleotide_shifts[shift-1], useful_bytes);
					rightshift8(nucleotide_shifts[shift], useful_bytes, 2);
				}

				// Compute limits of the superkmer (sequence induced by the minimizer)
				int skmer_start = mini_pos - k + m;
				// Write on block per kmer inside of the superkmer
				uint first_nucl = (4 - (seq_nucl % 4)) % 4;
				for (int kmer_idx=max(0, skmer_start) ; kmer_idx<min((int)nb_kmers, (int)mini_pos+1) ; kmer_idx++) {
					uint shift_idx = (nb_kmers - kmer_idx - 1) % 4;
					uint first_byte = (first_nucl + kmer_idx + shift_idx) / 4;

					out_section.write_compacted_sequence_without_mini(
						nucleotide_shifts[shift_idx] + first_byte,
						k - m,
						mini_pos - kmer_idx,
						data + kmer_idx * data_size
					);
				}

				if (skmer_start > 0 or mini_pos > nb_kmers - 1) {
					// Prepare sequence with minimizer
					seq_nucl = k + nb_kmers - 1;
					useful_bytes = seq_nucl%4==0 ? seq_nucl/4 : seq_nucl/4+1;
					in_section.add_minimizer(nb_kmers, nucleotide_shifts[0], mini_pos);
					// Compute all the shifts
					for (uint shift=1 ; shift<4 ; shift++) {
						memcpy(nucleotide_shifts[shift], nucleotide_shifts[shift-1], useful_bytes);
						rightshift8(nucleotide_shifts[shift], useful_bytes, 2);
					}

					first_nucl = (4 - (seq_nucl % 4)) % 4;
					// Save the kmers before superkmer
					for (int kmer_idx=0 ; kmer_idx<skmer_start ; kmer_idx++) {
						// Compute kmer coordinates
						uint shift_idx = (nb_kmers - kmer_idx - 1) % 4;
						uint first_byte = (first_nucl + kmer_idx + shift_idx) / 4;

						// copy the kmer and the data
						uint8_t * kmer = new uint8_t[useful_bytes];
						memcpy(kmer, nucleotide_shifts[shift_idx] + first_byte, useful_bytes);
						uint8_t * data_cpy = new uint8_t[data_size];
						memcpy(data_cpy, data + kmer_idx * data_size, data_size);

						// save kmer and data
						saved_kmers.push_back(kmer);
						saved_data.push_back(data_cpy);
					}

					// Save the kmers after superkmer
					for (uint kmer_idx=mini_pos+1 ; kmer_idx<nb_kmers ; kmer_idx++) {
						// Compute kmer coordinates
						uint shift_idx = (nb_kmers - kmer_idx - 1) % 4;
						uint first_byte = (first_nucl + kmer_idx + shift_idx) / 4;

						// copy the kmer and the data
						uint8_t * kmer = new uint8_t[useful_bytes];
						memcpy(kmer, nucleotide_shifts[shift_idx] + first_byte, useful_bytes);
						uint8_t * data_cpy = new uint8_t[data_size];
						memcpy(data_cpy, data + kmer_idx * data_size, data_size);

						// save kmer and data
						saved_kmers.push_back(kmer);
						saved_data.push_back(data_cpy);
					}
				}

			}

			in_section.close();
			out_section.close();

			if (saved_kmers.size() > 0) {
				// Final step: Write saved kmers into a raw section
				// Open sections
				Section_Raw raw_section(&outfile);

				// Rewrite block per block
				uint64_t k = outfile.global_vars["k"];
				for (uint kmer_idx=0 ; kmer_idx<saved_kmers.size() ; kmer_idx++) {
					raw_section.write_compacted_sequence(
							saved_kmers[kmer_idx],
							k,
							saved_data[kmer_idx]
					);

					delete[] saved_kmers[kmer_idx];
					delete[] saved_data[kmer_idx];
				}

				raw_section.close();
			}
		}

		section_type = infile.read_section_type();
	}

	for (uint i=0 ; i<4 ; i++)
		delete[] nucleotide_shifts[i];
	delete[] data;
	infile.close();
	outfile.close();
}
