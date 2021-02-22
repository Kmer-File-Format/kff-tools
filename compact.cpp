#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <algorithm>

#include "sequences.hpp"
#include "compact.hpp"


using namespace std;


Compact::Compact() {
	input_filename = "";
	output_filename = "";

	load_mem_size = 1;
	loading_memory = new uint8_t[load_mem_size];
}

void Compact::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("compact", "Read a kff file and try to compact the kmers from minimizer sections. The available ram must be sufficent to load a complete minimizer section into memory");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "Input kff file to compact.");
	input_option->required();

	CLI::Option * out_option = subapp->add_option("-o, --outfile", output_filename, "Kff to write (must be different from the input)");
	out_option->required();
}

void Compact::compact(string input, string output) {
	// Read the encoding of the first file and push it as outcoding
	Kff_file infile(input, "r");

	// Write header of the output
	Kff_file outfile(output, "w");
	outfile.write_encoding(infile.encoding);
	// Set metadata
	uint8_t * metadata = new uint8_t[infile.metadata_size];
	infile.read_metadata(metadata);
	outfile.write_metadata(infile.metadata_size, metadata);
	delete[] metadata;

	long buffer_size = 1048576; // 1 MB
	char buffer[1048576];

	// Read section by section
	char section_type = infile.read_section_type();
	while(not infile.fs.eof()) {
		switch (section_type) {
			// Write the variables that change from previous sections (possibly sections from other input files)
			case 'v':
			{
				// Read the values
				Section_GV in_sgv(&infile);
				Section_GV out_sgv(&outfile);

				// Verify the presence and value of each variable in output
				for (auto& tuple : in_sgv.vars) {
					out_sgv.write_var(tuple.first, tuple.second);
				}

				in_sgv.close();
				out_sgv.close();
			}
			break;

			// copy the sequence section from input to output
			case 'r':
			{
				// Analyse the section size to prepare copy
				cerr << "WARNING: Raw sections are not compacted !" << endl;
				auto begin_byte = infile.fs.tellp();
				if (not infile.jump_next_section()) {
					cerr << "Error inside of the input file." << endl;
					cerr << "Impossible to jump over the raw section at byte " << begin_byte << endl;
					exit(1);
				}
				auto end_byte = infile.fs.tellp();
				long size = end_byte - begin_byte;
				infile.fs.seekp(begin_byte);

				// copy from input to output
				while (size > 0) {
					size_t to_copy = size > buffer_size ? buffer_size : size;

					infile.fs.read(buffer, to_copy);
					outfile.fs.write(buffer, to_copy);

					size -= to_copy;
				}
			}
			case 'm':
			{
				uint m = infile.global_vars["m"];
				// Open the input minimizer section
				Section_Minimizer sm(&infile);
				uint8_t * mini = new uint8_t[(m+3)/4];
				memcpy(mini, sm.minimizer, (m+3)/4);
				// Load all the blocks
				this->kmer_nbs.clear();
				this->mini_pos.clear();
				this->loadSectionBlocks(sm, infile);
				// Compute paths from left to right for each available sequence
				vector<vector<uint> > paths = this->link_kmers(sm.nb_blocks, infile);
				// Assemble and write the paths
				// close the section
				sm.close();

				delete[] mini;
			}
			break;

			default:
				cerr << "Unknown section type " << section_type << " in file " << input_filename << endl;
				exit(2);
		}

		// Prepare next section
		section_type = infile.read_section_type();
	}

	infile.close();
	outfile.close();
}

void Compact::loadSectionBlocks(Section_Minimizer & sm, Kff_file & infile) {
	cout << "--- load ---" << endl;
	uint k = infile.global_vars["k"];
	uint m = infile.global_vars["m"];
	uint max_kmers = infile.global_vars["max"];
	uint data_size = infile.global_vars["data_size"];
	
	uint skmer_nucl_bytes = (max_kmers + k - 1 - m + 3) / 4;
	uint skmer_data_bytes = max_kmers * data_size;
	uint max_block_size = skmer_nucl_bytes + skmer_data_bytes;

	// Realloc
	uint max_size = max_block_size * sm.nb_blocks;
	if (load_mem_size <  max_size) {
		delete[] loading_memory;
		loading_memory = new uint8_t[max_size];
		load_mem_size = max_size;
	}

	// Read all the sequences
  for (uint64_t i=0 ; i<sm.nb_blocks ; i++) {
  	// Save sequences and data
    uint64_t minimizer_position;
    uint64_t nb_kmers = sm.read_compacted_sequence_without_mini(
    		this->loading_memory + i * max_block_size,
    		this->loading_memory + i * max_block_size + skmer_nucl_bytes,
    		minimizer_position);

    // Save values
    this->kmer_nbs.push_back(nb_kmers);
    this->mini_pos.push_back(minimizer_position);

    // Tmp print
    cout << nb_kmers << " " << minimizer_position << endl;
    cout << (uint)(this->loading_memory + i * max_block_size)[0] << " ";
    cout << (uint)(this->loading_memory + i * max_block_size)[1] << " ";
    cout << endl;
    cout << "minimizer position: " << minimizer_position << endl;
  }
}

vector<vector<uint> > Compact::link_kmers(uint nb_kmers, Kff_file & infile) {
	cout << "--- Create links between sequences ---" << endl;
	// Usefull variables
	uint k = infile.global_vars["k"];
	uint m = infile.global_vars["m"];
	uint max_kmers = infile.global_vars["max"];
	uint data_size = infile.global_vars["data_size"];
	cout << "k " << k << " m " << m << endl;
	
	uint skmer_nucl_bytes = (max_kmers + k - 1 - m + 3) / 4;
	uint skmer_data_bytes = max_kmers * data_size;
	uint max_block_size = skmer_nucl_bytes + skmer_data_bytes;

	uint8_t * sub_seq = new uint8_t[k/4 + 1];

	// map <prefix_size, map<64bits_prefix_value, vector<sequence_idx> > >
	unordered_map<uint, unordered_map<uint64_t, vector<uint> > > prefix_bins;
	vector<uint> present_bins;

	// Distribute kmers into positional bins.
	// Bin x means x nucleotides inside of the prefix
	uint kmer_idx = 0;
	cout << "prefix computation " << endl;
	for (uint & mini_idx : this->mini_pos) {
		cout << "mini idx " << mini_idx << endl;
		// Create absent bin
		if (prefix_bins.find(mini_idx) == prefix_bins.end()) {
			cout << "new bin " << mini_idx << endl;
			prefix_bins[mini_idx] = unordered_map<uint64_t, vector<uint> >();
			present_bins.push_back(mini_idx);
		}

		// Compute the prefix
		uint8_t * seq = loading_memory + kmer_idx * max_block_size;
		subsequence(seq,
				k - 1 + this->kmer_nbs[kmer_idx] - m,
				sub_seq, 0, (k-1)-m-1);
		cout << (uint)sub_seq[0] << " " << (uint)sub_seq[1] << endl;
		cout << "pref param " << (k - 1 - m) << endl;
		// Prefix to int value
		uint64_t prefix_val = seq_to_uint(sub_seq, k - m - 1);
		cout << "pref val " << prefix_val << endl;

		// Insert the prefix
		prefix_bins[mini_idx][prefix_val].push_back(kmer_idx);
		kmer_idx += 1;
	}
	cout << endl;
	cout << "before sort " << present_bins.size() << endl;
	for (auto & i : present_bins) cout << i << " ";
	cout << endl;
	sort(present_bins.begin(), present_bins.end());
	cout << "before reverse " << present_bins.size() << endl;
	for (auto & i : present_bins) cout << i << " ";
	cout << endl;
	reverse(present_bins.begin(), present_bins.end());
	cout << "after all " << present_bins.size() << endl;
	for (auto & i : present_bins) cout << i << " ";
	cout << endl;

	// Create paths of overlaping sequences
	// Each sequence is present in only one path
	vector<vector<uint> > paths;
	vector<uint> current_path;

	while(present_bins.size() > 0) {
		// Get a start point
		unordered_map<uint64_t, vector<uint> > & first_idx_bin = prefix_bins[present_bins[0]];
		uint64_t first_prefix = first_idx_bin.begin()->first;
		vector<uint> & first_kmers_idx = first_idx_bin.begin()->second;

		cout << "first [" << first_prefix << "] " << first_kmers_idx.size() << endl;

		// Start a path
		current_path.push_back(first_kmers_idx[0]);
		cout << "first in path " << current_path[0] << endl;
		first_kmers_idx.erase(first_kmers_idx.begin());
		if (first_kmers_idx.size() == 0) {
			first_idx_bin.erase(first_prefix);
			if (first_idx_bin.size() == 0) {
				prefix_bins.erase(present_bins[0]);
				present_bins.erase(present_bins.begin());
			}
		}

		// Loop over sequences to concatenate
		while (true) {
			uint last_idx = current_path.back();
			uint last_kmer_nb = this->kmer_nbs[last_idx];
			uint last_mini_idx = this->mini_pos[last_idx];
			uint8_t * last_seq = loading_memory + last_idx * max_block_size;
			uint last_seq_size = k - 1 + last_kmer_nb - m;

			// Compute the next bin to look for
			int suff_idx = (int)last_mini_idx - (int)last_kmer_nb;
			cout << "suff idx " << suff_idx << endl;
			if (suff_idx < 0 or prefix_bins.find((uint)suff_idx) == prefix_bins.end()) {
				cerr << "No sequence with the good minimizer idx" << endl;
				break;
			}

			// Extract the (k-1)-suffix of the last element of the path
			cout << "suffix " << (last_seq_size - 1) << endl;
			// TODO : HERE
			subsequence(last_seq, last_seq_size, sub_seq, suff_idx, last_seq_size-1);
			uint64_t suff_val = seq_to_uint(sub_seq, last_seq_size - 1 - suff_idx);

			// Look for the elements with the good (k-1)-prefix
			if (prefix_bins[suff_idx].find(suff_val) == prefix_bins[suff_idx].end()) {
				cerr << "No sequence with the good prefix" << endl;
				break;
			}
			// Iterate over candidates
			bool valid_candidate = false;
			uint compaction_idx;
			vector<uint> & candidates = prefix_bins[suff_idx][suff_val];
			if (last_seq_size - 1 - suff_idx <= 32) {
				valid_candidate = true;
				compaction_idx = 0;
				current_path.push_back(candidates[0]);
			} else {
				cerr << "TODO: implement verification for suffix larger than 32 nucleotides" << endl;
				exit(1);
			}

			// Remove concatenated from candidates
			if (valid_candidate) {
				candidates.erase(candidates.begin() + compaction_idx);
				// Remove the (k-1) mer from possibilities if no remaining prefix
				if (candidates.size() == 0) {
					prefix_bins[suff_idx].erase(suff_val);

					// Remove the minimizer indice for kmer if there is no available kmer in this bin.
					if (prefix_bins[suff_idx].size() == 0) {
						prefix_bins.erase(suff_idx);
						remove(present_bins.begin(), present_bins.end(), suff_idx);
					}
				}
			}
		}

		// Add the path
		paths.push_back(current_path);
		current_path = vector<uint>();
	}

	delete[] sub_seq;

	cout << "paths " << paths.size() << " " << paths[0].size() << endl << endl;
	return paths;
}

void Compact::exec() {
	this->compact(input_filename, output_filename);
}
