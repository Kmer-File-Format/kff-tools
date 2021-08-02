#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <queue>

#include "encoding.hpp"
#include "sequences.hpp"
#include "compact.hpp"
#include "merge.hpp"


using namespace std;


Compact::Compact() {
	input_filename = "";
	output_filename = "";
}


Compact::~Compact() {
}


void Compact::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("compact", "Read a kff file and try to compact the kmers from minimizer sections. The available ram must be sufficent to load a complete minimizer section into memory.");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "Input kff file to compact.");
	input_option->required();
	CLI::Option * out_option = subapp->add_option("-o, --outfile", output_filename, "Kff to write (must be different from the input)");
	out_option->required();
}


void Compact::exec() {
	Kff_file infile(input_filename, "r");
	Kff_file outfile(output_filename, "w");

	// TODO: Write encoding
	// TODO: Set the flags
	// TODO: Write header

	// Compute file size
	long current_pos = infile.fs.tellp();
	infile.fs.seekg(0, infile.fs.end);
	long file_size = infile.fs.tellp();
	infile.fs.seekp(current_pos);

	while (infile.fs.tellp() != file_size - 3) {
		char section_type = infile.read_section_type();

		if (section_type == 'v') {
			Section_GV isgv(&infile);
			Section_GV osgv(&outfile);
			for (auto & p : isgv.vars)
				osgv.write_var(p.first, p.second);
			isgv.close();
			osgv.close();
		}
		else if (section_type == 'i') {
			Section_Index si(&infile);
			si.close();
		}
		else if (section_type == 'r') {
			Section_Raw sr(&infile);
			sr.close();
		}
		else if (section_type == 'm') {
			uint k = outfile.global_vars["k"];
			uint m = outfile.global_vars["m"];

			// Rewrite a value section if max is not sufficently large
			if (outfile.global_vars["max"] < (k-m) * 2) {
				unordered_map<string, uint64_t> values(outfile.global_vars);
				Section_GV sgv(&outfile);

				for (auto & p : values)
					if (p.first != "max")
						sgv.write_var(p.first, p.second);
				sgv.write_var("max", (k-m)*2);

				sgv.close();
			}

			// Compact and save the kmers
			Section_Minimizer sm(&infile);
			this->compact_section(sm, outfile);
			sm.close();
		}
	}

	infile.close();
	outfile.close();
}

void Compact::compact_section(Section_Minimizer & ism, Kff_file & outfile) {
	// General variables
	uint k = outfile.global_vars["k"];
	uint m = outfile.global_vars["m"];
	uint data_size = outfile.global_vars["data_size"];
	uint kmer_bytes = (k - m + 3) / 4;

	// Buffers
	uint8_t * seq_buffer = new uint8_t[(2 * k + 3) / 4];
	uint8_t * data_buffer = new uint8_t[2 * k * data_size];
	uint64_t mini_pos;

	// Kmer storing space
	uint kmer_buffer_size = 1 << 15;
	uint next_free = 0;
	uint8_t * kmers = (uint8_t *)malloc(kmer_buffer_size);
	vector<vector<uint8_t *> > kmers_per_index(k-m+1);

	// 1 - Load the input section
	for (uint n=0 ; n<ism.nb_blocks ; n++) {
		// Read sequence
		uint nb_kmers = ism.read_compacted_sequence_without_mini(
			seq_buffer, data_buffer, mini_pos);

		// Add kmer by index
		for (uint kmer_idx=0 ; kmer_idx<nb_kmers ; kmer_idx++) {
			uint kmer_pos = k - m - mini_pos + kmer_idx;

			// Realloc if needed
			if (kmer_buffer_size - next_free < kmer_bytes + data_size) {
				kmer_buffer_size *= 2;
				kmers = (uint8_t *) realloc((void *)kmers, kmer_buffer_size);
			}

			// Copy kmer sequence
			subsequence(seq_buffer, k - m + nb_kmers - 1, kmers + next_free, kmer_idx, kmer_idx + k - m - 1);
			// Copy data array
			memcpy(kmers + next_free + kmer_bytes, data_buffer + kmer_idx * data_size, data_size);
			// Update
			kmers_per_index[kmer_pos].push_back(kmers + next_free);
			next_free += kmer_bytes + data_size;
		}
	}

	// 2 - Compact kmers
	vector<pair<uint8_t *, uint8_t *> > to_compact = this->greedy_assembly(kmers_per_index, k, m);

	Section_Minimizer osm(&outfile);
	osm.write_minimizer(ism.minimizer);
	// TODO: Write compaction from to_compact
	osm.close();

	// Cleaning
	free(kmers);
	delete[] seq_buffer;
	delete[] data_buffer;
}

vector<pair<uint8_t *, uint8_t *> > Compact::greedy_assembly(vector<vector<uint8_t *> > & kmers, const uint k, const uint m) {
	uint nb_kmers = k - m + 1;
	vector<pair<uint8_t *, uint8_t *> > assembly;
	uint8_t encoding[] = {0, 1, 3, 2};
	Stringifyer strif(encoding);

	// Index kmers from the 0th set
	for (uint8_t * kmer : kmers[0])
		assembly.emplace_back(nullptr, kmer);

	for (uint i=0 ; i<nb_kmers-1 ; i++) {
		// Index kmers in ith set
		unordered_map<uint64_t, vector<uint8_t *> > index;
		
		for (uint8_t * kmer : kmers[i]) {
			// Get the suffix
			uint64_t val = subseq_to_uint(kmer, nb_kmers-1, 1, nb_kmers-1);
			// Add a new vector for this value
			if (index.find(val) == index.end())
				index[val] = vector<uint8_t *>();
			// Add the kmer to the value list
			index[val].push_back(kmer);
		}

		// link kmers from (i+1)th set to ith kmers.
		for (uint8_t * kmer : kmers[i+1]) {
			uint64_t val = subseq_to_uint(kmer, nb_kmers-1, 0, nb_kmers-2);

			if (index.find(val) == index.end()) {
				// No kmer available for matching
				assembly.emplace_back(nullptr, kmer);
			} else {
				bool chaining_found = false;
				uint candidate_pos = 0;
				// verify complete matching for candidates kmers
				for (uint8_t * candidate : index[val]) {
					// If the kmers can be assembled
					if (sequence_compare(
								kmer, nb_kmers-1, 0, nb_kmers-2,
								candidate, nb_kmers-1, 1, nb_kmers-1
							) == 0) {
						// Update status
						chaining_found = true;
						assembly.emplace_back(kmer, candidate);

						cout << strif.translate(kmer, nb_kmers-1) << " -> " << strif.translate(candidate, nb_kmers-1) << endl;

						// remove candidate from list
						index[val].erase(index[val].begin()+candidate_pos);
						// Quit candidate searching
						break;
					}

					candidate_pos += 1;
				}
				// If no assembly possible, create a new superkmer
				if (not chaining_found)
					assembly.emplace_back(nullptr, kmer);
			}
		}
	}

	return assembly;
}

