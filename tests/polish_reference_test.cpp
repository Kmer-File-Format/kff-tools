// polish_sort() compared with the original algorithm (copied below), on random kmer matrices
// and random colinear links between consecutive columns.

#include <algorithm>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "lest.hpp"
#include "compact.hpp"
#include "skmers.hpp"

using namespace std;

namespace {

/** Original polish_sort() (kff-tools 83fff4d), num_kmer_nucl = k - m */
vector<vector<uint8_t *> > reference_polish_sort(const vector<vector<uint8_t *> > & matrix , const vector<vector<pair<uint64_t, uint64_t> > > & pairs, size_t num_kmer_nucl) {
	// Index the kmer columns
	unordered_map<uint8_t *, size_t> columns;
	for (size_t col=0 ; col<matrix.size() ; col++)
		for (uint8_t * kmer : matrix[col])
			columns[kmer] = col;

	// Unify kmers into skmers
	unordered_map<uint8_t *, vector<uint8_t *> *> rev_skmers_index;
	unordered_set<vector<uint8_t *> *> skmers;
	for (size_t col_idx = 0 ; col_idx<pairs.size() ; col_idx++) {
		// Get the current column
		const vector<pair<uint64_t, uint64_t> > & column = pairs[col_idx];

		// for each new link, try to extend previous skmers
		for (const pair<uint64_t, uint64_t> & link : column) {
			uint8_t * lkmer = matrix[col_idx][link.first];
			uint8_t * rkmer = matrix[col_idx+1][link.second];

			// Previous kmer not found
			if (rev_skmers_index.find(lkmer) == rev_skmers_index.end()) {
				// Create a new superkmer
				vector<uint8_t *> * new_sk = new vector<uint8_t *>();
				skmers.insert(new_sk);
				new_sk->push_back(lkmer); new_sk->push_back(rkmer);
				// Insert it in the rev index
				rev_skmers_index[rkmer] = new_sk;
			}
			// Previous kmer found: extend it
			else {
				// Get the skmer
				vector<uint8_t *> * sk = rev_skmers_index[lkmer];
				// Extend it
				sk->push_back(rkmer);
				// change the reference of the sk in the index
				rev_skmers_index.erase(lkmer);
				rev_skmers_index[rkmer] = sk;
			}
		}
	}


	// Index skmers
	unordered_map<uint8_t *, vector<uint8_t *> *> sk_index;
	for(auto & kv: rev_skmers_index) {
		vector<uint8_t *> * sk = kv.second;
		// Index by kmer
		for (uint8_t * kmer: *sk)
			sk_index[kmer] = sk;
	}

	// Index lonely kmers
	for(size_t col=0 ; col<matrix.size() ; col++) {
		for (uint8_t * kmer : matrix[col]) {
			if (sk_index.find(kmer) == sk_index.end()) {
				sk_index[kmer] = new vector<uint8_t *>{kmer};
				skmers.insert(sk_index[kmer]);
			}
		}
	}

	// Init the kmer per superkmer counters with the first row
	unordered_map<vector<uint8_t *> *, uint64_t> sk_counts;
	for (const vector<uint8_t *> & col: matrix) {
		if (col.size() == 0) continue;

		auto & sk = sk_index[col[0]];
		if (sk_counts.find(sk) == sk_counts.end())
			sk_counts[sk] = 0;
		sk_counts[sk] += 1;
	}

	bool matrix_consumed = false;
	// Row indexes for each or the columns in the kmer matrix.
	vector<size_t> current_kmers(matrix.size(), 0);
	// Memory to strore interleaves
	vector<uint8_t *> memory;
	// Interleaves of current pointed kmers
	vector<interleved_t> current_interleaves;

	// Fill the first interleaves
	for (size_t i=0 ; i<matrix.size() ; i++) {
		memory.push_back(new uint8_t[(num_kmer_nucl + 3) / 4]);
		current_interleaves.push_back((interleved_t){0,0,0});
		if (matrix[i].size() > 0) {
			current_interleaves[i] = interleaved(matrix[i][0], memory[i],
												num_kmer_nucl, num_kmer_nucl - i);
		}
	}

	vector<vector<uint8_t *> > result;

	// Add iteratively the skmers into the sorted list
	while(not matrix_consumed) {
		vector<interleved_t> candidates;
		// Construct the list of compatible kmers
		for (size_t i=0 ; i<matrix.size() ; i++) {
			// Add the interleaved kmer with also a position interleaved.
			// The one in the middle is first in the list and then alternate left-right
			size_t position = matrix.size()/2;
			if (i % 2 == 0)
				position += i/2;
			else
				position -= 1 + i/2;

			if (current_kmers[position] < matrix[position].size()) {
				// Get the related skmer
				uint8_t * kmer = matrix[position][current_kmers[position]];
				vector<uint8_t *> * sk = sk_index[kmer];

				// Add the skmer only if all the related kmers are present
				if ((sk_counts.find(sk) != sk_counts.end()) and (sk_counts[sk] == sk->size())) {
					candidates.push_back(current_interleaves[position]);
				}
			}
		}

		// Select the min kmer and add the related skmer to the output
		interleved_t selected = min_interleaved(candidates.begin(), candidates.end());
		candidates = vector<interleved_t>();
		size_t col = selected.suf_size;
		uint8_t * selected_kmer = matrix[col][current_kmers[col]];
		vector<uint8_t *> * selected_sk = sk_index[selected_kmer];

		result.push_back(*selected_sk);

		// remove selected kmers from the counts
		sk_counts.erase(selected_sk);

		// Jump over selected kmers in the matrix
		for (uint8_t * kmer : *selected_sk) {
			col = columns[kmer];
			current_kmers[col] += 1;

			size_t row = current_kmers[col];
			if (row < matrix[col].size()) {
				// Update the interleaved value
				current_interleaves[col] = interleaved(matrix[col][row], memory[col],
													num_kmer_nucl, num_kmer_nucl - col);
				// Update the sk counts
				auto sk = sk_index[matrix[col][row]];
				if (sk_counts.find(sk) == sk_counts.end())
					sk_counts[sk] = 0;
				sk_counts[sk] += 1;
			}
		}

		// Is it over
		matrix_consumed = true;
		for (size_t i=0 ; i<matrix.size() ; i++)
			if (current_kmers[i] < matrix[i].size()) {
				matrix_consumed = false;
				break;
			}
	}

	// free memory
	for (uint8_t * mem : memory)
		delete[](mem);
	for (vector<uint8_t *> * sk : skmers)
		delete(sk);

	return result;
}


/** Random strictly increasing links between 2 columns, as given by a colinear chain */
vector<pair<uint64_t, uint64_t> > random_links(mt19937_64 & rnd, uint64_t size1, uint64_t size2) {
	vector<pair<uint64_t, uint64_t> > links;
	if (size1 == 0 or size2 == 0)
		return links;
	const uint64_t density = rnd() % 4;
	uint64_t i = rnd() % 2, j = rnd() % 2;
	while (i < size1 and j < size2) {
		if (density == 0 or rnd() % 4 < density)
			links.emplace_back(i, j);
		i += 1 + (rnd() % 3 == 0 ? rnd() % 3 : 0);
		j += 1 + (rnd() % 3 == 0 ? rnd() % 3 : 0);
	}
	return links;
}

} // namespace

const lest::test module[] = {

	CASE("Compact::polish_sort gives the skmers of the original algorithm") {
		mt19937_64 rnd(51);
		for (int test = 0; test < 3000; test++) {
			const uint nb_nucl = 1 + rnd() % (test % 7 == 0 ? 70 : 20);
			Compact comp;
			comp.m = 1 + rnd() % 10;
			comp.k = comp.m + nb_nucl;
			const size_t nb_bytes = (nb_nucl + 3) / 4;

			// Random kmers (padding bits included), few distinct values in some tests (ties)
			const uint64_t nb_values = test % 3 == 0 ? 1 + rnd() % 4 : 1000000;
			vector<vector<uint8_t> > values(min<uint64_t>(nb_values, 64));
			for (vector<uint8_t> & v : values) {
				v.resize(nb_bytes);
				for (uint8_t & b : v) b = rnd();
			}
			vector<vector<uint8_t> > storage;
			vector<vector<uint8_t *> > matrix(nb_nucl + 1);
			vector<size_t> sizes(nb_nucl + 1);
			size_t total = 0;
			for (size_t & s : sizes) {
				s = rnd() % 5 == 0 ? 0 : rnd() % (test % 10 == 0 ? 300 : 12);
				total += s;
			}
			// At least 1 kmer: the original algorithm crashes on an empty matrix (compaction
			// never gives it an empty section)
			if (total == 0) {
				sizes[rnd() % sizes.size()] = 1;
				total = 1;
			}
			storage.reserve(total);
			for (size_t col = 0; col <= nb_nucl; col++)
				for (size_t row = 0; row < sizes[col]; row++) {
					if (nb_values < 1000000)
						storage.push_back(values[rnd() % values.size()]);
					else {
						storage.emplace_back(nb_bytes);
						for (uint8_t & b : storage.back()) b = rnd();
					}
					matrix[col].push_back(storage.back().data());
				}

			vector<vector<pair<uint64_t, uint64_t> > > links;
			for (size_t col = 0; col < nb_nucl; col++)
				links.push_back(random_links(rnd, matrix[col].size(), matrix[col + 1].size()));

			const vector<vector<uint8_t *> > expected = reference_polish_sort(matrix, links, nb_nucl);
			const vector<vector<uint8_t *> > skmers = comp.polish_sort(matrix, links);
			const bool same_skmers = skmers == expected;
			EXPECT( same_skmers );
		}
	},
};

extern lest::tests & specification();

MODULE( specification(), module )
