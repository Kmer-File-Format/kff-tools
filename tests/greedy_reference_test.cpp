// greedy_assembly() compared with the version of kff-tools 5cb3b8f (copied below), on random
// columns with runs of equal kmers, kmers sharing their suffix and overlapping kmers.

#include <algorithm>
#include <random>
#include <vector>

#include "lest.hpp"
#include "compact.hpp"
#include "sequences.hpp"

using namespace std;

namespace {

vector<pair<uint8_t *, uint8_t *> > reference_greedy_assembly(vector<vector<uint8_t *> > & kmers, uint k, uint m) {
	uint nb_nucl = k - m;
	vector<pair<uint8_t *, uint8_t *> > assembly;

	// Index kmers from the 0th set
	for (uint8_t * kmer : kmers[0]) {
		assembly.emplace_back(nullptr, kmer);
	}

	// Kmers of the ith set sorted by suffix value. The sort is stable: kmers with the same
	// value stay in their order of the set (the candidates are tested in this order).
	struct Candidate {
		uint64_t value;
		uint8_t * kmer;
		bool used;
	};
	vector<Candidate> index;

	for (uint i=0 ; i<nb_nucl ; i++) {
		// Index kmers in ith set
		index.clear();
		for (uint8_t * kmer : kmers[i]) {
			// Get the suffix
			uint64_t val = subseq_to_uint(kmer, nb_nucl, 1, nb_nucl-1);
			index.push_back({val, kmer, false});
		}
		stable_sort(index.begin(), index.end(), [](const Candidate & a, const Candidate & b) {
			return a.value < b.value;
		});

		// link kmers from (i+1)th set to ith kmers.
		for (uint8_t * kmer : kmers[i+1]) {
			uint64_t val = subseq_to_uint(kmer, nb_nucl, 0, nb_nucl-2);
			auto candidate = lower_bound(index.begin(), index.end(), val, [](const Candidate & c, uint64_t v) {
				return c.value < v;
			});

			bool chaining_found = false;
			// verify complete matching for the remaining candidates kmers
			for ( ; candidate != index.end() and candidate->value == val ; candidate++) {
				if (candidate->used)
					continue;
				// If the kmers can be assembled
				if (sequence_compare(
							kmer, nb_nucl, 0, nb_nucl-2,
							candidate->kmer, nb_nucl, 1, nb_nucl-1
						) == 0) {
					// Update status
					chaining_found = true;
					assembly.emplace_back(candidate->kmer, kmer);
					// remove candidate from list
					candidate->used = true;
					// Quit candidate searching
					break;
				}
			}
			// If no assembly possible, create a new superkmer
			if (not chaining_found) {
				assembly.emplace_back(nullptr, kmer);
			}
		}
	}

	// Index last kmers without compaction
	int assembly_idx = assembly.size()-1;
	for (auto it=kmers[nb_nucl].end() ; it>kmers[nb_nucl].begin() ; it--) {
		uint8_t * kmer = *(it-1);

		if (assembly_idx < 0 or kmer != assembly[assembly_idx].second) {
			assembly.emplace_back(nullptr, kmer);
		} else {
			assembly_idx--;
		}
	}

	return assembly;
}


/** Nucleotides of a kmer (right aligned in its bytes, random padding bits, 8 more bytes) */
vector<uint8_t> greedy_kmer(mt19937_64 & rnd, const vector<uint8_t> & nucleotides) {
	const size_t nb_nucl = nucleotides.size();
	const size_t offset = (4 - nb_nucl % 4) % 4;
	vector<uint8_t> kmer((nb_nucl + 3) / 4 + 8);
	for (uint8_t & b : kmer) b = rnd();
	for (size_t i = 0; i < nb_nucl; i++) {
		const size_t pos = offset + i;
		kmer[pos / 4] &= ~(3 << (2 * (3 - pos % 4)));
		kmer[pos / 4] |= nucleotides[i] << (2 * (3 - pos % 4));
	}
	return kmer;
}

} // namespace

const lest::test module[] = {

	CASE("Compact::greedy_assembly gives the links of the original version") {
		mt19937_64 rnd(58);
		Compact comp;
		for (int test = 0; test < 1500; test++) {
			const size_t nb_nucl = 1 + rnd() % (test % 4 == 0 ? 70 : 20);
			comp.m = 1 + rnd() % 10;
			comp.k = comp.m + nb_nucl;
			// Few distinct sequences: equal kmers, shared suffixes, overlaps between columns
			vector<vector<uint8_t> > pool;
			for (size_t i = 1 + rnd() % (test % 3 == 0 ? 4 : 30); i > 0; i--) {
				vector<uint8_t> nucl(nb_nucl);
				for (uint8_t & x : nucl) x = rnd() % 4;
				pool.push_back(nucl);
				if (rnd() % 3 == 0) {
					nucl[0] = rnd() % 4;
					pool.push_back(nucl);
				}
				// a difference far from the end (same 32 last nucleotides)
				if (nb_nucl > 34 and rnd() % 3 == 0) {
					nucl[rnd() % (nb_nucl - 33)] ^= 1;
					pool.push_back(nucl);
				}
			}
			vector<vector<uint8_t> > storage;
			storage.reserve(20000);
			vector<vector<uint8_t *> > matrix(nb_nucl + 1);
			const size_t max_run = test % 5 == 0 ? 20 : 3;
			for (size_t col = 0; col <= nb_nucl; col++) {
				for (size_t r = rnd() % (test % 10 == 0 ? 40 : 8); r > 0; r--) {
					vector<uint8_t> nucl;
					if (col > 0 and rnd() % 4 != 0) {
						// successor of a kmer of the pool (or of the previous column)
						const vector<uint8_t> & left = pool[rnd() % pool.size()];
						nucl.assign(left.begin() + 1, left.end());
						nucl.push_back(rnd() % 4);
						if (rnd() % 2 == 0)
							pool.push_back(nucl);
					} else
						nucl = pool[rnd() % pool.size()];
					for (size_t x = 1 + rnd() % max_run; x > 0 and storage.size() < 20000; x--) {
						storage.push_back(greedy_kmer(rnd, nucl));
						matrix[col].push_back(storage.back().data());
					}
				}
				if (test % 2 == 1)
					shuffle(matrix[col].begin(), matrix[col].end(), rnd);
			}

			vector<vector<uint8_t *> > reference_matrix = matrix;
			const vector<pair<uint8_t *, uint8_t *> > expected = reference_greedy_assembly(reference_matrix, comp.k, comp.m);
			const bool same_links = comp.greedy_assembly(matrix) == expected;
			EXPECT( same_links );
		}
	},
};

extern lest::tests & specification();

MODULE( specification(), module )
