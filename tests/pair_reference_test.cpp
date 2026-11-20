// pair_kmers() compared with the version of kff-tools 3124a11 (copied below), on random columns
// with many overlapping kmers and hash collisions.

#include <algorithm>
#include <random>
#include <vector>

#include "lest.hpp"
#include "compact.hpp"
#include "sequences.hpp"

using namespace std;

namespace {

vector<pair<uint64_t, uint64_t> > reference_pair_kmers(const vector<uint8_t *> & column1, const vector<uint8_t *> & column2, uint k, uint m) {
	const uint nb_nucl = k - m;
    vector<pair<uint64_t, uint64_t> > pairs;
	pairs.reserve(max(column1.size(), column2.size()));

	// Index the second column by their prefix hash: (hash, kmer position) sorted by hash, then
	// by position (the candidates of a hash are tested in the column order)
	vector<pair<uint64_t, uint64_t> > index;
	index.reserve(column2.size());
	for (uint64_t idx=0 ; idx<column2.size() ; idx++) {
		// Get the hash corresponding to the k-m-1 prefix
		index.emplace_back(subseq_to_uint(column2[idx], nb_nucl, 0, nb_nucl-2), idx);
	}
	// Explicit comparison: this file specialises std::operator< for pairs (second first)
	auto by_hash = [](const pair<uint64_t, uint64_t> & a, const pair<uint64_t, uint64_t> & b) {
		return a.first != b.first ? a.first < b.first : a.second < b.second;
	};
	sort(index.begin(), index.end(), by_hash);


	// Looks for suffix matches of the first column
	for (uint64_t idx=0 ; idx<column1.size() ; idx++) {
		uint8_t * kmer = column1[idx];

		// Get the hash corresponding to the k-m-1 suffix
		uint64_t hash = subseq_to_uint(kmer, nb_nucl, 1, nb_nucl-1);

		// Test each candidate with the same hash
		auto candidate = lower_bound(index.begin(), index.end(), make_pair(hash, uint64_t(0)), by_hash);
		for ( ; candidate != index.end() and candidate->first == hash ; candidate++) {
			uint64_t candidate_vect_idx = candidate->second;
            if (sequence_compare(
						column2[candidate_vect_idx], nb_nucl, 0, nb_nucl-2,
						kmer, nb_nucl, 1, nb_nucl-1
					) == 0) {
				pairs.emplace_back(idx, candidate_vect_idx);
			}
		}
	}

	return pairs;
}


/** Nucleotide i of a right aligned kmer */
uint8_t get_nucleotide(const vector<uint8_t> & kmer, uint nb_nucl, uint i) {
	const uint pos = (4 - nb_nucl % 4) % 4 + i;
	return (kmer[pos / 4] >> (2 * (3 - pos % 4))) & 3;
}

void set_nucleotide(vector<uint8_t> & kmer, uint nb_nucl, uint i, uint8_t nucl) {
	const uint pos = (4 - nb_nucl % 4) % 4 + i;
	kmer[pos / 4] &= ~(3 << (2 * (3 - pos % 4)));
	kmer[pos / 4] |= nucl << (2 * (3 - pos % 4));
}

} // namespace

const lest::test module[] = {

	CASE("Compact::pair_kmers gives the pairs of the original version") {
		mt19937_64 rnd(54);
		for (int test = 0; test < 2000; test++) {
			const uint nb_nucl = 1 + rnd() % (test % 5 == 0 ? 70 : 36);
			Compact comp;
			comp.m = 1 + rnd() % 10;
			comp.k = comp.m + nb_nucl;
			const size_t nb_bytes = (nb_nucl + 3) / 4;

			// Kmers of the first column (random padding bits), few distinct values in some tests.
			// 8 more bytes: the original version reads after the kmer when k - m = 1 (in the
			// kmer buffer, the data and the minimizer position follow the kmer).
			const size_t nb1 = rnd() % (test % 10 == 0 ? 400 : 30);
			const size_t nb_values = test % 3 == 0 ? 1 + rnd() % 5 : 1000000;
			vector<vector<uint8_t> > storage;
			storage.reserve(2 * nb1 + 64);
			vector<uint8_t *> column1, column2;
			for (size_t i = 0; i < nb1; i++) {
				storage.emplace_back(nb_bytes + 8);
				for (uint8_t & b : storage.back()) b = rnd();
				if (nb_values < 1000000 and i >= nb_values)
					storage.back() = storage[rnd() % nb_values];
				column1.push_back(storage.back().data());
			}
			// Second column: shifted kmers of the first one (they overlap), some of them with a
			// difference far from the end (same 32 last nucleotides), and random kmers
			const size_t nb2 = rnd() % (nb1 + 20);
			for (size_t j = 0; j < nb2; j++) {
				storage.emplace_back(nb_bytes + 8);
				vector<uint8_t> & kmer = storage.back();
				for (uint8_t & b : kmer) b = rnd();
				if (nb1 > 0 and rnd() % 4 != 0) {
					const vector<uint8_t> & left = storage[rnd() % nb1];
					for (uint i = 0; i + 1 < nb_nucl; i++)
						set_nucleotide(kmer, nb_nucl, i, get_nucleotide(left, nb_nucl, i + 1));
					if (nb_nucl > 34 and rnd() % 3 == 0)
						set_nucleotide(kmer, nb_nucl, rnd() % (nb_nucl - 33), rnd() % 4);
				}
				column2.push_back(kmer.data());
			}

			const vector<pair<uint64_t, uint64_t> > expected = reference_pair_kmers(column1, column2, comp.k, comp.m);
			const bool same_pairs = comp.pair_kmers(column1, column2) == expected;
			EXPECT( same_pairs );
		}
	},
};

extern lest::tests & specification();

MODULE( specification(), module )
