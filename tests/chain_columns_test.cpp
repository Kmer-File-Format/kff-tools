// Compact::chain_columns() compared with colinear_chaining(pair_kmers()), on random columns with
// runs of equal kmers (d x e candidate pairs), kmers sharing their suffix and overlapping kmers.

#include <algorithm>
#include <random>
#include <vector>

#include "lest.hpp"
#include "compact.hpp"

using namespace std;

namespace {

/** Nucleotides of a kmer (right aligned in its bytes, random padding bits) */
vector<uint8_t> make_kmer(mt19937_64 & rnd, const vector<uint8_t> & nucleotides) {
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

vector<uint8_t> random_nucleotides(mt19937_64 & rnd, size_t nb) {
	vector<uint8_t> n(nb);
	for (uint8_t & x : n) x = rnd() % 4;
	return n;
}

} // namespace

const lest::test module[] = {

	CASE("Compact::chain_columns gives the chain of colinear_chaining(pair_kmers())") {
		mt19937_64 rnd(56);
		Compact comp;
		for (int test = 0; test < 1500; test++) {
			const size_t nb_nucl = 2 + rnd() % (test % 4 == 0 ? 63 : 20);
			comp.m = 1 + rnd() % 10;
			comp.k = comp.m + nb_nucl;
			// Few distinct nucleotide sequences: many equal kmers and shared suffixes
			const size_t nb_distinct = 1 + rnd() % (test % 3 == 0 ? 4 : 40);
			vector<vector<uint8_t> > pool;
			for (size_t i = 0; i < nb_distinct; i++) {
				pool.push_back(random_nucleotides(rnd, nb_nucl));
				// same suffix, other first nucleotide
				if (rnd() % 3 == 0) {
					vector<uint8_t> variant = pool.back();
					variant[0] = rnd() % 4;
					pool.push_back(variant);
				}
			}

			vector<vector<uint8_t> > storage;
			storage.reserve(4000);
			vector<uint8_t *> column1, column2;
			const size_t max_run = test % 5 == 0 ? 40 : 4;
			const size_t nb_runs1 = rnd() % (test % 10 == 0 ? 60 : 15);
			for (size_t r = 0; r < nb_runs1; r++) {
				const vector<uint8_t> & nucl = pool[rnd() % pool.size()];
				for (size_t x = 1 + rnd() % max_run; x > 0; x--) {
					storage.push_back(make_kmer(rnd, nucl));
					column1.push_back(storage.back().data());
				}
			}
			// Second column: successors of kmers of the pool (shifted), in runs, and random kmers
			const size_t nb_runs2 = rnd() % (test % 10 == 0 ? 60 : 15);
			for (size_t r = 0; r < nb_runs2; r++) {
				vector<uint8_t> nucl;
				if (rnd() % 5 != 0) {
					const vector<uint8_t> & left = pool[rnd() % pool.size()];
					nucl.assign(left.begin() + 1, left.end());
					nucl.push_back(rnd() % 4);
				} else
					nucl = random_nucleotides(rnd, nb_nucl);
				for (size_t x = 1 + rnd() % max_run; x > 0; x--) {
					storage.push_back(make_kmer(rnd, nucl));
					column2.push_back(storage.back().data());
				}
			}
			// Sorted columns (as in compaction: equal kmers next to each other) or not
			if (test % 2 == 0) {
				auto by_bytes = [nb_nucl](const uint8_t * a, const uint8_t * b) {
					// compare the nucleotides only (padding bits ignored)
					const size_t offset = (4 - nb_nucl % 4) % 4;
					for (size_t i = 0; i < nb_nucl; i++) {
						const size_t pos = offset + i;
						const int na = (a[pos / 4] >> (2 * (3 - pos % 4))) & 3;
						const int nb = (b[pos / 4] >> (2 * (3 - pos % 4))) & 3;
						if (na != nb) return na < nb;
					}
					return false;
				};
				stable_sort(column1.begin(), column1.end(), by_bytes);
				stable_sort(column2.begin(), column2.end(), by_bytes);
			} else {
				shuffle(column1.begin(), column1.end(), rnd);
			}

			vector<pair<uint64_t, uint64_t> > candidates = comp.pair_kmers(column1, column2);
			const vector<pair<uint64_t, uint64_t> > expected = comp.colinear_chaining(candidates);
			vector<pair<uint64_t, uint64_t> > chain;
			comp.chain_columns(column1, column2, chain);
			const bool same_chain = chain == expected;
			EXPECT( same_chain );
		}
	},
};

extern lest::tests & specification();

MODULE( specification(), module )
