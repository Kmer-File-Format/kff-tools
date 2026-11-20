// compute_scores() compared with the original algorithm: update_score() called
// pair by pair in the (first, second) order.

#include <algorithm>
#include <random>
#include <set>
#include <vector>

#include "lest.hpp"
#include "colinear.hpp"
#include "compact.hpp"

using namespace std;

namespace {

/** Random distinct pairs. When dense, few distinct first values: many pairs
 *  share the same first and get equal scores (low complexity sections). */
vector<PairInt> random_pairs(mt19937_64 & rnd, size_t nb, uint64_t nb_first, uint64_t nb_second) {
	set<PairInt> pairs;
	while (pairs.size() < nb)
		pairs.insert(PairInt(rnd() % nb_first, rnd() % nb_second));
	vector<PairInt> result(pairs.begin(), pairs.end());
	shuffle(result.begin(), result.end(), rnd);
	return result;
}

} // namespace

const lest::test module[] = {

	CASE("Colinear::compute_scores gives the same scores and chains as update_score pair by pair") {
		mt19937_64 rnd(40);
		for (int test = 0; test < 600; test++) {
			const size_t nb = 1 + rnd() % (test % 10 == 0 ? 3000 : 200);
			const uint64_t nb_first = 1 + rnd() % (test % 3 == 0 ? 8 : nb + 5);
			const uint64_t nb_second = 1 + rnd() % (test % 4 == 0 ? 8 : nb + 5);
			if (nb > nb_first * nb_second)
				continue;
			vector<PairInt> pairs = random_pairs(rnd, nb, nb_first, nb_second);

			// Original algorithm
			vector<PairInt> ref_pairs = pairs;
			Colinear reference(ref_pairs);
			sort(ref_pairs.begin(), ref_pairs.end(), [](const PairInt & a, const PairInt & b) {
				return a.first != b.first ? a.first < b.first : a.second < b.second;
			});
			for (PairInt & p : ref_pairs)
				reference.update_score(p);

			// Tested function
			vector<PairInt> new_pairs = pairs;
			Colinear colinear(new_pairs);
			colinear.compute_scores(new_pairs);

			// Orders of the pairs: tree order after the constructor, (first, second) after
			// compute_scores, as with the original std::sort calls
			vector<PairInt> tree_order = pairs;
			sort(tree_order.begin(), tree_order.end(), [](const PairInt & a, const PairInt & b) {
				return a.second != b.second ? a.second < b.second : a.first > b.first;
			});
			vector<PairInt> constructor_order = pairs;
			Colinear sorting(constructor_order);
			const bool same_tree_order = constructor_order == tree_order;
			EXPECT( same_tree_order );
			vector<PairInt> processing_order = pairs;
			sort(processing_order.begin(), processing_order.end(), [](const PairInt & a, const PairInt & b) {
				return a.first != b.first ? a.first < b.first : a.second < b.second;
			});
			const bool same_processing_order = new_pairs == processing_order;
			EXPECT( same_processing_order );

			const bool same_scores = colinear.nc_tree == reference.nc_tree;
			const bool same_chains = colinear.chains == reference.chains;
			const bool same_longest = colinear.longest_chain() == reference.longest_chain();
			EXPECT( same_scores );
			EXPECT( same_chains );
			EXPECT( same_longest );

			// Scores and ranks not packed (more than 2^32 pairs)
			vector<PairInt> wide_pairs = pairs;
			Colinear wide(wide_pairs);
			wide.compute_scores(wide_pairs, false);
			const bool same_wide = wide.nc_tree == reference.nc_tree and wide.chains == reference.chains;
			EXPECT( same_wide );
		}
	},
	CASE("colinear_longest_chain gives the longest chain of Colinear") {
		mt19937_64 rnd(55);
		for (int test = 0; test < 3000; test++) {
			const size_t nb = 1 + rnd() % (test % 10 == 0 ? 2000 : 100);
			const uint64_t nb_first = 2 + rnd() % (test % 3 == 0 ? 8 : nb + 5);
			const uint64_t nb_second = 2 + rnd() % (test % 4 == 0 ? 8 : nb + 5);
			if (nb > nb_first * nb_second)
				continue;
			vector<PairInt> pairs = random_pairs(rnd, nb, nb_first, nb_second);
			// Sparse values (comparison sorts instead of counting sorts)
			if (test % 5 == 0)
				for (PairInt & p : pairs) {
					p.first = p.first * 100000007ULL + 3;
					p.second = p.second * 1000000009ULL + 7;
				}

			vector<PairInt> ref_pairs = pairs;
			Colinear reference(ref_pairs);
			reference.compute_scores(ref_pairs);
			const vector<PairInt> expected = reference.longest_chain();

			// Any order of the pairs (sorted by the function)
			vector<PairInt> given = pairs;
			if (test % 2 == 0)
				sort(given.begin(), given.end(), [](const PairInt & a, const PairInt & b) {
					return a.first != b.first ? a.first < b.first : a.second < b.second;
				});
			vector<PairInt> given_wide = given;
			const bool same_chain = colinear_longest_chain(given) == expected;
			const bool same_chain_wide = colinear_longest_chain(given_wide, false) == expected;
			EXPECT( same_chain );
			EXPECT( same_chain_wide );
		}
	},

	CASE("Compact::colinear_chaining gives the longest chain of Colinear") {
		mt19937_64 rnd(53);
		Compact comp;
		for (int test = 0; test < 4000; test++) {
			// Candidates as pair_kmers() gives them: sorted by first, then second. Some of them
			// increasing on both sides, possibly with a few disorders.
			const size_t nb = 1 + rnd() % (test % 10 == 0 ? 500 : 30);
			vector<PairInt> candidates;
			if (test % 2 == 0) {
				uint64_t first = rnd() % 3, second = rnd() % 3;
				for (size_t i = 0; i < nb; i++) {
					candidates.emplace_back(first, second);
					first += 1 + rnd() % 2;
					second += (test % 4 == 0 and rnd() % 50 == 0) ? 0 : 1 + rnd() % 2;
				}
			} else {
				// nb_first * nb_second > nb: enough distinct pairs
				candidates = random_pairs(rnd, nb, 2 + nb / (1 + rnd() % 3), nb + 2);
			}
			sort(candidates.begin(), candidates.end(), [](const PairInt & a, const PairInt & b) {
				return a.first != b.first ? a.first < b.first : a.second < b.second;
			});

			vector<PairInt> pairs = candidates;
			Colinear colinear(pairs);
			colinear.compute_scores(pairs);
			const vector<PairInt> expected = colinear.longest_chain();

			vector<PairInt> given = candidates;
			const bool same_chain = comp.colinear_chaining(given) == expected;
			EXPECT( same_chain );
		}
	},
};


extern lest::tests & specification();

MODULE( specification(), module )
