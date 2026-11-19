// compute_scores() compared with the original algorithm: update_score() called
// pair by pair in the (first, second) order.

#include <algorithm>
#include <random>
#include <set>
#include <vector>

#include "lest.hpp"
#include "colinear.hpp"

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

			const bool same_scores = colinear.nc_tree == reference.nc_tree;
			const bool same_chains = colinear.chains == reference.chains;
			const bool same_longest = colinear.longest_chain() == reference.longest_chain();
			EXPECT( same_scores );
			EXPECT( same_chains );
			EXPECT( same_longest );
		}
	},
};

extern lest::tests & specification();

MODULE( specification(), module )
