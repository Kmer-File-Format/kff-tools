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
		}
	},
};

extern lest::tests & specification();

MODULE( specification(), module )
