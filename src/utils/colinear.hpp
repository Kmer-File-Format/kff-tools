#include <vector>
#include <unordered_map>
#include <cstdint>


#ifndef COLINEAR_H
#define COLINEAR_H

// Redefine pair operators
typedef  std::pair<uint64_t, uint64_t> PairInt;
typedef std::pair<PairInt, uint64_t> Node;

struct hash_PairInt {
    // The 2 values are mixed (splitmix64 finalizer): with first ^ second, the pairs of close
    // indices had very few distinct hashes and the lookups became linear.
    size_t operator()(const PairInt& p) const
    {
        uint64_t x = p.first * 0x9E3779B97F4A7C15ULL + p.second;
        x ^= x >> 30;
        x *= 0xBF58476D1CE4E5B9ULL;
        x ^= x >> 27;
        x *= 0x94D049BB133111EBULL;
        x ^= x >> 31;
        return x;
    }
};


/** Longest chain of pairs increasing on both values, the one that Colinear::longest_chain() gives
 * after Colinear::compute_scores(), without the tree of Colinear.
 * @param pairs Distinct pairs, sorted in (first, second) order by the function.
 * @param packed Use scores and ranks packed in 64 bits when there are less than 2^32 pairs
 * (false: always 2 values, for the tests).
 * @return The chain, in increasing order.
 **/
std::vector<PairInt> colinear_longest_chain(std::vector<PairInt> & pairs, bool packed=true);


class Colinear
{
private:
	std::unordered_map<PairInt, uint64_t, hash_PairInt> tree_idxs;
public:
	// Non coliding tree
	std::vector<Node> nc_tree;
	// Previous compatible nodes (one by leaf)
	std::vector<uint64_t> chains;

	/** Create the structure for kmer colinear chaining.
	 * @param pairs kmer pair coordinates. This list is modified by the algorithm.
	 **/
	Colinear(std::vector<PairInt> & pairs);
	~Colinear();

	/** Compute the scores from the pairs used for construction. This function is out of the constructor as it can take a while to be computed.
	 * @param pairs The pairs here should be the same used as construction time.
	 * @param packed Use scores and ranks packed in 64 bits when there are less than 2^32 pairs
	 * (false: always 2 values, for the tests).
	 **/
	void compute_scores(std::vector<PairInt> & pairs, bool packed=true);

	/** Get the longest list of compatible kmer pairs.
	 * The compute_scores function needs to be called before extracting the longest chain.
	 * @return List of compatible pairs
	 **/
	std::vector<PairInt> longest_chain();

	// --- Construction functions ---

	/** Create a tree structure inside of the vector named tree using sorted_nodes as leaf keys.
	 * @param tree Modified vector that will contains the tree structure
	 * @param sorted_nodes Sorted nodes to construct the tree.
	 **/
	void init_tree(std::vector<Node> & tree, std::vector<PairInt> & sorted_nodes);

	/** Update the score of a node from the tree
	 * @param pair Pair to update
	 **/
	void update_score(PairInt & pair);
	/** Fills tree_idxs (pair -> leaf), only needed by update_score() */
	void index_leaves();

	/**
	 **/
	void propagate_up(uint64_t idx, uint64_t score);

	// --- Tree functions ---
	uint64_t up_tree(uint64_t idx);
	uint64_t left_subtree(uint64_t idx);
	uint64_t right_subtree(uint64_t idx);
	bool is_left(uint64_t idx);
	bool is_right(uint64_t idx);
};


#endif

