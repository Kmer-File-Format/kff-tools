#include <vector>
#include <unordered_map>
#include <cstdint>


#ifndef COLINEAR_H
#define COLINEAR_H

// Redefine pair operators
typedef  std::pair<uint64_t, uint64_t> PairInt;
typedef std::pair<PairInt, uint64_t> Node;

struct hash_PairInt {
    size_t operator()(const PairInt& p) const
    {
        if (p.first != p.second) {
            return p.first ^ p.second;             
        }
        return p.first;
    }
};


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
	 **/
	void compute_scores(std::vector<PairInt> & pairs);

	/** Get the longest list of compatible kmer pairs.
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

