#include <cmath>
#include <algorithm>

#include "colinear.hpp"



#include <iostream>


using namespace std;


/* Stable counting sort of elements by a key in [0, max_key]. */
template <typename T, typename Key>
static void counting_sort(vector<T> & elements, uint64_t max_key, Key key) {
	vector<uint64_t> starts(max_key + 2, 0);
	for (const T & e : elements)
		starts[key(e) + 1] += 1;
	for (uint64_t k=1 ; k<starts.size() ; k++)
		starts[k] += starts[k - 1];
	vector<T> sorted(elements.size());
	for (const T & e : elements)
		sorted[starts[key(e)]++] = e;
	elements.swap(sorted);
}

/* The pair values are kmer positions: counting sorts are used when they are not too sparse. */
static bool dense_values(const vector<PairInt> & pairs, uint64_t & max_first, uint64_t & max_second) {
	max_first = max_second = 0;
	for (const PairInt & p : pairs) {
		max_first = std::max(max_first, p.first);
		max_second = std::max(max_second, p.second);
	}
	const uint64_t limit = 8 * pairs.size() + (1 << 20);
	return max_first <= limit and max_second <= limit;
}


Colinear::Colinear(vector<PairInt> & pairs)
{
	// Insert the pairs as leves in the compatible tree with the right order:
	// second value, then first value in decreasing order
	uint64_t max_first, max_second;
	if (pairs.size() > 1 and dense_values(pairs, max_first, max_second)) {
		// Same order in linear time: by decreasing first value, then (stable) by second value
		counting_sort(pairs, max_first, [max_first](const PairInt & p) { return max_first - p.first; });
		counting_sort(pairs, max_second, [](const PairInt & p) { return p.second; });
	} else {
		sort(pairs.begin(), pairs.end(), 
			[](const PairInt & a, const PairInt & b) -> bool
			{ 
				if (a.second != b.second)
					return a.second < b.second; 
				else
					return a.first > b.first;
			}
		);
	}
	this->init_tree(this->nc_tree, pairs);

	// Create a list of previous compatible Node.
	this->chains = vector<uint64_t>(pairs.size(), this->nc_tree.size());
};


Colinear::~Colinear()
{

};


void Colinear::compute_scores(vector<PairInt> & pairs)
{
	// Same result as calling update_score() for each pair in the (first, second) order, in
	// O(log n) per pair.
	//
	// When update_score() processes a pair (a, b), the scored leaves on its left are the pairs
	// (f, s) with f < a and s < b (compatible) and the previous pairs of its group (a, s) with
	// s < b (incompatible). update_score() looks, among the leaves with the max score M, for the
	// leftmost compatible one; if all of them are incompatible, the pair gets the chain of the
	// rightmost one. The pairs of a group are therefore processed together: their scores are
	// added to the tree at the end of the group, so that the tree only contains compatible
	// leaves, and the group max is tracked apart.
	// Leaves of the pairs in the (first, second) order. The leaves are sorted by second value
	// (then decreasing first value): a stable sort by first value is enough.
	const uint64_t no_leaf = this->nc_tree.size();
	vector<uint64_t> leaves(this->chains.size());
	for (size_t i=0 ; i<leaves.size() ; i++)
		leaves[i] = 2 * i;
	uint64_t max_first = 0;
	for (uint64_t leaf : leaves)
		max_first = std::max(max_first, this->nc_tree[leaf].first.first);
	if (max_first <= 8 * leaves.size() + (1 << 20))
		counting_sort(leaves, max_first, [this](uint64_t leaf) { return this->nc_tree[leaf].first.first; });
	else
		stable_sort(leaves.begin(), leaves.end(), [this](uint64_t a, uint64_t b) {
			return this->nc_tree[a].first.first < this->nc_tree[b].first.first;
		});
	// The pairs are left sorted in this order, as before
	pairs.resize(leaves.size());
	for (size_t i=0 ; i<leaves.size() ; i++)
		pairs[i] = this->nc_tree[leaves[i]].first;

	// The scores are computed in a segment tree in breadth first order (score of the leaf of
	// rank r, i.e. leaf 2r of nc_tree, in seg[seg_size + r]): its upper levels stay in the
	// cache, unlike the in-order layout of nc_tree, which is filled at the end.
	size_t seg_size = 1;
	while (seg_size < leaves.size())
		seg_size <<= 1;
	vector<uint64_t> seg(2 * seg_size, 0);
	// Max score of the leaves on the left of the rank, and leftmost leaf with this score
	auto left_max_leaf = [&seg, seg_size](uint64_t rank, uint64_t & max_leaf) -> uint64_t {
		uint64_t max_score = 0;
		uint64_t max_node = 0;
		// The left siblings met when going up cover the ranks on the left, further and further
		for (uint64_t node=seg_size+rank ; node>1 ; node>>=1)
			if ((node & 1) and seg[node - 1] > 0 and seg[node - 1] >= max_score) {
				max_score = seg[node - 1];
				max_node = node - 1;
			}
		if (max_score == 0)
			return 0;
		while (max_node < seg_size)
			max_node = seg[2 * max_node] == max_score ? 2 * max_node : 2 * max_node + 1;
		max_leaf = 2 * (max_node - seg_size);
		return max_score;
	};

	vector<pair<uint64_t, uint64_t> > group;  // leaf, score
	size_t group_start = 0;
	while (group_start < leaves.size())
	{
		size_t group_end = group_start;
		while (group_end < leaves.size() and this->nc_tree[leaves[group_end]].first.first == this->nc_tree[leaves[group_start]].first.first)
			group_end += 1;

		group.clear();
		// Max score in the group and rightmost leaf with that score (the group leaves are in
		// increasing tree order)
		uint64_t group_max = 0;
		uint64_t group_max_leaf = no_leaf;
		for (size_t i=group_start ; i<group_end ; i++)
		{
			uint64_t leaf = leaves[i];
			uint64_t max_leaf = no_leaf;
			uint64_t left_score = left_max_leaf(leaf / 2, max_leaf);
			uint64_t score = 1;

			if (left_score > 0 and left_score >= group_max)
			{
				// Leftmost compatible leaf with the max score
				this->chains[leaf/2] = max_leaf;
				score += left_score;
			}
			else if (group_max > 0)
			{
				// All the leaves with the max score are in the group
				this->chains[leaf/2] = this->chains[group_max_leaf/2];
				if (this->chains[leaf/2] < no_leaf)
					score += seg[seg_size + this->chains[leaf/2] / 2];
			}

			group.emplace_back(leaf, score);
			if (score >= group_max)
			{
				group_max = score;
				group_max_leaf = leaf;
			}
		}

		for (const pair<uint64_t, uint64_t> & leaf_score : group)
			for (uint64_t node=seg_size+leaf_score.first/2 ; node>=1 and seg[node]<leaf_score.second ; node>>=1)
				seg[node] = leaf_score.second;
		group_start = group_end;
	}

	// Scores of nc_tree: the leaves, then each internal node is the max of its children (as
	// propagate_up() makes it)
	for (size_t rank=0 ; rank<leaves.size() ; rank++)
		this->nc_tree[2 * rank].second = seg[seg_size + rank];
	for (uint64_t level=1 ; (uint64_t(1) << level) - 1 < this->nc_tree.size() ; level++)
		for (uint64_t i=(uint64_t(1) << level) - 1 ; i<this->nc_tree.size() ; i+=(uint64_t(1) << (level + 1)))
			this->nc_tree[i].second = std::max(this->nc_tree[i - (uint64_t(1) << (level - 1))].second,
			                                   this->nc_tree[i + (uint64_t(1) << (level - 1))].second);
}




void Colinear::index_leaves()
{
	// The leaves after the pairs only fill the tree
	for (uint64_t leaf=0 ; leaf<2*this->chains.size() ; leaf+=2)
		this->tree_idxs[this->nc_tree[leaf].first] = leaf;
}


void Colinear::update_score(PairInt & pair)
{
	if (this->tree_idxs.empty())
		this->index_leaves();
	uint64_t pair_idx = this->tree_idxs[pair];
	uint64_t idx = pair_idx;
	uint64_t max_score = 1; // One node alone like this one at the moment
	vector<uint64_t> max_idxs;
	
	// cout << pair.first << "-" << pair.second << endl;

	// Go up to the root to find the max non coliding value
	uint64_t up = this->up_tree(idx);
	while (up < this->nc_tree.size())
	{
		// If left node, compare max
		if (this->is_right(idx))
		{
			uint64_t left = this->left_subtree(up);
			Node & down_left = this->nc_tree[left];

			// If better or similar score
			if (down_left.second >= max_score)
			{
				if (down_left.second > max_score)
					max_idxs.clear();
				max_idxs.push_back(left);
				max_score = down_left.second;
			}
		}

		// Go one step up
		idx = up;
		up = this->up_tree(idx);
	}

	// cout << "max_idxs " << max_idxs.size() << endl << endl;

	// Explore the maxs to get the previous compatible
	while (max_idxs.size() > 0)
	{
		uint64_t idx = max_idxs.back();
		max_idxs.pop_back();
		uint64_t level = __builtin_ctz(~idx);

		// Go from the max idx node to one of the leaves with that max value
		while (level != 0)
		{
			uint64_t left_idx = this->left_subtree(idx);
			Node & left_node = this->nc_tree[left_idx];

			uint64_t right_idx = this->right_subtree(idx);
			Node & right_node = this->nc_tree[right_idx];

			// cout << idx << " <" << left_idx << " " << right_idx << ">" << endl;

			// Go for both subtrees
			if (left_node.second == max_score and right_node.second == max_score)
			{
				idx = left_idx;
				max_idxs.push_back(right_idx);
			}
			// Go left subtree
			else if (left_node.second == max_score)
			{
				idx = left_idx;
			}
			// Go right subtree
			else
			{
				idx = right_idx;
			}

			level = __builtin_ctz(~idx);			
		}

		Node & leaf = this->nc_tree[idx];
		PairInt & leaf_pair = leaf.first;
		// If the node is compatible
		if (pair.first != leaf_pair.first)
		{
			this->chains[pair_idx/2] = idx;
			// No need to search for other links, they will be equivalent in score at most.
			break;
		}
		// If not, register the value and continue
		else
		{
			// To be compatible, get the compatible parent of the current non compatible pair.
			this->chains[pair_idx/2] = this->chains[idx/2];
			// No break: continue to look for compatibles with the same score
		}
	}

	// Compute the new score
	uint64_t score = 1;
	// If there is a compatible predecessor
	if (this->chains[pair_idx/2] < this->nc_tree.size())
	{
		Node & prev_node = this->nc_tree[this->chains[pair_idx/2]];
		score += prev_node.second;
	}
	this->propagate_up(pair_idx, score);
};


vector<PairInt> Colinear::longest_chain()
{
	uint64_t max_score = 0;
	uint64_t max_idx = 0;

	// Search for the max values
	for (uint64_t i(0) ; i<this->nc_tree.size() ; i+=2)
	{
		Node & node = this->nc_tree[i];
		if (node.second > max_score)
		{
			max_score = node.second;
			max_idx = i;
		}
	}

	// Get the full chain of pairs
	vector<PairInt> chain(max_score, PairInt(0, 0));
	for (uint64_t score=max_score, idx=max_idx ; score>0 ; score--)
	{
		chain[score-1] = this->nc_tree[idx].first;
		idx = this->chains[idx/2];
	}

	return chain;
};


// ---------------------- Tree functions ----------------------


/** Create a tree structure inside of the vector named tree using sorted_nodes as leaf keys.
 * @param tree Modified vector that will contains the tree structure
 * @param sorted_nodes Sorted nodes to construct the tree.
 **/
void Colinear::init_tree(vector<Node> & tree, vector<PairInt> & sorted_nodes)
{
	// Compute the array size needed to store the tree
	uint64_t bits_needed = static_cast<uint64_t>(ceil(log2(sorted_nodes.size())));
	uint64_t tree_size = (1u << (bits_needed+1)) - 1u;
	tree.resize(tree_size);

	// Insert the pairs as root in the tree
	uint64_t tree_idx = 0;
	for (PairInt & p : sorted_nodes)
	{
		tree[tree_idx] = Node(p, 0);
		tree_idx += 2;
	}

	// Init intermediate nodes
	// For each level in the tree
	for (uint64_t level(1) ; level<=bits_needed ; level++)
	{
		// For each node from the current level
		for (uint64_t i((1 << level) - 1) ; i<tree_size ; i += (1 << (level + 1)))
		{
			// Get the right child
			uint64_t right_child_idx = i + (1 << (level - 1));
			Node & right_child = tree[right_child_idx];
			// Init the node
			tree[i] = Node(right_child.first, 0);
		}
	}
}


void Colinear::propagate_up(uint64_t idx, uint64_t score)
{
	while (idx < this->nc_tree.size())
	{
		Node & node = this->nc_tree[idx];
		// Nothing more to update
		if (score <= node.second)
			break;

		// Update the score
		node.second = score;
		// Go up on level
		idx = this->up_tree(idx);
	}
}


uint64_t Colinear::up_tree(uint64_t idx)
{
	// Lvl is equivalent to number of trailing 1s
	uint64_t level = __builtin_ctz(~idx);

	// Right part of the new index (one level up)
	uint64_t right_slice = (1u << (level + 1)) - 1;

	// Left part keeps the same bytes (everything except level and l/r bits)
	uint64_t left_mask = ~((1u << (level + 2)) - 1);
	uint64_t left_slice = idx & left_mask;

	return left_slice | right_slice;
};

uint64_t Colinear::left_subtree(uint64_t idx)
{
	// Lvl is equivalent to number of trailing 1s
	uint64_t level = __builtin_ctz(~idx);

	// Right part of the new index (one level down)
	uint64_t right_slice = (1u << (level - 1)) - 1;

	// Left part keeps the same bytes (everything except level and l/r bits)
	uint64_t left_mask = ~((1u << (level + 1)) - 1);
	uint64_t left_slice = idx & left_mask;

	return left_slice | right_slice;
};

uint64_t Colinear::right_subtree(uint64_t idx)
{
	// Lvl is equivalent to number of trailing 1s
	uint64_t level = __builtin_ctz(~idx);

	// Right part of the new index (one level down)
	uint64_t right_slice = (1u << (level - 1)) - 1;

	// Middle slice that indicates the rightness of the node
	uint64_t middle_slice = 1 << level;

	// Left part keeps the same bytes (everything except level and l/r bits)
	uint64_t left_mask = ~((1u << (level + 1)) - 1);
	uint64_t left_slice = idx & left_mask;

	return left_slice | middle_slice | right_slice;
};

bool Colinear::is_right(uint64_t idx)
{
	// Lvl is equivalent to number of trailing 1s
	uint64_t level = __builtin_ctz(~idx);

	uint64_t rightness = (idx >> (level + 1)) & 0b1;

	return rightness == 1;
};

bool Colinear::is_left(uint64_t idx)
{
	// Lvl is equivalent to number of trailing 1s
	uint64_t level = __builtin_ctz(~idx);

	uint64_t rightness = (idx >> (level + 1)) & 0b1;

	return rightness == 0;
};
