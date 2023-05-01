#include <cmath>
#include <algorithm>

#include "colinear.hpp"



#include <iostream>


using namespace std;


Colinear::Colinear(vector<PairInt> & pairs)
{
	// Insert the pairs as root in the tree with the right order
	sort(pairs.begin(), pairs.end(), 
    	[](const PairInt & a, const PairInt & b) -> bool
		{ 
			if (a.second != b.second)
				return a.second < b.second; 
			else
				return a.first < b.first;
		}
	);
	this->init_tree(this->nc_tree, pairs);

	// Update the pair scores one by one
	sort(pairs.begin(), pairs.end());
	for (PairInt & p : pairs)
	{
		// Update score
		this->update_score(p);
	}
};


Colinear::~Colinear()
{

};


void Colinear::update_score(PairInt & pair)
{
	uint64_t idx = this->tree_idxs[pair];
	uint64_t max_score = 0;
	uint64_t max_idx = idx;
	
	// Go up to the root to find the max non coliding value
	uint64_t up = this->up_tree(tree_idx);
	while (up < this->nc_tree.size())
	{
		// If left node, compare max
		if (this->is_left(idx))
		{
			Node & up_node = this->nc_tree[up];
		}

		// Go one step up
		idx = up;
		up = this->up_tree(tree_idx);
	}
};


vector<PairInt> Colinear::longest_chain()
{
	return vector<PairInt>();
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
		this->tree_idxs[p] = tree_idx;
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


uint64_t Colinear::up_tree(uint64_t idx)
{
	// Lvl is equivalent to number of trailing 1s
	uint64_t level = __builtin_ctz(~current_idx);

	// Right part of the new index (one level up)
	uint64_t right_slice = (1u << (level + 1)) - 1;

	// Left part keeps the same bytes (everything except level and l/r bits)
	uint64_t left_mask = ~((1u << (level + 2)) - 1);
	uint64_t left_slice = current_idx & left_mask;

	return left_slice | right_slice;
};

uint64_t Colinear::left_subtree(uint64_t idx)
{
	// Lvl is equivalent to number of trailing 1s
	uint64_t level = __builtin_ctz(~current_idx);

	// Right part of the new index (one level down)
	uint64_t right_slice = (1u << (level - 1)) - 1;

	// Left part keeps the same bytes (everything except level and l/r bits)
	uint64_t left_mask = ~((1u << (level + 1)) - 1);
	uint64_t left_slice = current_idx & left_mask;

	return left_slice | right_slice;
};

uint64_t Colinear::right_subtree(uint64_t idx)
{
	// Lvl is equivalent to number of trailing 1s
	uint64_t level = __builtin_ctz(~current_idx);

	// Right part of the new index (one level down)
	uint64_t right_slice = (1u << (level - 1)) - 1;

	// Middle slice that indicates the rightness of the node
	uint64_t middle_slice = 1 << level;

	// Left part keeps the same bytes (everything except level and l/r bits)
	uint64_t left_mask = ~((1u << (level + 1)) - 1);
	uint64_t left_slice = current_idx & left_mask;

	return left_slice | middle_slice | right_slice;
};

bool Colinear::is_right(uint64_t idx)
{
	// Lvl is equivalent to number of trailing 1s
	uint64_t level = __builtin_ctz(~current_idx);

	uint64_t rightness = (idx >> (level + 1)) & 0b1;

	return rightness == 1;
};

bool Colinear::is_left(uint64_t idx)
{
	// Lvl is equivalent to number of trailing 1s
	uint64_t level = __builtin_ctz(~current_idx);

	uint64_t rightness = (idx >> (level + 1)) & 0b1;

	return rightness == 0;
};
