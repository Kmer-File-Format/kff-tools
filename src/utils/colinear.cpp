#include <cmath>
#include <algorithm>

#include "colinear.hpp"



#include <iostream>


using namespace std;


Colinear::Colinear(vector<PairInt> & pairs)
{
	// Insert the pairs as leves in the compatible tree with the right order
	sort(pairs.begin(), pairs.end(), 
    	[](const PairInt & a, const PairInt & b) -> bool
		{ 
			if (a.second != b.second)
				return a.second < b.second; 
			else
				return a.first > b.first;
		}
	);
	this->init_tree(this->nc_tree, pairs);

	// Create a list of previous compatible Node.
	this->chains = vector<uint64_t>(pairs.size(), this->nc_tree.size());
};


Colinear::~Colinear()
{

};


void Colinear::compute_scores(vector<PairInt> & pairs)
{
	// Update the pair scores one by one
	sort(pairs.begin(), pairs.end());
	for (PairInt & p : pairs)
	{
		// Update score
		this->update_score(p);
	}
}


void Colinear::update_score(PairInt & pair)
{
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
