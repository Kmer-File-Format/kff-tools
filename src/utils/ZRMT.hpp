#include <vector>


#ifndef ZRMT_H
#define ZRMT_H

using namespace std;


template <class K>
class ZeroRangeMaxTree {
public:
	uint64_t max_real_leaf;
	vector<pair<K, uint64_t> > tree;
	uint64_t next_2pow;

	ZeroRangeMaxTree(const vector<K> keys) {
		this->max_real_leaf = keys.size();
		// Compute the size for a perfect binary tree
		this->next_2pow = ceil(log2(keys.size()));
		uint64_t max_leaf_idx = (1 << next_2pow) - 1;
		// Resize the datastructure to have one space between each key (for the tree internal nodes)
		this->tree.resize((1 << (next_2pow + 1)) - 1);
		// Init the leaves with the keys
		for (uint idx=0 ; idx<keys.size() ; idx++)
			this->tree[idx*2] = pair<K, uint64_t>(keys[idx], 0);
		for (uint idx=keys.size() ; idx<= max_leaf_idx ; idx++)
			this->tree[idx*2] = pair<K, uint64_t>(NULL, 0);
		// Init internal nodes
		for (uint64_t lvl=1 ; lvl<=next_2pow ; lvl++) {
			for (uint64_t idx=(1<<lvl)-1 ; idx<this->tree.size() ; idx += 1 << (lvl + 1)) {
				this->tree[idx] = pair<K, uint64_t>(this->tree[idx+(1<<(lvl-1))].first, 0);
			}
		}
	}

	/** Find the position of key in the tree
	 */
	uint64_t find(K key) {
		uint64_t begin = 0;
		uint64_t end = (this->tree.size() + 1) / 2;

		while (begin <= end) {
			uint64_t middle = (begin + end) / 2;
			if (this->tree[middle*2].first == key)
				return middle*2;
			else if (this->tree[middle*2].first < key) {
				begin = middle + 1;
			} else {
				end = middle - 1;
			}
		}

		return this->tree.size();
	}

	/** Update the value linked to key
	 */
	void update(K key, uint64_t val) {
		uint64_t current_position = this->find(key);
		// Position in original vector
		uint64_t key_position = current_position / 2;
		// Update the leaf
		this->tree[current_position].second = val;
		// Update internal values
		for (uint lvl=0 ; lvl<this->next_2pow ; lvl++) {
			uint64_t shift = 1 << lvl;
			if (key_position % 2 == 0) {
				current_position += shift;
			} else {
				current_position -= shift;
			}
			this->tree[current_position].second = val;
			key_position >>= 1;
		}
	}

	/** Range max query with left boundary at 0
	 */
	uint64_t zero_range(K key) {
		// Key positions
		uint64_t current_position = this->find(key);
		uint64_t key_position = current_position / 2;
		
		// Init
		uint64_t max = this->tree[current_position].second;

		// Find max value up the tree
		for (uint lvl=1 ; lvl<=next_2pow ; lvl++) {
			uint64_t shift = 1 << (lvl-1);

			// Up in the tree from left child
			if (key_position % 2 == 0) {
				// No change from left
				current_position += shift;
			}
			// Up in the tree from right child
			else {
				current_position -= shift;
				// Change the max if left child is better
				if (this->tree[current_position - shift].second > max) {
					max = this->tree[current_position - shift].second;
				}
			}
			key_position >>= 1;
		}

		return max;
	}


	/** Search for the first key that corresponds to some max value
	 */
	K first_max_key(uint64_t value) {
		uint64_t current_position = this->tree.size() / 2;

		for (uint lvl=next_2pow ; lvl>0 ; lvl--) {
			uint shift = 1 << (lvl - 1);

			// If possible to go down left, GO
			if (this->tree[current_position - shift].second >= value)
				current_position -= shift;
			else
				current_position += shift;
		}

		// Is the max here ?
		if (this->tree[current_position].second == value)
			return this->tree[current_position].first;
		// Not found
		else
			return NULL;
	}


	void print() {
		for (const auto &p : this->tree)
            cout << "(" << p.first << "," << p.second << ") ";
        cout << endl;
	}
};

#endif