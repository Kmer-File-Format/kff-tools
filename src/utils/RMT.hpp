#include <vector>


#ifndef ZRMT_H
#define ZRMT_H

using namespace std;


template <class K>
class RangeMaxTree {
public:
	uint64_t max_real_leaf;
	vector<pair<K, uint64_t> > tree;
	uint64_t next_2pow;

	RangeMaxTree(const vector<K> keys) {
		this->max_real_leaf = keys.size()-1;
		// Compute the size for a perfect binary tree
		this->next_2pow = ceil(log2(keys.size()));
		uint64_t max_leaf_idx = (1 << next_2pow) - 1;
		// Resize the datastructure to have one space between each key (for the tree internal nodes)
		this->tree.resize((1 << (next_2pow + 1)) - 1);
		// Init the leaves with the keys
		for (uint idx=0 ; idx<keys.size() ; idx++)
			this->tree[idx*2] = pair<K, uint64_t>(keys[idx], 0);
		K max_key = keys[keys.size()-1];
		for (uint idx=keys.size() ; idx<= max_leaf_idx ; idx++) {
			this->tree[idx*2] = pair<K, uint64_t>(max_key, 0);
			this->tree[idx*2].second = 0;
		}
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
		// uint64_t end = (this->tree.size() + 1) / 2;
		uint64_t end = this->max_real_leaf;

		while (begin <= end) {
			uint64_t middle = (begin + end) / 2;
			if (this->tree[middle*2].first == key) {
				return middle*2;
			}
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
			if (val > this->tree[current_position].second) {
				this->tree[current_position].second = val;
				key_position >>= 1;
			} else
				return;
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


	/** Range max query
	 */
	uint64_t range(K left, K right) {
		if (right < left)
			return 0;
		if (left == right)
			return this->tree[this->find(left)].second;

		// Key positions
		uint64_t current_left = this->find(left);
		uint64_t left_position = current_left / 2;
		uint64_t current_right = this->find(right);
		uint64_t right_position = current_right / 2;
		uint64_t diff = left_position ^ right_position;
		uint64_t leftmost_common_bit = floor(log2(diff));
		
		// Init
		uint64_t left_max = this->tree[current_left].second;
		uint64_t right_max = this->tree[current_right].second;

		// Find max value up the tree
		for (uint lvl=0 ; lvl<leftmost_common_bit ; lvl++) {
			uint64_t shift = 1 << lvl;

			// --- Left max computation ---
			// Left position up in the tree from right child
			if (left_position % 2 == 1)
				current_left -= shift;
			// Left position up in the tree from left child
			else {
				current_left += shift;
				// Change the max if left child is better
				if (this->tree[current_left + shift].second > this->tree[current_left - shift].second) {
					left_max = this->tree[current_left + shift].second;
				}
			}
			left_position >>= 1;

			// --- Right max computation ---
			// Right position up in the tree from right child
			if (right_position % 2 == 0)
				current_right += shift;
			// Left position up in the tree from left child
			else {
				current_right -= shift;
				// Change the max if left child is better
				if (this->tree[current_right - shift].second > this->tree[current_right + shift].second) {
					right_max = this->tree[current_right - shift].second;
				}
			}
			right_position >>= 1;
		}

		// Combine max
		return max(left_max, right_max);
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
			return K();
	}


	/** Search for the first key after the bound with the minimal value value
	 */
	K bounded_first_max_key(uint64_t value, K bound) {
		uint64_t current_position = this->find(bound);
		uint64_t position = current_position/2;

		bool found_value = this->tree[current_position].second == value;
		uint64_t lvl = 0;
		// Go up in the tree from the bound ; until we reach the value searched.
		while (lvl <= this->next_2pow and not found_value) {
			// One level up
			uint64_t shift = 1 << lvl;
			if (position % 2 == 0)
				current_position += shift;
			else
				current_position -= shift;
			position >>= 1;
			lvl += 1;

			// Verify right subtree
			if (this->tree[current_position + shift].second >= value) {
				current_position += shift;
				found_value = true;
				lvl -= 1;
			}
		}

		// Verify the lvl
		if (lvl > this->next_2pow)
			return K();

		// Go down leftmost leave equaling the value
		while (lvl != 0) {
			uint64_t shift = 1 << lvl;
			if (this->tree[current_position - shift].second >= value) {
				current_position -= shift;
			} else {
				current_position += shift;
			}
			lvl -= 1;
		}

		return this->tree[current_position].first;
	}


	void print() {
		for (const auto &p : this->tree)
            cout << "(" << p.first << "," << p.second << ") ";
        cout << endl;
	}
};

#endif