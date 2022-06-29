#include <iostream>
#include "skmers.hpp"

using namespace std;


uint8_t get_nucl(const uint8_t * seq, size_t position) {
	size_t byte_idx = position / 4;
	uint8_t nucl = seq[byte_idx];

	size_t nucl_shift = 3 - (position % 4);
	nucl >>= nucl_shift * 2;

	return nucl & 0b11;
}


void set_nucl(uint8_t * seq, size_t position, uint8_t nucl) {
	size_t byte_idx = position / 4;
	size_t nucl_shift = 3 - (position % 4);

	// Remove previous bits
	uint8_t mask = 0b11 << (nucl_shift * 2);
	seq[byte_idx] &= ~mask;

	// Add new nucl
	nucl <<= nucl_shift * 2;
	seq[byte_idx] += nucl;
}


bool inf_interleaved(interleved_t i1, interleved_t i2) {
	size_t max_inter_i1_size = i1.pref_size < i1.suf_size ? 2 * i1.pref_size + 1 : 2 * i1.suf_size;
	size_t max_inter_i2_size = i2.pref_size < i2.suf_size ? 2 * i2.pref_size + 1 : 2 * i2.suf_size;
	size_t max_comparable = min(max_inter_i1_size, max_inter_i2_size);

	// Interleaved comparison
	for (size_t i=0 ; i<max_comparable ; i++) {
		uint8_t mask = 0;
		switch(i % 4) {
			case 0: mask = 0b11000000; break;
			case 1: mask = 0b00110000; break;
			case 2: mask = 0b00001100; break;
			case 3: mask = 0b00000011; break;
		}

		uint8_t n1 = i1.nucl[i/4] & mask;
		uint8_t n2 = i2.nucl[i/4] & mask;
		if (n1 != n2) {
			return n1 < n2;
		}
	}

	// If interleaved is equal and prefix/suffix are symetric in between the 2 interleaved, continue comparison
	if (i1.pref_size == i2.pref_size) {
		for (size_t i=max_comparable ; i<i1.pref_size+i1.suf_size ; i++) {
			uint8_t mask = 0;
			switch(i % 4) {
				case 0: mask = 0b11000000; break;
				case 1: mask = 0b00110000; break;
				case 2: mask = 0b00001100; break;
				case 3: mask = 0b00000011; break;
			}

			uint8_t n1 = i1.nucl[i/4] & mask;
			uint8_t n2 = i2.nucl[i/4] & mask;
			if (n1 != n2) {
				return n1 < n2;
			}
		}	
	}

	return false;
}


interleved_t interleaved(const uint8_t * skmer, uint8_t * interleaved_skmer, size_t size, size_t mini_position) {
	// Init the interleaved struct
	interleved_t interleaved;
	interleaved.nucl = interleaved_skmer;
	interleaved.pref_size = mini_position;
	interleaved.suf_size = size - mini_position;

	// Alias usefull variables
	size_t pref = interleaved.pref_size;
	size_t suf = interleaved.suf_size;

	size_t current_pref_idx = pref - 1;
	size_t current_suf_idx = 0;

	// Transform the skmer into interleaved
	for (size_t i=0 ; i<size ; i++) {
		// Position to extract the nucl
		size_t position;
		// No prefix
		if (pref == 0)
			position = mini_position + current_suf_idx++;
		// No suffix
		else if (suf == 0)
			position = current_pref_idx--;
		// look in the suffix
		else if (i % 2 == 0) {
			position = mini_position + current_suf_idx++;
			suf -= 1;
		}
		// Look in the prefix
		else {
			position = current_pref_idx--;
			pref -= 1;
		}

		// Extract nucl from skmer
		uint8_t nucl = get_nucl(skmer, position);
		// Insert nucl
		set_nucl(interleaved_skmer, i, nucl);
	}

	return interleaved;
}