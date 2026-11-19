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


namespace {

/** 128-bit integer */
struct U128 {
	uint64_t hi;
	uint64_t lo;
};

inline U128 shift_left(U128 x, uint64_t k) {
	if (k == 0) return x;
	if (k >= 128) return {0, 0};
	if (k >= 64) return {x.lo << (k - 64), 0};
	return {(x.hi << k) | (x.lo >> (64 - k)), x.lo << k};
}

inline U128 shift_right(U128 x, uint64_t k) {
	if (k == 0) return x;
	if (k >= 128) return {0, 0};
	if (k >= 64) return {0, x.hi >> (k - 64)};
	return {x.hi >> k, (x.lo >> k) | (x.hi << (64 - k))};
}

/** The k lowest bits of x */
inline U128 low_bits(U128 x, uint64_t k) {
	if (k >= 128) return x;
	if (k >= 64) return {k == 64 ? 0 : x.hi & (~uint64_t(0) >> (128 - k)), x.lo};
	return {0, k == 0 ? 0 : x.lo & (~uint64_t(0) >> (64 - k))};
}

/** Reverses the order of the 32 nucleotides (2-bit groups) of x */
inline uint64_t reverse_nucleotides(uint64_t x) {
	x = (x >> 32) | (x << 32);
	x = ((x >> 16) & 0x0000FFFF0000FFFFULL) | ((x & 0x0000FFFF0000FFFFULL) << 16);
	x = ((x >> 8) & 0x00FF00FF00FF00FFULL) | ((x & 0x00FF00FF00FF00FFULL) << 8);
	x = ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((x & 0x0F0F0F0F0F0F0F0FULL) << 4);
	x = ((x >> 2) & 0x3333333333333333ULL) | ((x & 0x3333333333333333ULL) << 2);
	return x;
}

/** The 16 nucleotides of v, each one followed by 2 null bits */
inline uint64_t spread_nucleotides(uint32_t v) {
	uint64_t x = v;
	x = (x | (x << 16)) & 0x0000FFFF0000FFFFULL;
	x = (x | (x << 8)) & 0x00FF00FF00FF00FFULL;
	x = (x | (x << 4)) & 0x0F0F0F0F0F0F0F0FULL;
	x = (x | (x << 2)) & 0x3333333333333333ULL;
	return x;
}

/** The 32 nucleotides of v, each one followed by 2 null bits */
inline U128 spread_nucleotides(uint64_t v) {
	return {spread_nucleotides(static_cast<uint32_t>(v >> 32)), spread_nucleotides(static_cast<uint32_t>(v))};
}

/** x shifted, 0 from 64 bits */
inline uint64_t shift_left64(uint64_t x, uint64_t k) { return k >= 64 ? 0 : x << k; }
inline uint64_t shift_right64(uint64_t x, uint64_t k) { return k >= 64 ? 0 : x >> k; }
/** The k lowest bits of x */
inline uint64_t low_bits64(uint64_t x, uint64_t k) { return k >= 64 ? x : x & ((uint64_t(1) << k) - 1); }


/** interleaved_value() when the nucleotides are in the first 32 positions: 64-bit operations */
void interleaved_value64(const uint8_t * kmer, uint64_t first_pos, uint64_t nb_nucl, uint64_t pref,
                                uint64_t & high, uint64_t & low) {
	const uint64_t nb_bytes = (first_pos + nb_nucl + 3) / 4;
	uint64_t value = 0;
	for (uint64_t b=0 ; b<nb_bytes ; b++)
		value = (value << 8) | kmer[b];
	value = low_bits64(value >> (2 * (4 * nb_bytes - first_pos - nb_nucl)), 2 * nb_nucl);

	const uint64_t suf = nb_nucl - pref;
	const uint64_t common = pref < suf ? pref : suf;
	const uint64_t suffix = low_bits64(value, 2 * suf);
	const uint64_t prefix = shift_right64(value, 2 * suf);
	const uint64_t rev_prefix = shift_right64(reverse_nucleotides(prefix), 2 * (32 - pref));

	// common <= 16: the heads have at most 32 bits
	const uint64_t suffix_head = shift_right64(suffix, 2 * (suf - common));
	const uint64_t prefix_head = shift_right64(rev_prefix, 2 * (pref - common));
	uint64_t result = (spread_nucleotides(static_cast<uint32_t>(suffix_head)) << 2)
	                  | spread_nucleotides(static_cast<uint32_t>(prefix_head));

	const uint64_t tail_size = nb_nucl - 2 * common;
	const uint64_t tail = low_bits64(suf > pref ? suffix : rev_prefix, 2 * tail_size);
	high = 0;
	low = shift_left64(result, 2 * tail_size) | tail;
}


} // namespace


void interleaved_value(const uint8_t * kmer, uint64_t first_pos, uint64_t nb_nucl, uint64_t pref,
                       uint64_t & high, uint64_t & low) {
	if (first_pos + nb_nucl <= 32) {
		interleaved_value64(kmer, first_pos, nb_nucl, pref, high, low);
		return;
	}

	// Nucleotides of the kmer, the first one in the most significant bits
	const uint64_t nb_bytes = (first_pos + nb_nucl + 3) / 4;
	U128 value = {0, 0};
	for (uint64_t b=0 ; b<nb_bytes ; b++) {
		value.hi = (value.hi << 8) | (value.lo >> 56);
		value.lo = (value.lo << 8) | kmer[b];
	}
	value = low_bits(shift_right(value, 2 * (4 * nb_bytes - first_pos - nb_nucl)), 2 * nb_nucl);

	const uint64_t suf = nb_nucl - pref;
	const uint64_t common = pref < suf ? pref : suf;
	// Suffix in order and prefix in reverse order: the nucleotides by distance to the minimizer,
	// the closest one in the most significant bits
	const U128 suffix = low_bits(value, 2 * suf);
	const U128 prefix = shift_right(value, 2 * suf);
	const U128 rev_prefix = shift_right({reverse_nucleotides(prefix.lo), reverse_nucleotides(prefix.hi)}, 2 * (64 - pref));

	// The common distances interleaved (S0 P0 S1 P1 ...), then the rest of the longest part
	const uint64_t suffix_head = shift_right(suffix, 2 * (suf - common)).lo;
	const uint64_t prefix_head = shift_right(rev_prefix, 2 * (pref - common)).lo;
	U128 result = shift_left(spread_nucleotides(suffix_head), 2);
	const U128 spread_prefix = spread_nucleotides(prefix_head);
	result.hi |= spread_prefix.hi;
	result.lo |= spread_prefix.lo;

	const uint64_t tail_size = nb_nucl - 2 * common;
	const U128 tail = low_bits(suf > pref ? suffix : rev_prefix, 2 * tail_size);
	result = shift_left(result, 2 * tail_size);
	high = result.hi | tail.hi;
	low = result.lo | tail.lo;
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


uint64_t interleaved_size(const uint64_t k, const uint64_t m, const uint64_t mini_pos)
{
	uint64_t pref_size = mini_pos;
	uint64_t suff_size = k - m - pref_size;

	if (pref_size < suff_size)
		return pref_size * 2 + 1;
	else
		return suff_size * 2;
}


interleved_t interleaved(const uint8_t * skmer, uint8_t * interleaved_skmer, size_t size, size_t mini_position) {
	// Init the interleaved struct
	interleved_t interleaved;
	interleaved.nucl = interleaved_skmer;
	interleaved.pref_size = mini_position;
	interleaved.suf_size = size - interleaved.pref_size;
	
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