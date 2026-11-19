#include <cmath>

#include "blocks.hpp"


static uint64_t bytes_from_bit_array(uint64_t bits_per_elem, uint64_t nb_elem) {
	if (bits_per_elem == 0 or nb_elem == 0)
		return 0;
	else
		return ((bits_per_elem * nb_elem - 1) / 8) + 1;
}


MinimizerBlockEncoder::MinimizerBlockEncoder(uint64_t k, uint64_t m, uint64_t max, uint64_t data_size)
	: k(k), m(m), max(max), data_size(data_size)
{
	// Same computations as the Section_Minimizer constructor of the API
	uint64_t nb_bits = static_cast<uint64_t>(ceil(log2(max)));
	this->nb_kmers_bytes = bytes_from_bit_array(nb_bits, 1);
	uint64_t mini_pos_bits = static_cast<uint8_t>(ceil(log2(k+max-1)));
	this->mini_pos_bytes = bytes_from_bit_array(mini_pos_bits, 1);
	this->minimizer_bytes = bytes_from_bit_array(2, m);
}


void MinimizerBlockEncoder::append(std::vector<uint8_t> & out, const uint8_t * seq, uint64_t seq_size,
                                   uint64_t mini_pos, const uint8_t * data) const {
	const uint64_t nb_kmers = seq_size - k + 1;

	// 1 - Number of kmers and minimizer position, big endian
	for (uint64_t b=this->nb_kmers_bytes ; b>0 ; b--)
		out.push_back(static_cast<uint8_t>(nb_kmers >> (8 * (b - 1))));
	for (uint64_t b=this->mini_pos_bytes ; b>0 ; b--)
		out.push_back(static_cast<uint8_t>(mini_pos >> (8 * (b - 1))));

	// 2 - Sequence without the minimizer, right aligned, padding bits set to 0
	const uint64_t no_mini_size = seq_size - m;
	const uint64_t suffix_size = no_mini_size - mini_pos;
	const size_t start = out.size();
	out.resize(start + bytes_from_bit_array(2, no_mini_size), 0);
	uint8_t * no_mini = out.data() + start;

	if (seq_size <= 32) {
		// 64 bits arithmetic
		uint64_t full = 0;
		for (uint64_t i=0 ; i<(seq_size + 3) / 4 ; i++)
			full = (full << 8) | seq[i];
		const uint64_t suffix = suffix_size == 0 ? 0 : full & ((uint64_t(1) << (2 * suffix_size)) - 1);
		const uint64_t prefix = mini_pos == 0 ? 0 : (full >> (2 * (suffix_size + m))) & ((uint64_t(1) << (2 * mini_pos)) - 1);
		uint64_t value = (suffix_size >= 32 ? 0 : prefix << (2 * suffix_size)) | suffix;
		for (uint64_t i=(no_mini_size + 3) / 4 ; i>0 ; i--) {
			no_mini[i - 1] = static_cast<uint8_t>(value);
			value >>= 8;
		}
	} else {
		const uint64_t seq_offset = (4 - seq_size % 4) % 4;
		uint64_t dst = (4 - no_mini_size % 4) % 4;
		auto copy = [&](uint64_t from, uint64_t nb) {
			for (uint64_t i=0 ; i<nb ; i++, dst++) {
				const uint64_t src = seq_offset + from + i;
				const uint8_t nucl = (seq[src / 4] >> (2 * (3 - src % 4))) & 0b11;
				no_mini[dst / 4] |= nucl << (2 * (3 - dst % 4));
			}
		};
		copy(0, mini_pos);
		copy(mini_pos + m, suffix_size);
	}

	// 3 - Data of the kmers
	out.insert(out.end(), data, data + bytes_from_bit_array(this->data_size * 8, nb_kmers));
}
