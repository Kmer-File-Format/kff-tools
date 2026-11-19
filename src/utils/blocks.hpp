#include <cstdint>
#include <vector>

#ifndef BLOCKS_H
#define BLOCKS_H

/**
 * Encodes the blocks of a minimizer section in memory, with exactly the bytes
 * that Section_Minimizer::write_compacted_sequence() writes in a file:
 *   number of kmers | minimizer position | sequence without the minimizer | data
 * The width of the 2 first fields is computed as in the API, from k, m and max.
 */
class MinimizerBlockEncoder {
public:
	uint64_t k;
	uint64_t m;
	uint64_t max;
	uint64_t data_size;
	// Field widths in bytes
	uint64_t nb_kmers_bytes;
	uint64_t mini_pos_bytes;
	uint64_t minimizer_bytes;

	MinimizerBlockEncoder(uint64_t k, uint64_t m, uint64_t max, uint64_t data_size);

	/** Appends to out the bytes of the block of the sequence seq (seq_size nucleotides,
	 * right aligned) where the minimizer starts at mini_pos.
	 * @param data data_size bytes per kmer of the sequence.
	 */
	void append(std::vector<uint8_t> & out, const uint8_t * seq, uint64_t seq_size,
	            uint64_t mini_pos, const uint8_t * data) const;
};

#endif
