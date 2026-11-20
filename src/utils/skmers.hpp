#include <cstdint>
#include <vector>
#include <iterator>


#ifndef SKMERS_H
#define SKMERS_H


typedef struct interleved_s {
	uint8_t * nucl;
	std::size_t pref_size;
	std::size_t suf_size;
} interleved_t;


/** Compute the size of the interleaved part of a kmer.
 * @param k kmer size
 * @param m minimizer size
 * @mini_pos position of the minimizer inside of the kmer (also prefix size)
 * @return Number of nucleotides that can be included in interleaved sequences.
 **/
uint64_t interleaved_size(const uint64_t k, const uint64_t m, const uint64_t mini_pos);


/** Translate a superkmer into its interleaved representation.
 * 
 * @param skmer A Byte array containing the superkmer without its minimizer
 * @param interleaved_skmer A Byte array already allocated to put the interleaved in.
 * @param size The number of nucleotides inside of skmer (so the size without the minimizer)
 * @param mini_position The position of the minimizer inside of the skmer
 * 
 * @return An interleaved struct containing the skmer and its prefixe/suffix sizes
 */
interleved_t interleaved(const uint8_t * skmer, uint8_t * interleaved_skmer, std::size_t size, std::size_t mini_position);


/** Nucleotides of a kmer in the interleaved order around its minimizer, as a 128-bit integer.
 * The interleaved order is the order of interleaved(): for each distance t to the minimizer, the
 * suffix nucleotide (index pref + t), then the prefix one (index pref - 1 - t), when they exist.
 *
 * @param kmer Byte array of the nucleotides (2 bits each, from the most significant bits of the
 * first byte).
 * @param first_pos Position of the first nucleotide of the kmer in the byte array.
 * @param nb_nucl Number of nucleotides of the kmer (first_pos + nb_nucl <= 64).
 * @param pref Number of nucleotides before the minimizer (pref <= nb_nucl).
 * @param high, low The value, aligned on the right: the last nucleotide is in the 2 lowest bits
 * of low.
 */
void interleaved_value(const uint8_t * kmer, uint64_t first_pos, uint64_t nb_nucl, uint64_t pref,
                       uint64_t & high, uint64_t & low);


/** Compare 2 interleaved values
 * 
 * @param i1 A skmer in the interleaved format
 * @param i2 A skmer in the interleaved format
 * 
 * @return True if i1 is strictly befor i2 in the interleaved order
 */
bool inf_interleaved(interleved_t i1, interleved_t i2);


/** Selecte the min interleaved along an iterator
 * 
 * @param begin Start of the iterator to explore
 * @param end End of the iterator to explore
 * 
 * @return The first element in the interleaved order
 */
template<typename Iterator>
interleved_t min_interleaved(Iterator begin, Iterator end) {
	interleved_t min_i = *begin;
	for (Iterator it = ++begin; it != end; it++) {
		if (inf_interleaved(*it, min_i))
			min_i = *it;
	}
	return min_i;
}




#endif