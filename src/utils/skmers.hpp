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