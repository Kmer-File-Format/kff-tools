#include <cstdint>
#include <vector>


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


/** Selecte the min interleaved in a list
 * 
 * @param A list of superkmers interleaved
 * 
 * @return The first one in the interleaved order
 */
interleved_t min_interleaved(std::vector<interleved_t> skmers);


#endif