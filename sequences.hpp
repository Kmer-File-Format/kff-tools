#include <stdint.h>
#include <cstdlib>

/* TO READ BEFORE ANY USAGE !!

This file provide function for dna sequences represented as a uint8_t array.
The arrays are big endian represented.
The begining of the prefix is present at the byte 0.
Byte 0 also have useless bits at the begining if the sequence size is not 4 divisable.

*/

#ifndef SEQ_H
#define SEQ_H

/** Extract a subsequence of sequence.
  * @param sequence Original sequence
  * @param seq_size Size in nucleotides of the sequence
  * @param extracted A memory space already allocated by the user to copy the subsequence.
  * This memory must be 1 byte larger than the requiered space for the subsequence.
  * @param begin_nucl first nucleotide to extract (between 0 and seq_size-1)
  * @param end_nucl last nucleotide to extract (between 0 and seq_size-1)
  */
void subsequence(const uint8_t * sequence, const uint seq_size, uint8_t * extracted, const uint begin_nucl, const uint end_nucl);


/** Translate a sequence to an 64 bits integer. If the sequence length is more than 32, then only
	* the last 32 nucleotides are used for the conversion.
  *
  * @param seq A nucleotide sequence to translate
  * @param seq_size Size in nucleotides of the sequence
  */
uint64_t seq_to_uint(const uint8_t * seq, uint seq_size);


// ----- Usefull binary functions -----
void leftshift8(uint8_t * bitarray, size_t length, size_t bitshift);
void rightshift8(uint8_t * bitarray, size_t length, size_t bitshift);
uint8_t fusion8(uint8_t left_bits, uint8_t right_bits, size_t merge_index);

#endif
