#include <stdint.h>
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>
#include <utility>

#include "encoding.hpp"
#include "kff-cpp-api/kff_io.hpp"

/* TO READ BEFORE ANY USAGE !!

This file provide function for dna sequences represented as a uint8_t array.
The arrays are big endian represented.
The begining of the prefix is present at the byte 0.
Byte 0 also have useless bits at the begining if the sequence size is not 4 divisable.

*/

#ifndef SEQ_H
#define SEQ_H

/** Generic object to read sequence files.
 * Each call to next_sequence should return a sequence
 **/
class SequenceStream {
public:
  /** Read the next sequence in the file.
    * If this sequence is larger than max_seq_size, it will fill seq with the max_seq_size first
    * nucleotides and return the full sequence size. remaining nucleotides are discarded.
    * @param seq Pointer that will be filled with the sequence array. The array is erased during the
    * next call to the method.
    * @return size of the sequence that have been read.
    **/
  virtual uint next_sequence(uint8_t * & seq, uint8_t * & data) = 0;
};

/** Read txt sequence file (1 sequence per line)
 */
class TxtSeqStream : public SequenceStream {
private:
  std::fstream fs;
  Binarizer bz;

  uint buffer_size;
  uint8_t * buffer;

public:
  TxtSeqStream(const std::string filename, const uint8_t encoding[4]) 
      : fs(filename, std::fstream::in)
      , bz(encoding)
      , buffer_size(1024)
      , buffer(new uint8_t[1024])
  {};
  ~TxtSeqStream() {
    this->fs.close();
    delete[] buffer;
  }
  uint next_sequence(uint8_t * & seq, uint8_t * & data);
};


/** Read txt sequence file (1 sequence per line)
 */
class KffSeqStream : public SequenceStream {
private:
public:
  Kff_reader reader;
  KffSeqStream(const std::string filename)
      : reader(filename)
  {};
  uint next_sequence(uint8_t * & seq, uint8_t * & data);
};


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

/** Translate a binarized uint sequence into an binarized sequence array.
  * @param seq Sequence stored in a uint (ie max 32 nucleotides)
  * @param bin_seq A Byte array to store the sequence. Must be allocated.
  * @param size Number of nucleotides in the sequence.
  */
void uint_to_seq(uint seq, uint8_t * bin_seq, uint size);


// ----- Minimizer search related functions -----

/** Compute all the candidates hash values of the sequence and return them into a vactor.
 * @param seq binarized sequence
 * @param size seq size in nucleotides
 * @param k kmer size
 * @param m minimizer size max = 31
 * @return vector containing all the hashed m-size windows
 **/
std::vector<uint64_t> compute_mini_candidates(const uint8_t * seq, const uint size, const uint k, const uint m);
/** Compute all the minimizers of a sequence.
 * @return All pair minimizer position/
 **/
std::vector<std::pair<int, uint64_t> > compute_minizers(const uint8_t * seq, const uint size, const uint k, const uint m);
/** Compute all the superkmers
 * @return All the begin/end pair positions
 **/
std::vector<std::pair<uint, uint> > compute_skmers(const uint seq_size, const uint k, const uint m, std::vector<std::pair<int, uint64_t> > & minimizers);
/** Compute all the superkmers
 **/
std::vector<std::pair<uint, uint> > compute_skmers(const uint8_t * seq, const uint size, const uint k, const uint m);

/** Search for the minimizer inside of a sequence (forward only)
  * @param seq binarized sequence.
  * @param size size in nucleotide of the sequence
  * @param m Minimizer size (max 32)
  * @param minimizer Modified during the execution to store the minimizer.
  * @param minimizer_position Modified during the execution to store the position of the minimizer.
  */
void search_mini(uint8_t * seq, const uint size, const uint m, uint & minimizer, uint & minimizer_position);


// ----- Usefull binary functions -----
void leftshift8(uint8_t * bitarray, size_t length, size_t bitshift);
void rightshift8(uint8_t * bitarray, size_t length, size_t bitshift);
uint8_t fusion8(uint8_t left_bits, uint8_t right_bits, size_t merge_index);

#endif
