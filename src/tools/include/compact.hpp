#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>

#include "CLI/CLI.hpp"
#include "kfftools.hpp"


#ifndef COMPACT_H
#define COMPACT_H

class Compact: public KffTool {
private:
	std::string input_filename;
	std::string output_filename;


public:
	uint k;
	uint m;
	uint data_size;
	uint mini_pos_size;
	uint bytes_compacted;
	uint offset_idx;
	bool sorted;

	uint8_t * kmer_buffer;
	uint64_t buffer_size;
	uint64_t next_free;

	Compact();
	~Compact();


	/** Write a minimizer section containing all the superkmers given as paths. The process preserve
	 * the order of the superkmers.
	 * 
	 * @param paths List of all the kmer paths that represent virtual superkmers.
	 * @param sm Section Minimizer to fill.
	 * @param data_size Size in bytes of each data.
	 **/
	void write_paths(const std::vector<std::vector<uint8_t *> > & paths, Section_Minimizer & sm, const uint data_size);

	/** Take a list of kmer pair to assemble and return a list of paths from left to right of
	 * virtual superkmers.
	 * 
	 * @param to_compact Pairs of kmers to compact into paths
	 * 
	 * @return A list of virtual superkmer paths (path = kmer list from left to right).
	 **/
	std::vector<std::vector<uint8_t *> > pairs_to_paths(const std::vector<std::pair<uint8_t *, uint8_t *> > & to_compact);

	/** Add a kmer, its data and its minimizer position in the compaction buffer.
	 * WARNING: All the requested values must be set in the compact object (ie k, m, data_size, mini_pos_size, bytes_compacted)
	 * @param seq kmer to add. Must be compacted without minimizer.
	 * @param data data related to the kmer to add (can be a nullpointer on data_size=0)
	 * @param mini_pos Minimizer position inside of the kmer.
	 * 
	 * @return The position of the kmer (in bytes) inside of the buffer)
	 **/
	long add_kmer_to_buffer(const uint8_t * seq, const uint8_t * data, const uint64_t mini_pos);

	/** Extract a kmer minimizer position from the kmer buffer.
	 * 
	 * @param pos The kmer position in the buffer
	 * 
	 * @return The minimizer position inside of the kmer
	 **/
	uint mini_pos_from_buffer(const uint8_t * kmer) const;
	
	std::vector<std::vector<uint8_t *> > prepare_kmer_matrix(Section_Minimizer & sm);

	/** Return the result of the comparison between kmers in the buffer.
	 * WARNING: The comparator assumes that the minimizers are at the same place in the words
	 * 
	 * @param kmer1 First kmer in the buffer
	 * @param kmer2 Second kmer in the buffer
	 * 
	 * @return A negative number is the first kmer is smaller, 0 if both are equal, a positive number
	 * otherwise.
	 **/
	int interleaved_compare_kmers(const uint8_t * kmer1, const uint8_t * kmer2) const;

	/** Sort each column of the matrix using a kmer order. The input columns of the matrix are
	 * modified during this process.
	 * 
	 * @param kmer_matrix A matrix where all the kmers of the same column share the same minimizer 
	 * position.
	 **/
	void sort_matrix(std::vector<std::vector<uint8_t *> > & kmer_matrix);
	
	/** Take a succesive pair of columns of the sorted matrix and output the kmer
	 * pairs that are overlaping.
	 * 
	 * @param column1 column of the matrix for left kmers
	 * @param column2 column of the matrix for right kmers
	 * @return A vector of all overlaping pairs. A pair corresponds to both kmer positions in
	 * their original vector.
	 **/
	std::vector<std::pair<uint64_t, uint64_t> > pair_kmers(const std::vector<uint8_t *> & column1, const std::vector<uint8_t *> & column2) const;

	/** Performs a Longest increasing subsequence on a sorted vector of potential kmer overlaps.
	 * The goal here is to select the maximum number of links (to maximize the compaction) preserving
	 * the kmer order for each column of the matrix (order on the kmers that share the same minimizer
	 * position)
	 * 
	 * @param candidates Sorted list of overlaping candidate. This list should be sorted by the first
	 * kmer order, then the second kmer for equalities.
	 * 
	 * @return The list of selected links. All other links are removed to keep the order.
	 **/
	std::vector<std::pair<uint64_t, uint64_t> > colinear_chaining(const std::vector<std::pair<uint64_t, uint64_t> > & candidates) const;

	/** From the kmer matrix and the list of all the colinear chained pairs of kmers, generate the ordered list of superkmers.
	 * 
	 * @param matrix The full kmer matrix 2D sorted
	 * @param colinear_chainings sorted pairs of kmers to assemble (number of matrix columns - 1 lists)
	 * @return The list of sorted linked kmers. They can be compacted in a skmers style.
	 **/
	std::vector<std::vector<uint8_t *> > polish_sort(const std::vector<std::vector<uint8_t *> > & matrix , const std::vector<std::vector<std::pair<uint64_t, uint64_t> > > & pairs) const;

	/** Assemble all the kmers into sorted virtual superkmers.
	 * The algorithm garanty that the compaction is optimal (ie. it not possible to save more space 
	 * keeping the order). This compaction is not necessary the only one that is optimal.
	 * The matrix given in parameter will be sorted during this step.
	 * 
	 * @param Matrix of kmer positions in the memory buffer. There are k-m+1 vectors in the matrix.
	 * Each of this vectors contains all the kmers that share the same minimizer position.
	 * 
	 * @return Each pair of linked kmers. If a kmer is not linked one of the two values is a nullpointer.
	 **/
	std::vector<std::vector<uint8_t *> > sorted_assembly(std::vector<std::vector<uint8_t *> > & positions);
	/** Assemble all the kmers into virtual superkmers.
	 * The algorithm garanty that the compaction is optimal (ie. it not possible to save more space).
	 * This compaction is not necessary the only one that is optimal.
	 * 
	 * @param Matrix of kmer positions in the memory buffer. There are k-m+1 vectors in the matrix.
	 * Each of this vectors contains all the kmers that share the same minimizer position.
	 * 
	 * @return Each pair of linked kmers. If a kmer is not linked one of the two values is a nullpointer.
	 **/
	std::vector<std::pair<uint8_t *, uint8_t *> >  greedy_assembly(std::vector<std::vector<uint8_t *> > & kmers);


	void cli_prepare(CLI::App * subapp);
	/** Read a Section_Raw and write a bucketized and compacted file of the kmers.
	  * @param insection Section to bucketize then compact.
	  * @param prefix Prefix of the output file.
	  *
	  * @return Name of the file containing the result.
	  */
	void exec();
	void compact_section(Section_Minimizer & ism, Kff_file & outfile);

};

#endif