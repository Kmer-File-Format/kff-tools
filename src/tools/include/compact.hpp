#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>

#include "CLI/CLI.hpp"
#include "kfftools.hpp"
#include "blocks.hpp"


#ifndef COMPACT_H
#define COMPACT_H

/** A minimizer section loaded in memory, then the bytes of its compacted version */
struct CompactSection {
	// Variables of the input section
	uint64_t k, m, max, data_size;
	// Variables of the output file when the compacted section is written
	uint64_t out_k, out_m, out_max, out_data_size;
	// Bytes of the minimizer, as read
	std::vector<uint8_t> minimizer;
	// Blocks: number of kmers, minimizer position and position of the sequence (followed by
	// the data) in bytes
	struct Block {
		uint64_t nb_kmers;
		uint64_t mini_pos;
		uint64_t offset;
	};
	std::vector<Block> blocks;
	std::vector<uint8_t> bytes;
	// Bytes of the compacted section: type, minimizer, number of blocks and blocks
	std::vector<uint8_t> output;
};

/** Superkmers as paths of kmers of the buffer, stored one after the other: the kmers of the path
 * i are kmers[ends[i-1]] to kmers[ends[i] - 1] (ends[-1] = 0). */
struct Paths {
	std::vector<uint8_t *> kmers;
	std::vector<uint64_t> ends;

	size_t size() const { return ends.size(); }
	void clear() { kmers.clear(); ends.clear(); }
	/** Ends the current path (the kmers added since the end of the previous one) */
	void close_path() { ends.push_back(kmers.size()); }
};

/** Symbols of chain_columns(): the k-m-1 nucleotides on 64 bits (up to 32 nucleotides) or 128
 * bits, with a position */
struct ChainSymbol64 {
	uint64_t lo;
	uint64_t pos;
	bool same(const ChainSymbol64 & o) const { return lo == o.lo; }
	bool operator<(const ChainSymbol64 & o) const { return lo != o.lo ? lo < o.lo : pos < o.pos; }
	static ChainSymbol64 make(uint64_t, uint64_t lo, uint64_t pos) { return ChainSymbol64{lo, pos}; }
	uint64_t hash() const { return lo * 0x9E3779B97F4A7C15ULL; }
};

struct ChainSymbol128 {
	uint64_t hi;
	uint64_t lo;
	uint64_t pos;
	bool same(const ChainSymbol128 & o) const { return hi == o.hi and lo == o.lo; }
	bool operator<(const ChainSymbol128 & o) const {
		if (hi != o.hi) return hi < o.hi;
		if (lo != o.lo) return lo < o.lo;
		return pos < o.pos;
	}
	static ChainSymbol128 make(uint64_t hi, uint64_t lo, uint64_t pos) { return ChainSymbol128{hi, lo, pos}; }
	uint64_t hash() const { return (lo ^ (hi * 0xBF58476D1CE4E5B9ULL)) * 0x9E3779B97F4A7C15ULL; }
};

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
	// Number of sections compacted in parallel
	uint nb_threads;

	uint8_t * kmer_buffer;
	uint64_t buffer_size;
	uint64_t next_free;

	// Sort keys of a column (sort_matrix(), up to 32 nucleotides)
	struct Keyed64 {
		uint64_t key;
		uint8_t * kmer;
	};
	std::vector<Keyed64> sort_keys, sort_keys_tmp;
	// Keys of the columns sorted by the last sort_matrix() (valid up to 32 nucleotides)
	std::vector<std::vector<uint64_t> > column_keys;
	bool column_keys_valid = false;

	// Buffers of write_path(), kept from path to path
	std::vector<uint8_t> skmer_bytes;
	std::vector<uint8_t> skmer_data;
	// Buffers of sorted_assembly(), kept from section to section: candidate links of 2
	// columns, index of pair_kmers() and links of each column to the next one
	std::vector<std::pair<uint64_t, uint64_t> > candidates;
	std::vector<std::pair<uint64_t, uint64_t> > pair_index;
	std::vector<std::vector<std::pair<uint64_t, uint64_t> > > links;
	// Kmer matrix of the section and number of kmers per column, kept from section to section
	std::vector<std::vector<uint8_t *> > matrix;
	std::vector<uint64_t> column_counts;

    Compact(std::string input_filename, std::string output_filename, bool sorted = false);
	Compact();
	~Compact();


	/** Encode the blocks of a minimizer section containing all the superkmers given as paths.
	 * The process preserve the order of the superkmers.
	 * 
	 * @param paths List of all the kmer paths that represent virtual superkmers.
	 * @param encoder Encoder of the blocks of the output section (gives the data size).
	 * @param output The blocks are appended to this array.
	 * @return The number of blocks written.
	 **/
	uint64_t write_paths(const std::vector<std::vector<uint8_t *> > & paths, const MinimizerBlockEncoder & encoder, std::vector<uint8_t> & output);
	uint64_t write_paths(const Paths & paths, const MinimizerBlockEncoder & encoder, std::vector<uint8_t> & output);
	/** Encode the block of one path of size kmers */
	void write_path(uint8_t * const * path, size_t size, const MinimizerBlockEncoder & encoder, std::vector<uint8_t> & output);

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
	std::vector<std::vector<uint8_t *> > prepare_kmer_matrix(const CompactSection & section);
	/** Same, the matrix is filled in matrix (its columns are reused) */
	void prepare_kmer_matrix(const CompactSection & section, std::vector<std::vector<uint8_t *> > & matrix);

	/** Load the blocks of a minimizer section into memory.
	 * @param outfile The output file: its current variables are the ones of the compacted section.
	 **/
	void load_section(Section_Minimizer & ism, Kff_file & outfile, CompactSection & section) const;

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

	/** Key of a kmer of the buffer (k - m <= 64) for the order of interleaved_compare_kmers():
	 * the nucleotides sorted by distance to the minimizer (the suffix nucleotide first at each
	 * distance), 2 bits each, as a 128-bit integer (high, low).
	 **/
	void interleaved_key(const uint8_t * kmer, uint64_t & high, uint64_t & low) const;

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
	/** Same, the pairs are written in pairs (index: buffer of the function) */
	void pair_kmers(const std::vector<uint8_t *> & column1, const std::vector<uint8_t *> & column2,
	                std::vector<std::pair<uint64_t, uint64_t> > & pairs, std::vector<std::pair<uint64_t, uint64_t> > & index) const;

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
	std::vector<std::pair<uint64_t, uint64_t> > colinear_chaining(std::vector<std::pair<uint64_t, uint64_t> > & candidates) const;
	/** Same chain as colinear_chaining(pair_kmers(column1, column2)), computed without the
	 * candidate pairs (2 <= k - m <= 64): the work is about the size of the chains instead of the
	 * number of candidates (equal kmers make d x e candidates).
	 **/
	void chain_columns(const std::vector<uint8_t *> & column1, const std::vector<uint8_t *> & column2,
	                   std::vector<std::pair<uint64_t, uint64_t> > & chain);
	/** Buffers of chain_columns() */
	struct ChainChange {
		uint64_t row, column, level, prev;
	};
	struct ChainBuffers {
		// Groups of the second column (same symbol): group of each column, first column of each
		// group in columns, columns of the groups in increasing order, group of each row
		std::vector<uint64_t> group_of_column, group_start, group_columns, group_of_row, table_groups;
		std::vector<uint64_t> thresholds, level_head;
		std::vector<ChainChange> changes;
		std::vector<std::pair<uint64_t, uint64_t> > pending;
	};
	ChainBuffers chain_buffers;
	std::vector<ChainSymbol64> chain_xs64, chain_index64;
	std::vector<ChainSymbol128> chain_xs128, chain_index128;
	/** Same, the chain is written in chain (candidates can be modified) */
	void colinear_chaining(std::vector<std::pair<uint64_t, uint64_t> > & candidates, std::vector<std::pair<uint64_t, uint64_t> > & chain) const;

	/** From the kmer matrix and the list of all the colinear chained pairs of kmers, generate the ordered list of superkmers.
	 * 
	 * @param matrix The full kmer matrix 2D sorted
	 * @param colinear_chainings sorted pairs of kmers to assemble (number of matrix columns - 1 lists)
	 * @return The list of sorted linked kmers. They can be compacted in a skmers style.
	 **/
	std::vector<std::vector<uint8_t *> > polish_sort(const std::vector<std::vector<uint8_t *> > & matrix , const std::vector<std::vector<std::pair<uint64_t, uint64_t> > > & pairs) const;
	/** Buffers of polish_sort(), kept from section to section */
	struct PolishBuffers {
		std::vector<uint64_t> col_start, prev_kmer, sk_of, sk_first, sk_col, sk_size, sk_counts, head_sk;
		std::vector<uint8_t> head_ready;
		std::vector<size_t> current_kmers, candidate_order;
		std::vector<std::pair<uint64_t, uint64_t> > current_values;
	};
	PolishBuffers polish_buffers;
	/** Same, the superkmers are appended to paths (buffers of the object) */
	void polish_sort(const std::vector<std::vector<uint8_t *> > & matrix , const std::vector<std::vector<std::pair<uint64_t, uint64_t> > > & pairs, Paths & paths);
	void polish_sort(const std::vector<std::vector<uint8_t *> > & matrix , const std::vector<std::vector<std::pair<uint64_t, uint64_t> > > & pairs, Paths & paths, PolishBuffers & buffers) const;
	/** keys: interleaved values of the kmers (keys[col][row], 64 bits) instead of computing them, or
	 * nullptr. They are the keys of sort_matrix() when k - m <= 32 is a multiple of 4. */
	void polish_sort(const std::vector<std::vector<uint8_t *> > & matrix , const std::vector<std::vector<std::pair<uint64_t, uint64_t> > > & pairs, Paths & paths, PolishBuffers & buffers,
	                 const std::vector<std::vector<uint64_t> > * keys) const;

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
	void sorted_assembly(std::vector<std::vector<uint8_t *> > & positions, Paths & paths);
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
	/** Compact a section loaded in memory: fills section.output */
	void compact_section(CompactSection & section);

};


#endif