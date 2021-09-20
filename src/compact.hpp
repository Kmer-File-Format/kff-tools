#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>

#include "CLI11.hpp"
#include "kfftools.hpp"


#ifndef COMPACT_H
#define COMPACT_H

class Compact: public KffTool {
private:
	std::string input_filename;
	std::string output_filename;

	uint8_t * kmer_buffer;
	uint64_t buffer_size;
	uint64_t next_free;

	void write_paths(const std::vector<std::vector<uint8_t *> > & paths, Section_Minimizer & sm, const uint k, const uint m, const uint data_size);
	std::vector<std::pair<uint8_t *, uint8_t *> >  greedy_assembly(uint8_t * kmer_buffer, std::vector<std::vector<long> > & positions, const uint k, const uint m);
	std::vector<std::vector<uint8_t *> > pairs_to_paths(const std::vector<std::pair<uint8_t *, uint8_t *> > & to_compact);

public:
	Compact();
	~Compact();
	void cli_prepare(CLI::App * subapp);
	/** Read a Section_Raw and write a bucketized and compacted file of the kmers.
	  * @param insection Section to bucketize then compact.
	  * @param prefix Prefix of the output file.
	  *
	  * @return Name of the file containing the result.
	  */
	void exec();
	void compact_section(Section_Minimizer & ism, Kff_file & outfile);

	std::vector<std::vector<long> > prepare_kmer_matrix(Section_Minimizer & sm);
};

#endif