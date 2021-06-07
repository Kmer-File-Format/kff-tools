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

	bool split;
	uint m;

	uint load_mem_size;
	uint8_t * loading_memory;
	std::vector<uint> kmer_nbs;
  std::vector<uint> mini_pos;

  uint8_t * kmer_buffer;
  uint8_t * skmer_buffer;
	uint8_t * data_buffer;
	
	std::unordered_map<std::string, uint64_t> saved_variables;

	void write_variables(std::unordered_map<std::string, uint64_t> & variables, Kff_file & file);

	void loadSectionBlocks(Section_Minimizer & ms, Kff_file & infile);
	std::vector<std::vector<uint> > link_kmers(uint nb_kmers, Kff_file & infile);
	void compact_and_save(std::vector<std::vector<uint> > paths, Kff_file & outfile, uint8_t * minimizer, uint input_max_kmers);

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
	std::string bucketize(Kff_file & infile, std::string & prefix, uint m);
	void compact(std::string input, std::string output);
	void exec();
};

#endif