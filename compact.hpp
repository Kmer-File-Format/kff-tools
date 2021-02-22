#include <string>
#include <iostream>
#include <vector>

#include "CLI11.hpp"
#include "kfftools.hpp"


#ifndef COMPACT_H
#define COMPACT_H

class Compact: public KffTool {
private:
	std::string input_filename;
	std::string output_filename;

	uint load_mem_size;
	uint8_t * loading_memory;
	std::vector<uint> kmer_nbs;
  std::vector<uint> mini_pos;

	void loadSectionBlocks(Section_Minimizer & ms, Kff_file & infile);
	std::vector<std::vector<uint> > link_kmers(uint nb_kmers, Kff_file & infile);

public:
	Compact();
	void cli_prepare(CLI::App * subapp);
	void compact(std::string input, std::string output);
	void exec();
};

#endif