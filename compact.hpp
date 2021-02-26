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

  uint8_t * kmer_buffer;
  uint8_t * skmer_buffer;
	uint8_t * data_buffer;

	void loadSectionBlocks(Section_Minimizer & ms, Kff_file & infile);
	std::vector<std::vector<uint> > link_kmers(uint nb_kmers, Kff_file & infile);
	void compact_and_save(std::vector<std::vector<uint> > paths, Kff_file & outfile, uint8_t * minimizer, uint input_max_kmers);

public:
	Compact();
	~Compact();
	void cli_prepare(CLI::App * subapp);
	void compact(std::string input, std::string output);
	void exec();
};

#endif