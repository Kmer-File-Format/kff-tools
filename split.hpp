// #include <experimental/filesystem>
#include <string>
#include <iostream>

#include "CLI11.hpp"
#include "kfftools.hpp"
#include "kmer_file_format/C++/kff_io.hpp"


#ifndef SPLIT_H
#define SPLIT_H

class Split: public KffTool {
private:
	std::string input_filename;
	std::string output_dirname;
	
	std::ifstream input_file;
	std::istream * in;
public:
	Split();
	void cli_prepare(CLI::App * subapp);
	void init();
};


#endif