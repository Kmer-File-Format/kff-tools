#include <string>
#include <iostream>
#include <vector>

#include "CLI11.hpp"
#include "kfftools.hpp"
#include "kmer_file_format/C++/kff_io.hpp"


#ifndef MERGE_H
#define MERGE_H

class Merge: public KffTool {
private:
	std::vector<std::string> input_filenames;
	std::string output_filename;

public:
	Merge();
	void cli_prepare(CLI::App * subapp);
	void exec();
};

#endif