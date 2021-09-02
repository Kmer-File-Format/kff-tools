#include <string>
#include <iostream>
#include <vector>

#include "CLI11.hpp"
#include "kfftools.hpp"


#ifndef MERGE_H
#define MERGE_H

class Merge: public KffTool {
private:
	std::vector<std::string> input_filenames;
	std::string output_filename;

public:
	Merge();
	void cli_prepare(CLI::App * subapp);
	void merge(const std::vector<std::string> inputs, std::string output);
	void merge(const std::vector<Kff_file *> & inputs, std::string output);
	void exec();
};

#endif