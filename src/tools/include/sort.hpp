#include <string>
#include <iostream>
#include <vector>

#include "CLI/CLI.hpp"
#include "kfftools.hpp"


#ifndef SORT_H
#define SORT_H

class Sort: public KffTool {
private:
	std::string input_filename;
	std::string output_filename;

public:
	Sort();
	void cli_prepare(CLI::App * subapp);
	void sort(std::string input, std::string output);
	void exec();
};

#endif
