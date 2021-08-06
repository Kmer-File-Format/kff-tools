// #include <experimental/filesystem>
#include <string>
#include <iostream>

#include "CLI11.hpp"
#include "kfftools.hpp"


#ifndef SPLIT_H
#define SPLIT_H

class Split: public KffTool {
private:
	std::string input_filename;
	std::string output_dirname;

public:
	Split();
	void cli_prepare(CLI::App * subapp);
	void exec();
};


#endif