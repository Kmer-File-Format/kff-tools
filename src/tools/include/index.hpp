#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>

#include "CLI/CLI.hpp"
#include "kfftools.hpp"


#ifndef INDEX_H
#define INDEX_H

class Index: public KffTool {
private:
	std::string input_filename;
	std::string output_filename;

public:
	Index();

	void cli_prepare(CLI::App * subapp);
	void exec();

};

#endif