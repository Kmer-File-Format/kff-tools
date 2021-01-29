#include <string>
#include <iostream>
#include <vector>

#include "CLI11.hpp"
#include "kfftools.hpp"


#ifndef CONVERT_H
#define CONVERT_H

class Convert: public KffTool {
private:
	std::string input_filename;
	std::string output_filename;
	bool verbose;

public:
	Convert();
	void cli_prepare(CLI::App * subapp);
	void exec();
};

#endif