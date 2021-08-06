#include <string>
#include <iostream>
#include <vector>

#include "CLI11.hpp"
#include "kfftools.hpp"


#ifndef DATARM_H
#define DATARM_H

class DataRm: public KffTool {
private:
	std::string input_filename;
	std::string output_filename;

public:
	DataRm();
	void cli_prepare(CLI::App * subapp);
	void exec();
};

#endif