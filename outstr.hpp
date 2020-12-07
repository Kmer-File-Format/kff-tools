#include <string>
#include <iostream>
#include <vector>

#include "CLI11.hpp"
#include "kfftools.hpp"


#ifndef OUTSTR_H
#define OUTSTR_H

class Outstr: public KffTool {
private:
	std::string input_filename;

public:
	Outstr();
	void cli_prepare(CLI::App * subapp);
	void exec();
};

#endif