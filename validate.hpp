#include <string>
#include <iostream>
#include <vector>

#include "CLI11.hpp"
#include "kfftools.hpp"


#ifndef VALID_H
#define VALID_H

class Validate: public KffTool {
private:
	std::string input_filename;
	bool verbose;

public:
	Validate();
	void cli_prepare(CLI::App * subapp);
	void exec();
};

#endif