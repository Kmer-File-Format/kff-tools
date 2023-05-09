#include <string>
#include <iostream>
#include <vector>

#include "CLI/CLI.hpp"
#include "kfftools.hpp"


#ifndef DISJOIN_H
#define DISJOIN_H

class Disjoin: public KffTool {
private:
	std::string input_filename;
	std::string output_filename;

public:
	Disjoin();
	void cli_prepare(CLI::App * subapp);
	void exec();
};

#endif