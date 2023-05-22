#include <string>
#include <iostream>
#include <vector>

#include "CLI/CLI.hpp"
#include "kfftools.hpp"


#ifndef SHUFFLE_H
#define SHUFFLE_H

class Shuffle: public KffTool {
private:
	std::string input_filename;
	std::string output_filename;

public:
	Shuffle();
	void cli_prepare(CLI::App * subapp);
	void shuffle(std::string input, std::string output);
	void exec();
};

#endif
