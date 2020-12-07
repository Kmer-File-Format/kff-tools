#include <string>
#include <iostream>
#include <vector>

#include "CLI11.hpp"
#include "kfftools.hpp"


#ifndef TRANSLATE_H
#define TRANSLATE_H

class Translate: public KffTool {
private:
	std::string input_filename;
	std::string output_filename;
	std::string encoding_str;

public:
	Translate();
	void cli_prepare(CLI::App * subapp);
	void exec();
};

#endif