#include <string>
#include <iostream>
#include <vector>

#include "CLI11.hpp"
#include "kfftools.hpp"


#ifndef INSTR_H
#define INSTR_H

class Instr: public KffTool {
private:
	std::string input_filename;
	std::string output_filename;

	uint data_size;
	uint k;
	uint max_kmerseq;

	char delimiter;
  char data_delimiter;

	void monofile();
	void multifile();

public:
	Instr();
	void cli_prepare(CLI::App * subapp);
	void exec();
};

#endif