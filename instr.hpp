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
	std::string output_prefix;
	uint m;
	uint data_size;
	bool split;

	void search_mini(uint8_t * bin, const uint k, uint & minimizer, uint & minimizer_position);
	void monofile();
	void multifile();

public:
	Instr();
	void cli_prepare(CLI::App * subapp);
	void exec();
};

#endif