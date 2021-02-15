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
	std::string output_prefix;
	uint m;
	bool verbose;

	void search_mini(uint8_t * bin, const uint k, uint & minimizer, uint & minimizer_position);
	void monofile();
	void multifile();

public:
	Convert();
	void cli_prepare(CLI::App * subapp);
	void exec();
};

#endif