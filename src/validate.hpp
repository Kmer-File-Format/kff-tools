#include <string>
#include <iostream>
#include <vector>

#include "CLI11.hpp"
#include "kfftools.hpp"
#include "encoding.hpp"


#ifndef VALID_H
#define VALID_H

class Validate: public KffTool {
private:
	std::string input_filename;
	bool index_only;
	bool verbose;
	Stringifyer strif;

public:
	Validate();
	void cli_prepare(CLI::App * subapp);
	void exec();

	bool is_valid_r_section(Kff_file & infile);
	bool is_valid_m_section(Kff_file & infile);
	bool is_valid_i_section(Kff_file & infile);
};

#endif