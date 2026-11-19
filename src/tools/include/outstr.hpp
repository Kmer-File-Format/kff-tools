#include <cstdint>
#include <string>
#include <iostream>
#include <vector>

#include "CLI/CLI.hpp"
#include "kfftools.hpp"


#ifndef OUTSTR_H
#define OUTSTR_H

/**
 * Writes the data of a kmer as outstr prints it:
 * - nothing when data_size is 0;
 * - the decimal value of the big endian bytes when data_size < 8 (computed on
 *   32 bits, as unsigned int);
 * - "[b0][b1]..." with the decimal value of each byte otherwise.
 *
 * @param out Destination of at most outstr_data_chars(data_size) characters.
 * @return The position after the last written character.
 */
char * outstr_write_data(const uint8_t * data, size_t data_size, char * out);
/** Maximal number of characters written by outstr_write_data */
size_t outstr_data_chars(size_t data_size);

class Outstr: public KffTool {
private:
	std::string input_filename;
	bool revcomp;

public:
	Outstr();
	void cli_prepare(CLI::App * subapp);
	void exec();
};

#endif