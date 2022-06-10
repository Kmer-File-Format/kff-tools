#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>

#include "CLI/CLI.hpp"
#include "kfftools.hpp"


#ifndef BUCKET_H
#define BUCKET_H

class Bucket: public KffTool {
private:
	std::string input_filename;
	std::string output_filename;

  uint8_t * kmer_buffer;
	uint8_t * data_buffer;
	uint8_t * copy_buffer;

	uint complement[4];

public:
	uint m;
	bool singleside;

	Bucket(uint8_t m=2, bool revcomp=true);
	~Bucket();
	void cli_prepare(CLI::App * subapp);
	void exec();

	void bucketize(uint8_t * kmer, uint8_t * data, uint64_t k);
	uint get_minimizer(const uint8_t * kmer, const uint64_t k, int & position) const;
};

#endif