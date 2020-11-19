#include <iostream>

#include "CLI11.hpp"
#include "split.hpp"


using namespace std;


int parse_args(int argc, char** argv) {
	// Main command
	CLI::App app{"kff-tools is a software for kff file manipulations. For more details on kff format, please refer to https://github.com/yoann-dufresne/kmer_file_format"};
	app.require_subcommand(1);
	CLI::Option * help =	app.get_help_ptr();

	// Define the subcommands
	CLI::App * split_app = app.add_subcommand("split", "Split a kff file into one file per section.");
	Split split;
	split.cli_prepare(split_app);

	// Parsing and return status
	try {
    app.parse(argc, argv);
	} catch (const CLI::ParseError &e) {
    auto val = app.exit(e);
		if (!help->empty()) {
			exit(0);
		} else if (val != 0) {
			exit(val);
		}
	}

	return 0;
}


int main(int argc, char** argv) {
	int ret_val;
	if ((ret_val = parse_args(argc, argv)) != 0) {
		cerr << "Error during parsing" << endl;
		return ret_val;
	}

	return 0;
}
