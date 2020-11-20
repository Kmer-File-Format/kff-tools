#include <iostream>

#include "kfftools.hpp"
#include "CLI11.hpp"
#include "split.hpp"


using namespace std;


KffTool * parse_args(int argc, char** argv) {
	// Main command
	CLI::App app{"kff-tools is a software for kff file manipulations. For more details on kff format, please refer to https://github.com/yoann-dufresne/kmer_file_format"};
	app.require_subcommand(1);
	CLI::Option * help =	app.get_help_ptr();

	// Define the subcommands
	CLI::App * split_app = app.add_subcommand("split", "Split a kff file into one file per section.");
	KffTool * split = new Split();
	split->cli_prepare(split_app);

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

	// Determine wich app was parsed
	KffTool * tool = nullptr;
	CLI::App * parsed_app = nullptr;
	if (split_app->parsed()) {
		tool = split;
		parsed_app = split_app;
	}

	// Sub-app help function
	if (not parsed_app->get_help_ptr()->empty()) {
		exit(0);
	}

	return tool;
}


int main(int argc, char** argv) {
	KffTool * tool = parse_args(argc, argv);
	tool->exec();

	return 0;
}
