#include <iostream>
#include "CLI11.hpp"


using namespace std;


int parse_args(int argc, char** argv) {
	CLI::App app{"kff-tools is a software for kff file manipulations. For more details on kff format, please refer to https://github.com/yoann-dufresne/kmer_file_format"};
	app.require_subcommand(1);
	CLI::Option * help =	app.get_help_ptr();

	try {
    app.parse(argc, argv);
	} catch (const CLI::ParseError &e) {
		if (!help->empty()) {
			exit(0);
		}
    return app.exit(e);
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
