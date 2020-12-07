#include <iostream>
#include <vector>

#include "kfftools.hpp"
#include "CLI11.hpp"

#include "split.hpp"
#include "merge.hpp"
#include "translate.hpp"
#include "outstr.hpp"
#include "datarm.hpp"


using namespace std;


KffTool * parse_args(int argc, char** argv, vector<KffTool *> tools) {
	// Main command
	CLI::App app{"kff-tools is a software for kff file manipulations. For more details on kff format, please refer to https://github.com/Kmer-File-Format/kff-reference"};
	app.require_subcommand(1);
	CLI::Option * help =	app.get_help_ptr();

	// Subcommands prepare
	for (KffTool * tool : tools) {
		tool->cli_prepare(&app);
	}

	// Parsing and return status if wrong
	try {
    app.parse(argc, argv);
	} catch (const CLI::ParseError &e) {
    auto val = app.exit(e);
		if (val != 0) {
			exit(val);
		}
	}

	// Help detection
	if (!help->empty()) {
		return nullptr;
	}

	// Read the command line return
	for (KffTool * tool : tools) {
		if (tool->subapp->parsed()) {
			// Help on tool triggered
			if (not tool->subapp->get_help_ptr()->empty()) {
				return nullptr;
			} else {
				return tool;
			}
		}
	}

	return nullptr;
}


int main(int argc, char** argv) {
	// Prepare tools
	vector<KffTool *> tools;
	tools.push_back(new Split());
	tools.push_back(new Merge());
	tools.push_back(new Translate());
	tools.push_back(new Outstr());
	tools.push_back(new DataRm());

	// Get the one selected
	KffTool * tool = parse_args(argc, argv, tools);

	if (tool != nullptr)
		tool->exec();

	for (KffTool * tool : tools)
		delete tool;

	return 0;
}
