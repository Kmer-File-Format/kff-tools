#include "split.hpp"


using namespace std;


Split::Split()
: in (& std::cin)
{
	// Paths
	this->input_filename = "";
	this->output_dirname = "./";
}

void Split::cli_prepare(CLI::App * subapp) {
	subapp->add_option("-i, --input", input_filename, "Input file in kff format");
	subapp->add_option("-o, --outdir", output_dirname, "Output directory where to put all the subkff files");
}

void Split::init() {
	// Open streams stream
	if (input_filename != "") {
    input_file.open(input_filename, std::fstream::in);
		in = static_cast<std::istream*>(&input_file);;
	}
}
