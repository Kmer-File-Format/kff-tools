#include "index.hpp"

using namespace std;


Index::Index() {
  input_filename = "";
  output_filename = "";
}


void Index::cli_prepare(CLI::App * app) {
  this->subapp = app->add_subcommand("index", "Copy a kff file and add an index of the sections");
  CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "Input kff file to index.");
  input_option->required();
  input_option->check(CLI::ExistingFile);
  CLI::Option * out_option = subapp->add_option("-o, --outfile", output_filename, "Indexed kff to write (must be different from the input)");
  out_option->required();
}


void Index::exec() {
  Kff_file infile(this->input_filename, "r");
  Kff_file outfile(this->output_filename, "w");

  // Set encoding
  outfile.write_encoding(infile.encoding);

  // Set flags
  outfile.set_uniqueness(infile.uniqueness);
  outfile.set_canonicity(infile.canonicity);
  outfile.set_indexation(true);

  // Copy metadata
  uint8_t * meta = new uint8_t[infile.metadata_size];
  infile.read_metadata(meta);
  outfile.write_metadata(infile.metadata_size, meta);
  delete[] meta;

  // Copy section by section
  char section_type = infile.read_section_type();
  while (infile.tellp() != infile.end_position) {
    if (section_type != 'i') {
      Section * current = SectionBuilder::build(&infile);
      current->copy(&outfile);
      delete current;
    }
  }

  infile.close();
  outfile.close();
}