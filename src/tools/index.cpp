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
  while (infile.tellp() != infile.end_position) {
    char section_type = infile.read_section_type();
    
    Section * current;
    // Read v sections
    if (section_type == 'v') {
      Section_GV * vsec = new Section_GV(&infile);
      current = vsec;
      // Remove footer variables from the v sections
      if (vsec->vars.find("footer_size") != vsec->vars.end())
        vsec->vars.erase("footer_size");
      if (vsec->vars.find("first_index") != vsec->vars.end())
        vsec->vars.erase("first_index");
    }
    // Read other type of section
    else {
      current = SectionBuilder::build(&infile);
    }
    
    // Copy only non i sections
    if (section_type != 'i') {
      current->copy(&outfile);
    }

    // Free memory from previous section
    delete current;
  }

  infile.close();
  outfile.close();
}