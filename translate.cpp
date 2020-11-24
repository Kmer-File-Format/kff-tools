#include <vector>
#include <string>

#include "translate.hpp"
#include "encoding.hpp"


using namespace std;


Translate::Translate() {
	input_filename = "";
	output_filename = "";
	encoding_str = "";
}


struct EncodingValidator : public CLI::Validator {
    EncodingValidator() {
        name_ = "ENCODING";
        func_ = [](const std::string &str) {
        	if(str.length() != 4)
            return string("The encoding must contains the 4 distinct letters A C G T.");
          else if (str.find("A") == string::npos or str.find("C") == string::npos or str.find("G") == string::npos or str.find("T") == string::npos)
          	return string("The encoding must contains the 4 distinct letters A C G T.");
          else
            return string();
        };
    }
};
const static EncodingValidator Lowercase;

void Translate::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("translate", "Translate a kff file from its encoding to another one.");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "The file to convert");
	input_option->required();

	CLI::Option * out_option = subapp->add_option("-o, --outfile", output_filename, "Translated output file.");
	out_option->required();

	CLI::Option * encoding = subapp->add_option("-e, --encoding", encoding_str, "A 4 letter string representing the encoding. For example AGTC represent the encoding where A=0, C=3, G=1, T=2.");
	encoding->required();
	encoding->check(EncodingValidator());
}

void Translate::exec() {
	// // Useful variables
	// long buffer_size = 1048576; // 1 MB
	// char buffer[1048576];
	// uint8_t global_encoding[4];

	// Prepare the destination encoding
	uint8_t dest_encoding[4];
	for (uint char_idx=0 ; char_idx<4 ; char_idx++) {
		char n = encoding_str[char_idx];
		switch(n) {
		case 'A':
			dest_encoding[0] = char_idx;
		break;
		case 'C':
			dest_encoding[1] = char_idx;
		break;
		case 'G':
			dest_encoding[2] = char_idx;
		break;
		case 'T':
			dest_encoding[3] = char_idx;
		break;
		}
	}

	// Read the encoding and prepare the translator
	Kff_file infile(input_filename, "r");
	infile.read_encoding();
	Translator translator(infile.encoding, dest_encoding);
	// translator.translate(uint8_t * sequence, size_t byte_length);

	// Write header of the output
	Kff_file outfile(output_filename, "w");
	outfile.write_encoding(
		dest_encoding[0],
		dest_encoding[1],
		dest_encoding[2],
		dest_encoding[3]
	);
	// Set metadata
	auto metadata_size = infile.size_metadata();
	uint8_t metadata[1024];
	infile.read_metadata(metadata_size, metadata);
	outfile.write_metadata(metadata_size, metadata);


	// // Append each file one by one
	// for (string in_filename : input_filenames) {
	// 	// Open the file
	// 	Kff_file infile(in_filename, "r");
	// 	infile.read_encoding();

	// 	// Encoding verification
	// 	for (uint i=0 ; i<4 ; i++) {
	// 		if (infile.encoding[i] != global_encoding[i]) {
	// 			cerr << "Wrong encoding for file " << in_filename << endl;
	// 			cerr << "Its nucleotide encoding is different from previous kff files." << endl;
	// 			cerr << "Please first use 'kff-tools translate' to have the same encoding" << endl;
	// 			exit(1);
	// 		}
	// 	}

	// 	// Jump over metadata
	// 	long metadata_size = static_cast<long>(infile.size_metadata());
	// 	infile.fs.seekp(infile.fs.tellp() + metadata_size);

	// 	// Read section by section
	// 	char section_type = infile.read_section_type();
	// 	while(not infile.fs.eof()) {
	// 		vector<string> to_copy;
	// 		long size, end_byte, begin_byte;

	// 		switch (section_type) {
	// 			// Write the variables that change from previous sections (possibly sections from other input files)
	// 			case 'v':
	// 			// Read the values
	// 			infile.open_section_GV();

	// 			// Verify the presence and value of each variable in output
	// 			for (auto& tuple : infile.global_vars) {
	// 				if (outfile.global_vars.find(tuple.first) == outfile.global_vars.end()
	// 						or outfile.global_vars[tuple.first] != tuple.second)
	// 					to_copy.push_back(tuple.first);
	// 			}
	// 			// Create a global variable section if needed
	// 			if (to_copy.size() > 0) {
	// 				Section_GV sgv = outfile.open_section_GV();
	// 				// Write variables
	// 				for (string s : to_copy)
	// 					sgv.write_var(s, infile.global_vars[s]);
	// 				sgv.close();
	// 			}
	// 			break;

	// 			// copy the sequence section from input to output
	// 			case 'r':
	// 			case 'm':
	// 			// Analyse the section size
	// 			begin_byte = infile.fs.tellp();
	// 			infile.jump_next_section();
	// 			end_byte = infile.fs.tellp();
	// 			size = end_byte - begin_byte;
	// 			infile.fs.seekp(begin_byte);

	// 			// Read from input and write into output
	// 			while (size > 0) {
	// 				size_t size_to_copy = size > buffer_size ? buffer_size : size;

	// 				infile.fs.read(buffer, size_to_copy);
	// 				outfile.fs.write(buffer, size_to_copy);

	// 				size -= size_to_copy;
	// 			}
	// 			break;

	// 			default:
	// 				cerr << "Unknown section type " << section_type << endl;
	// 				exit(2);
	// 		}

	// 		// Prepare next section
	// 		section_type = infile.read_section_type();
	// 	}

	// 	infile.close();
	// }

	// outfile.close();
}
