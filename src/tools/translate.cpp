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
	input_option->check(CLI::ExistingFile);

	CLI::Option * out_option = subapp->add_option("-o, --outfile", output_filename, "Translated output file.");
	out_option->required();

	CLI::Option * encoding = subapp->add_option("-e, --encoding", encoding_str, "A 4 letter string representing the encoding. For example AGTC represent the encoding where A=0, C=3, G=1, T=2.");
	encoding->required();
	encoding->check(EncodingValidator());
}

void Translate::exec() {
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

	const long buffer_size = 1048576; // 1 MB
	uint8_t buffer[1048576];

	// Read the encoding and prepare the translator
	Kff_file infile(input_filename, "r");
	Translator translator(infile.encoding, dest_encoding);

	// Write header of the output
	Kff_file outfile(output_filename, "w");
	outfile.set_indexation(false);
	outfile.write_encoding(
		dest_encoding[0],
		dest_encoding[1],
		dest_encoding[2],
		dest_encoding[3]
	);
	// Set metadata
	uint8_t * metadata = new uint8_t[infile.metadata_size];
	infile.read_metadata(metadata);
	outfile.write_metadata(infile.metadata_size, metadata);
	delete[] metadata;

	// Prepare sequence and data buffers
	uint8_t * nucleotides = new uint8_t[1];
	uint8_t * data = new uint8_t[1];

	// Read and write section per section
	char section_type = infile.read_section_type();
	while (infile.tellp() != infile.end_position) {
		// Read variables
		if (section_type == 'v') {
			// Load variables
			Section_GV isgv(&infile);
			Section_GV osgv(&outfile);

			bool nucl_buffer_changed = false;
			bool data_buffer_changed = false;
			for (auto var_tuple : isgv.vars) {
				osgv.write_var(var_tuple.first, var_tuple.second);

				if (var_tuple.first == "k" or var_tuple.first == "max")
					nucl_buffer_changed = true;
				if (var_tuple.first == "data_size" or var_tuple.first == "max")
					data_buffer_changed = true;
			}
			isgv.close();
			osgv.close();

			// Buffer updates
			if (nucl_buffer_changed) {
				delete[] nucleotides;
				uint max_nucl = outfile.global_vars["k"] + outfile.global_vars["max"] - 1;
				nucleotides = new uint8_t[max_nucl / 4 + 1];
			}
			if (data_buffer_changed) {
				delete[] data;
				uint nb_data = outfile.global_vars["max"];
				data = new uint8_t[nb_data * outfile.global_vars["data_size"]];
			}
		}
		// Pure copy
		else if (section_type == 'i') {
			// Read the section
			Section_Index si(&infile);
			si.close();
			// Compute size and go back to start
			long size = infile.tellp() - si.beginning;
			infile.jump(-size);
			// copy
			while (size > 0) {
				size_t size_to_copy = size > buffer_size ? buffer_size : size;

				infile.read(buffer, size_to_copy);
				outfile.write(buffer, size_to_copy);

				size -= size_to_copy;
			}
		}
		// translate a raw block
		else if (section_type == 'r') {
			// Open sections
			Section_Raw in_section(&infile);
			Section_Raw out_section(&outfile);

			// Translate block per block
			uint64_t k = outfile.global_vars["k"];
			for (uint i=0 ; i<in_section.nb_blocks ; i++) {
				uint nb_kmers = in_section.read_compacted_sequence(nucleotides, data);
				uint size = k + nb_kmers - 1;
				size = size % 4 == 0 ? size / 4 : size / 4 + 1;
				translator.translate(nucleotides, size);
				out_section.write_compacted_sequence(nucleotides, k + nb_kmers - 1, data);
			}

			in_section.close();
			out_section.close();
		}
		// Translate a minimizer block
		else if (section_type == 'm') {
			// Open sections
			Section_Minimizer in_section(&infile);
			Section_Minimizer out_section(&outfile);

			// translate and write the minimizer
			uint k = outfile.global_vars["k"];
			uint m = outfile.global_vars["m"];
			uint size = m % 4 == 0 ? m / 4 : m / 4 + 1;
			translator.translate(in_section.minimizer, size);
			out_section.write_minimizer(in_section.minimizer);

			// Translate block per block
			for (uint i=0 ; i<in_section.nb_blocks ; i++) {
				// Read
				uint64_t mini_pos;
				uint nb_kmers = in_section.read_compacted_sequence_without_mini(nucleotides, data, mini_pos);
				// Translate
				uint nucl_size = k - m + nb_kmers - 1;
				uint byte_size = nucl_size % 4 == 0 ? nucl_size / 4 : nucl_size / 4 + 1;
				translator.translate(nucleotides, byte_size);
				// Write
				out_section.write_compacted_sequence_without_mini(nucleotides, nucl_size, mini_pos, data);
			}

			in_section.close();
			out_section.close();
		} else {
			cerr << infile.tellp() << ": Unknown section " << section_type << endl;
			exit(1);
		}

		section_type = infile.read_section_type();
	}

	delete[] nucleotides;
	delete[] data;
	infile.close();
	outfile.close();
}
