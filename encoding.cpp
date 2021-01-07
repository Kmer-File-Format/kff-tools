#include "encoding.hpp"

using namespace std;



Translator::Translator(uint8_t source[4], uint8_t destination[4]) {
	// Nucleotide translation values
	uint8_t nucl_translation[4];
	for (uint i=0 ; i<4 ; i++)
		nucl_translation[source[i]] = destination[i] & 0b11;
	
	// Lookup translation for Bytes
	for (uint i=0 ; i<256 ; i++) {
		lookup[i] = 0;

		for (uint pos=0 ; pos<4 ; pos++) {
			// Get nucleotide
			uint8_t letter = (i >> (2*pos)) & 0b11;
			// Change encoding
			letter = nucl_translation[letter];
			// Writ in the lookup table
			lookup[i] = lookup[i] | (letter << (2*pos));
		}
	}
}

void Translator::translate(uint8_t * sequence, size_t byte_length) {
	// Translate Byte per Byte
	for (size_t idx=0 ; idx<byte_length ; idx++) {
		// Translate a Byte
		sequence[idx] = lookup[sequence[idx]];
	}
}


Stringifyer::Stringifyer(uint8_t encoding[4]) {
	// Nucleotide translation values
	string nucl_translation[4];
	nucl_translation[encoding[0]] = "A";
	nucl_translation[encoding[1]] = "C";
	nucl_translation[encoding[2]] = "G";
	nucl_translation[encoding[3]] = "T";
	
	// Lookup translation for Bytes
	for (uint i=0 ; i<256 ; i++) {
		lookup[i] = "";

		for (uint pos=0 ; pos<4 ; pos++) {
			// Get nucleotide
			uint8_t letter = (i >> (2*pos)) & 0b11;
			// Write in the lookup table
			lookup[i] = nucl_translation[letter] + lookup[i];
		}
	}
}


string Stringifyer::translate(uint8_t * sequence, size_t nucl_length) {
	string result = "";
	uint byte_length = nucl_length % 4 == 0 ? nucl_length / 4 : nucl_length / 4 + 1;

	// Prefix can be truncated
	result = lookup[sequence[0]];
	if (nucl_length % 4 != 0) {
		uint to_remove = 4 - (nucl_length % 4);
		result = result.substr(to_remove);
	}

	// Translate Byte per Byte
	for (size_t idx=1 ; idx<byte_length ; idx++) {
		// Translate a Byte
		uint8_t byte = sequence[idx];
		result += lookup[byte];
	}

	return result;
}