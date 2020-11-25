#include "encoding.hpp"


Translator::Translator(uint8_t source[4], uint8_t destination[4]) {
	// TODO: Need encoding verification !!

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
