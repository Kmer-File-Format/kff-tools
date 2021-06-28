#include "encoding.hpp"
#include "sequences.hpp"

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


RevComp::RevComp(const uint8_t encoding[4]) {
	this->reverse[encoding[0]] = encoding[3];
	this->reverse[encoding[1]] = encoding[2];
	this->reverse[encoding[2]] = encoding[1];
	this->reverse[encoding[3]] = encoding[0];

	for (uint i=0 ; i<256 ; i++) {
		uint8_t val = i;
		uint8_t rc_val = 0;

		for (uint j=0 ; j<4 ; j++) {
			rc_val <<= 2;
			rc_val += this->reverse[val & 0b11];
			val >>= 2;
		}
		this->translations[i] = rc_val;
	}
}

void RevComp::rev_comp(uint8_t * seq, const uint64_t seq_size) const {
	uint nb_bytes = (seq_size + 3) / 4;
	// reverse and translate each byte
	for (uint byte_idx=0 ; byte_idx<(nb_bytes+1)/2 ; byte_idx++) {
		uint8_t save = seq[byte_idx];
		seq[byte_idx] = seq[nb_bytes-1-byte_idx];
		seq[nb_bytes-1-byte_idx] = save;
	}

	uint8_t offset = (4 - (seq_size % 4)) % 4;
	rightshift8(seq, nb_bytes, offset*2);
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


string Stringifyer::translate(const uint8_t * sequence, const size_t nucl_length) const {
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



Binarizer::Binarizer(const uint8_t encoding[4]) {
	string l = "A";
	this->lookup[l] = 0 | (encoding[0] & 0b11);
	l = "C";
	this->lookup[l] = 0 | (encoding[1] & 0b11);
	l = "G";
	this->lookup[l] = 0 | (encoding[2] & 0b11);
	l = "T";
	this->lookup[l] = 0 | (encoding[3] & 0b11);

	const string letters[] = {"A", "C", "G", "T"};
	// First letter
	for (string letter1 : letters) {
		// Second letter
		for (string letter2 : letters) {
			this->lookup[letter1 + letter2] = (lookup[letter1] << 2) | lookup[letter2];
			// Third letter
			for (string letter3 : letters) {
				this->lookup[letter1 + letter2 + letter3] = (lookup[letter1] << 4) | (lookup[letter2] << 2) | lookup[letter3];
				// Last letter
				for (string letter4: letters) {
					this->lookup[letter1 + letter2 + letter3 + letter4] = (lookup[letter1] << 6) | (lookup[letter2] << 4) | (lookup[letter3] << 2) | lookup[letter4];
				}
			}
		}
	}
}


void Binarizer::translate(std::string sequence, uint8_t * binarized) {
	uint k = sequence.length();
	uint truncated = k % 4;
	uint remaining_bytes = k / 4;

	uint off_byte = 0;
	if (truncated > 0) {
		string prefix = sequence.substr(0, truncated);
		binarized[0] = this->lookup[prefix];
		off_byte = 1;
	}

	for (uint i=0 ; i<remaining_bytes ; i++) {
		string slice = sequence.substr(truncated + i * 4, 4);
		binarized[i + off_byte] = this->lookup[slice];
	}
}
