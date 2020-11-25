#include <cstdint>
#include <iostream>
#include <string>


/** 
	* Object used to translate compacted sequences from an encoding to another.
	* A translation Byte size lookup table is computed during the object creation
	* and used to fast translate sequences Byte per Byte.
	**/
class Translator {
private:
	uint8_t lookup[256];

public:
	/**
		* Constructor that construct a 256 Bytes lookup table for fast translation
		*
		* @param source The source encoding. A 4 cell array where each of the numbers
		* from 0 to 3 must be present.
		* @param destination The destination encoding. A 4 cell array where each of
		* the numbers from 0 to 3 must be present.
		**/
	Translator(uint8_t source[4], uint8_t destination[4]);
	/**
	  * Inplace translate the sequence from the source encoding to the destination.
	  *
	  * @param sequence 2-bit compacted sequence that will be translated regarding
	  * the encodings.
	  * @param byte_length Length in Bytes of the sequence.
	  **/
	void translate(uint8_t * sequence, size_t byte_length);
};


/** 
	* Object used to convert compacted sequences from 2 bits per nucleotide to string.
	* A translation Byte size lookup table is computed during the object creation
	* and used to fast transform sequences Byte per Byte.
	**/
class Stringifyer {
private:
	std::string lookup[256];

public:
	/**
		* Constructor that construct a 256 Bytes lookup table for fast conversion
		*
		* @param source The sequence encoding
		**/
	Stringifyer(uint8_t encoding[4]);
	/**
	  * read a 2-bits/nucl sequence and return a coresponding string
	  *
	  * @param sequence 2-bit compacted sequence that will be converted regarding
	  * the encodings.
	  * @param nucl_length Length in nucleotides of the sequence.
	  **/
	std::string translate(uint8_t * sequence, size_t nucl_length);
};
