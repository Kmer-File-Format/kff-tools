#include <cstdint>
#include <iostream>


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
