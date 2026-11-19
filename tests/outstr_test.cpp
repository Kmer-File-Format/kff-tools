// outstr data formatting compared with the original implementation.

#include <random>
#include <string>

#include "lest.hpp"
#include "outstr.hpp"

using namespace std;

namespace {

// format_data() of kff-tools 9b90270 (outstr.cpp)
string reference_format_data(uint8_t * data, size_t data_size) {
	if (data_size == 0)
		return "";
	else if (data_size < 8) {
		uint val = data[0];
		for (uint i=1 ; i<data_size ; i++) {
			val <<= 8;
			val += data[i];
		}
		return std::to_string(val);
	} else {
		string val = "";
		for (uint i=0 ; i<data_size ; i++) {
			val += "[" + std::to_string((uint)data[i]) + "]";
		}
		return val;
	}
}

} // namespace

const lest::test module[] = {

	CASE("outstr writes the data of the kmers as before") {
		mt19937 rnd(20);
		for (size_t data_size = 0; data_size <= 12; data_size++)
			for (int test = 0; test < 2000; test++) {
				vector<uint8_t> data(data_size + 1);
				for (uint8_t & byte : data)
					// Extreme values half of the time
					byte = test % 2 ? rnd() % 256 : (rnd() % 2 ? 0 : 255);
				string out(outstr_data_chars(data_size) + 4, '#');
				char * end = outstr_write_data(data.data(), data_size, &out[0]);
				string expected = reference_format_data(data.data(), data_size);
				EXPECT( string(&out[0], end) == expected );
				EXPECT( expected.size() <= outstr_data_chars(data_size) );
			}
	},
};

extern lest::tests & specification();

MODULE( specification(), module )
