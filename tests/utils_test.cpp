// Property tests of the sequence utilities against simple reference
// implementations working on strings.

#include <algorithm>
#include <array>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "lest.hpp"
#include "encoding.hpp"
#include "sequences.hpp"

using namespace std;

namespace {

typedef vector<uint8_t> Bytes;

const uint8_t DEFAULT_ENC[4] = {0, 1, 3, 2};

uint8_t code(char c, const uint8_t enc[4]) {
	switch (c) {
		case 'A': return enc[0];
		case 'C': return enc[1];
		case 'G': return enc[2];
		default: return enc[3];
	}
}

/** 2 bits per nucleotide, right aligned, padding bits set to 0 */
Bytes encode(const string & s, const uint8_t enc[4] = DEFAULT_ENC) {
	size_t pad = (4 - s.size() % 4) % 4;
	Bytes b((s.size() + 3) / 4, 0);
	for (size_t i = 0; i < s.size(); i++) {
		size_t p = pad + i;
		b[p / 4] |= code(s[i], enc) << (2 * (3 - p % 4));
	}
	return b;
}

/** Inverse of encode, padding bits ignored */
string decode(const uint8_t * b, size_t n, const uint8_t enc[4] = DEFAULT_ENC) {
	char letters[4];
	letters[enc[0]] = 'A'; letters[enc[1]] = 'C'; letters[enc[2]] = 'G'; letters[enc[3]] = 'T';
	size_t pad = (4 - n % 4) % 4;
	string s(n, '?');
	for (size_t i = 0; i < n; i++) {
		size_t p = pad + i;
		s[i] = letters[(b[p / 4] >> (2 * (3 - p % 4))) & 3];
	}
	return s;
}

/** Value of a sequence (2 bits per nucleotide, first nucleotide = most significant) */
uint64_t value(const string & s, const uint8_t enc[4] = DEFAULT_ENC) {
	uint64_t v = 0;
	for (char c : s)
		v = (v << 2) | code(c, enc);
	return v;
}

string revcomp(const string & s) {
	string r(s.rbegin(), s.rend());
	for (char & c : r)
		c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
	return r;
}

/** Random garbage in the padding bits of the first byte */
void dirty_padding(Bytes & b, size_t n, mt19937 & rnd) {
	size_t pad = (4 - n % 4) % 4;
	if (pad > 0)
		b[0] |= (rnd() & 0xff) & ~((1u << (2 * (4 - pad))) - 1) & 0xff;
}

string random_seq(mt19937 & rnd, size_t n) {
	string s(n, 'A');
	for (char & c : s)
		c = "ACGT"[rnd() % 4];
	return s;
}

vector<array<uint8_t, 4>> some_encodings() {
	vector<array<uint8_t, 4>> encs;
	array<uint8_t, 4> e = {0, 1, 2, 3};
	int i = 0;
	do {
		if (i++ % 5 == 0)
			encs.push_back(e);
	} while (next_permutation(e.begin(), e.end()));
	encs.push_back({0, 1, 3, 2});
	return encs;
}

} // namespace


const lest::test module[] = {

	CASE("Translator converts every nucleotide between two encodings") {
		mt19937 rnd(1);
		for (auto src : some_encodings())
			for (auto dst : some_encodings()) {
				Translator tr(src.data(), dst.data());
				string s = random_seq(rnd, 200);
				Bytes b = encode(s, src.data());
				tr.translate(b.data(), b.size());
				EXPECT( b == encode(s, dst.data()) );
			}
	},

	CASE("RevComp reverse complements sequences of any size and padding") {
		mt19937 rnd(2);
		for (auto enc : some_encodings()) {
			RevComp rc(enc.data());
			for (size_t n = 1; n <= 70; n++) {
				string s = random_seq(rnd, n);
				Bytes b = encode(s, enc.data());
				dirty_padding(b, n, rnd);
				rc.rev_comp(b.data(), n);
				EXPECT( b == encode(revcomp(s), enc.data()) );
			}
		}
	},

	CASE("RevComp reverses the data of the kmers") {
		RevComp rc(DEFAULT_ENC);
		for (uint64_t data_size = 1; data_size <= 3; data_size++)
			for (uint64_t nb_kmers = 1; nb_kmers <= 7; nb_kmers++) {
				Bytes data(data_size * nb_kmers), expected(data.size());
				for (size_t i = 0; i < data.size(); i++)
					data[i] = (uint8_t)i;
				for (uint64_t kmer = 0; kmer < nb_kmers; kmer++)
					for (uint64_t d = 0; d < data_size; d++)
						expected[(nb_kmers - 1 - kmer) * data_size + d] = data[kmer * data_size + d];
				rc.rev_data(data.data(), data_size, nb_kmers);
				EXPECT( data == expected );
			}
	},

	CASE("Stringifyer converts sequences of any size, padding bits ignored") {
		mt19937 rnd(3);
		for (auto enc : some_encodings()) {
			Stringifyer strif(enc.data());
			for (size_t n = 1; n <= 70; n++) {
				string s = random_seq(rnd, n);
				Bytes b = encode(s, enc.data());
				dirty_padding(b, n, rnd);
				EXPECT( strif.translate(b.data(), n) == s );
				if (n <= 32)
					EXPECT( strif.translate(value(s, enc.data()), n) == s );
			}
		}
	},

	CASE("Stringifyer::write gives the same characters as translate") {
		mt19937 rnd(12);
		for (auto enc : some_encodings()) {
			Stringifyer strif(enc.data());
			for (size_t n = 1; n <= 70; n++) {
				string s = random_seq(rnd, n);
				Bytes b = encode(s, enc.data());
				dirty_padding(b, n, rnd);
				string out(n + 8, '#');
				char * end = strif.write(b.data(), n, &out[0]);
				EXPECT( end == &out[0] + n );
				EXPECT( out.substr(0, n) == strif.translate(b.data(), n) );
				EXPECT( out.substr(n) == string(8, '#') );
			}
		}
	},

	CASE("Binarizer handles lower case letters and any size") {
		mt19937 rnd(4);
		for (auto enc : some_encodings()) {
			Binarizer bz(enc.data());
			for (size_t n = 1; n <= 70; n++) {
				string s = random_seq(rnd, n);
				string mixed = s;
				for (char & c : mixed)
					if (rnd() % 2)
						c = tolower(c);
				Bytes b((n + 3) / 4, 0xff);
				bz.translate(mixed, n, b.data());
				EXPECT( b == encode(s, enc.data()) );
			}
		}
	},

	CASE("Binarizer encodes unknown characters with the value 0") {
		uint8_t enc[4] = {2, 0, 3, 1};  // C has the value 0
		Binarizer bz(enc);
		Bytes b(2, 0xff);
		bz.translate(string("NAxGT"), 5, b.data());
		EXPECT( b == encode("CACGT", enc) );
	},

	CASE("subsequence extracts any slice of a sequence") {
		mt19937 rnd(5);
		for (size_t n = 1; n <= 40; n++) {
			string s = random_seq(rnd, n);
			Bytes b = encode(s);
			for (size_t begin = 0; begin < n; begin++)
				for (size_t end = begin; end < n; end++) {
					size_t len = end - begin + 1;
					Bytes out((len + 3) / 4 + 2, 0);
					subsequence(b.data(), n, out.data(), begin, end);
					EXPECT( decode(out.data(), len) == s.substr(begin, len) );
				}
		}
	},

	CASE("sequence_compare follows the encoding order") {
		const size_t max_len = 40;
		mt19937 rnd(6);
		for (int test = 0; test < 3000; test++) {
			size_t n1 = 1 + rnd() % 40, n2 = 1 + rnd() % 40;
			string s1 = random_seq(rnd, n1), s2 = random_seq(rnd, n2);
			size_t len = 1 + rnd() % min(max_len, min(n1, n2));
			size_t b1 = rnd() % (n1 - len + 1), b2 = rnd() % (n2 - len + 1);
			// Equal slices half of the time
			if (test % 2 == 0)
				s2.replace(b2, len, s1.substr(b1, len));
			Bytes e1 = encode(s1), e2 = encode(s2);
			int result = sequence_compare(e1.data(), n1, b1, b1 + len - 1, e2.data(), n2, b2, b2 + len - 1);
			string v1 = s1.substr(b1, len), v2 = s2.substr(b2, len);
			for (char & c : v1) c = '0' + code(c, DEFAULT_ENC);
			for (char & c : v2) c = '0' + code(c, DEFAULT_ENC);
			int expected = v1 < v2 ? -1 : (v1 == v2 ? 0 : 1);
			EXPECT( result == expected );
		}
		// Slices of different sizes are ordered by size
		Bytes e = encode("ACGTACGT");
		EXPECT( sequence_compare(e.data(), 8, 0, 2, e.data(), 8, 0, 3) == -1 );
		EXPECT( sequence_compare(e.data(), 8, 0, 3, e.data(), 8, 1, 2) == 1 );
	},

	CASE("seq_to_uint gives the value of the last 32 nucleotides") {
		mt19937 rnd(7);
		for (size_t n = 1; n <= 70; n++) {
			string s = random_seq(rnd, n);
			Bytes b = encode(s);
			EXPECT( seq_to_uint(b.data(), n) == value(s.substr(n > 32 ? n - 32 : 0)) );
		}
	},

	CASE("subseq_to_uint gives the value of a slice") {
		mt19937 rnd(8);
		for (size_t n = 1; n <= 45; n++) {
			string s = random_seq(rnd, n);
			Bytes b = encode(s);
			for (size_t begin = 0; begin < n; begin++)
				for (size_t end = begin; end < n; end++) {
					size_t first = end - begin + 1 > 32 ? end - 31 : begin;
					EXPECT( subseq_to_uint(b.data(), n, begin, end) == value(s.substr(first, end - first + 1)) );
				}
		}
	},

	CASE("uint_to_seq writes the bytes of a sequence value") {
		mt19937 rnd(9);
		for (size_t n = 1; n <= 16; n++) {
			string s = random_seq(rnd, n);
			Bytes b((n + 3) / 4, 0xff);
			uint_to_seq((uint)value(s), b.data(), n);
			EXPECT( b == encode(s) );
		}
	},

	CASE("SlidingWindow slides forward and reverse complement windows") {
		mt19937 rnd(10);
		for (auto enc : some_encodings())
			for (size_t w = 1; w <= 32; w += 3) {
				string s = random_seq(rnd, w + 30);
				Bytes b = encode(s, enc.data());
				SlidingWindow window(w, b.data(), s.size(), 0, enc.data());
				for (size_t pos = 0; pos + w <= s.size(); pos++) {
					if (pos > 0)
						window.next_char();
					EXPECT( window.fwd == value(s.substr(pos, w), enc.data()) );
					EXPECT( window.rev == value(revcomp(s.substr(pos, w)), enc.data()) );
				}
			}
	},

	CASE("Super-kmers are contiguous and all their kmers contain the minimizer") {
		mt19937 rnd(11);
		const uint ks[] = {5, 8, 11, 15, 21, 31};
		for (uint k : ks)
			for (uint m = 2; m < k && m <= 12; m += 3)
				for (int test = 0; test < 20; test++) {
					size_t n = k + rnd() % 150;
					string s = random_seq(rnd, n);
					// Low complexity sequences generate many equal minimizer candidates
					if (test % 4 == 0)
						for (char & c : s)
							c = "AC"[rnd() % 2];
					Bytes b = encode(s);
					MinimizerSearcher searcher(k, m, (uint8_t *)DEFAULT_ENC);
					vector<skmer> skmers = searcher.get_skmers(b.data(), n);

					EXPECT( skmers.front().start_position == 0u );
					EXPECT( skmers.back().stop_position == n - 1 );
					for (size_t i = 0; i < skmers.size(); i++) {
						const skmer & sk = skmers[i];
						EXPECT( sk.stop_position - sk.start_position + 1 >= k );
						if (i > 0)  // consecutive super-kmers overlap by k-1 nucleotides
							EXPECT( sk.start_position == skmers[i - 1].stop_position + 2 - k );

						// The minimizer occurrence
						bool reverse = sk.minimizer_position < 0;
						uint64_t mini_pos = reverse ? -sk.minimizer_position - 1 : sk.minimizer_position;
						string mmer = s.substr(mini_pos, m);
						EXPECT( value(reverse ? revcomp(mmer) : mmer) == sk.minimizer );

						// Every kmer of the super-kmer contains it and has it as minimal m-mer
						for (uint64_t kmer = sk.start_position; kmer + k <= sk.stop_position + 1; kmer++) {
							EXPECT( kmer <= mini_pos );
							EXPECT( mini_pos + m <= kmer + k );
							uint64_t min_value = UINT64_MAX;
							for (uint64_t p = kmer; p + m <= kmer + k; p++) {
								string candidate = s.substr(p, m);
								min_value = min(min_value, min(value(candidate), value(revcomp(candidate))));
							}
							EXPECT( min_value == sk.minimizer );
						}
					}
				}
	},
};

extern lest::tests & specification();

MODULE( specification(), module )
