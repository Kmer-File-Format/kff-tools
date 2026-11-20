// The interleaved keys used to sort the kmer matrix must give the order of
// Compact::interleaved_compare_kmers.

#include <algorithm>
#include <cstring>
#include <random>
#include <vector>

#include "lest.hpp"
#include "compact.hpp"
#include "sequences.hpp"

using namespace std;

const lest::test module[] = {

	CASE("Compact interleaved keys follow interleaved_compare_kmers") {
		mt19937_64 rnd(50);
		for (int config = 0; config < 400; config++) {
			const uint nb_nucl = 1 + rnd() % 64;
			Compact comp;
			comp.m = 1 + rnd() % 15;
			comp.k = comp.m + nb_nucl;
			comp.data_size = rnd() % 3;
			comp.mini_pos_size = 1 + rnd() % 2;
			comp.bytes_compacted = (nb_nucl + 3) / 4;
			const uint mini_pos = rnd() % (nb_nucl + 1);

			// Random kmers (random padding bits), some of them equal or almost equal
			vector<long> positions;
			vector<uint8_t> seq(comp.bytes_compacted), data(comp.data_size + 1);
			const int nb_kmers = 2 + rnd() % 200;
			for (int i = 0; i < nb_kmers; i++) {
				if (i == 0 or rnd() % 3) {
					for (uint8_t & b : seq) b = rnd();
				} else {
					// change a few nucleotides of the previous kmer (or none)
					for (int n = rnd() % 3; n > 0; n--)
						seq[rnd() % seq.size()] ^= 3 << (2 * (rnd() % 4));
					seq[0] ^= rnd() & 0xF0 & ~((1u << (2 * (4 - (4 - nb_nucl % 4) % 4))) - 1);
				}
				for (uint8_t & b : data) b = rnd();
				positions.push_back(comp.add_kmer_to_buffer(seq.data(), data.data(), mini_pos));
			}
			vector<uint8_t *> kmers;
			for (long pos : positions)
				kmers.push_back(comp.kmer_buffer + pos);

			// Same comparisons
			bool same_order = true;
			for (int test = 0; test < 400; test++) {
				uint8_t * a = kmers[rnd() % kmers.size()];
				uint8_t * b = kmers[rnd() % kmers.size()];
				uint64_t ah, al, bh, bl;
				comp.interleaved_key(a, ah, al);
				comp.interleaved_key(b, bh, bl);
				int expected = comp.interleaved_compare_kmers(a, b);
				int key = (ah != bh) ? (ah < bh ? -1 : 1) : (al != bl ? (al < bl ? -1 : 1) : 0);
				same_order = same_order and ((expected < 0) == (key < 0)) and ((expected == 0) == (key == 0));
			}
			EXPECT( same_order );

			// Same sorted column as std::sort with the original comparison: the same records at
			// each position (identical records can be swapped, the output is the same)
			vector<vector<uint8_t *> > matrix(1, kmers);
			comp.sort_matrix(matrix);
			vector<uint8_t *> expected_column = kmers;
			sort(expected_column.begin(), expected_column.end(), [&comp](const uint8_t * a, const uint8_t * b) {
				return comp.interleaved_compare_kmers(a, b) < 0;
			});
			const size_t record_size = comp.bytes_compacted + comp.data_size + comp.mini_pos_size;
			bool same_column = matrix[0].size() == expected_column.size();
			for (size_t j = 0; j < expected_column.size() and same_column; j++)
				same_column = memcmp(matrix[0][j], expected_column[j], record_size) == 0;
			EXPECT( same_column );
		}
	},

	CASE("extract_kmer gives the bytes of subsequence") {
		mt19937_64 rnd(61);
		bool same = true;
		uint64_t nb_tests = 0;
		for (uint64_t nb_nucl = 1; nb_nucl <= 28; nb_nucl++)
			for (uint64_t seq_size = nb_nucl; seq_size <= nb_nucl + 40; seq_size++)
				for (int rep = 0; rep < 3; rep++) {
					// Random bytes (padding bits included), 1 more byte for subsequence()
					vector<uint8_t> seq((seq_size + 3) / 4 + 1);
					for (uint8_t & b : seq) b = rnd();
					for (uint64_t begin = 0; begin + nb_nucl <= seq_size; begin++) {
						vector<uint8_t> expected((nb_nucl + 3) / 4 + 2, 0), result((nb_nucl + 3) / 4 + 2, 0);
						subsequence(seq.data(), seq_size, expected.data(), begin, begin + nb_nucl - 1);
						extract_kmer(seq.data(), seq_size, result.data(), begin, nb_nucl);
						same = same and equal(result.begin(), result.begin() + (nb_nucl + 3) / 4, expected.begin());
						nb_tests += 1;
					}
				}
		EXPECT( same );
		EXPECT( nb_tests > 30000u );
	},
};

extern lest::tests & specification();

MODULE( specification(), module )
