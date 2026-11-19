// The in-memory block encoder must give exactly the bytes written by the API.

#include <cstdio>
#include <fstream>
#include <iterator>
#include <random>
#include <string>
#include <vector>

#include "lest.hpp"
#include "blocks.hpp"
#include "kff_io.hpp"

using namespace std;

const lest::test module[] = {

	CASE("MinimizerBlockEncoder writes the bytes of Section_Minimizer::write_compacted_sequence") {
		mt19937_64 rnd(30);
		const string filename = "blocks_test_" + to_string(rnd()) + ".kff";
		for (int config = 0; config < 300; config++) {
			const uint64_t m = 1 + rnd() % 20;
			const uint64_t k = m + rnd() % (config % 3 == 0 ? 60 : 15);
			const uint64_t max = 1 + rnd() % (config % 5 == 0 ? 300 : 40);
			const uint64_t data_size = rnd() % 3;
			MinimizerBlockEncoder encoder(k, m, max, data_size);

			vector<uint8_t> expected_blocks;
			long section_start = 0;
			{
				Kff_file file(filename, "w");
				file.set_indexation(false);
				Section_GV sgv(&file);
				sgv.write_var("k", k);
				sgv.write_var("m", m);
				sgv.write_var("max", max);
				sgv.write_var("data_size", data_size);
				sgv.close();
				Section_Minimizer sm(&file);
				section_start = sm.beginning;
				for (int b = 0; b < 20; b++) {
					const uint64_t nb_kmers = 1 + rnd() % max;
					const uint64_t seq_size = nb_kmers + k - 1;
					const uint64_t mini_pos = rnd() % (seq_size - m + 1);
					vector<uint8_t> seq((seq_size + 3) / 4), data(nb_kmers * data_size + 1);
					for (uint8_t & byte : seq) byte = rnd();  // random padding bits included
					for (uint8_t & byte : data) byte = rnd();
					encoder.append(expected_blocks, seq.data(), seq_size, mini_pos, data.data());
					sm.write_compacted_sequence(seq.data(), seq_size, mini_pos, data.data());
				}
				sm.close();
				file.close();
			}
			ifstream is(filename, ios::binary);
			vector<uint8_t> bytes((istreambuf_iterator<char>(is)), istreambuf_iterator<char>());
			// Section: 'm' + minimizer + number of blocks (8 bytes), then the blocks, then "KFF"
			const size_t blocks_start = section_start + 1 + encoder.minimizer_bytes + 8;
			EXPECT( vector<uint8_t>(bytes.begin() + blocks_start, bytes.end() - 3) == expected_blocks );
		}
		remove(filename.c_str());
	},

	CASE("MinimizerBlockEncoder writes the bytes of write_compacted_sequence_without_mini") {
		mt19937_64 rnd(31);
		const string filename = "blocks_test_" + to_string(rnd()) + ".kff";
		for (int config = 0; config < 300; config++) {
			const uint64_t m = 1 + rnd() % 20;
			const uint64_t k = m + rnd() % (config % 3 == 0 ? 60 : 15);
			const uint64_t max = 1 + rnd() % (config % 5 == 0 ? 300 : 40);
			const uint64_t data_size = rnd() % 3;
			MinimizerBlockEncoder encoder(k, m, max, data_size);

			vector<uint8_t> expected_blocks;
			long section_start = 0;
			{
				Kff_file file(filename, "w");
				file.set_indexation(false);
				Section_GV sgv(&file);
				sgv.write_var("k", k);
				sgv.write_var("m", m);
				sgv.write_var("max", max);
				sgv.write_var("data_size", data_size);
				sgv.close();
				Section_Minimizer sm(&file);
				section_start = sm.beginning;
				for (int b = 0; b < 20; b++) {
					const uint64_t nb_kmers = 1 + rnd() % max;
					const uint64_t no_mini_size = nb_kmers + k - 1 - m;
					const uint64_t mini_pos = rnd() % (no_mini_size + 1);
					vector<uint8_t> seq((no_mini_size + 3) / 4 + 1), data(nb_kmers * data_size + 1);
					for (uint8_t & byte : seq) byte = rnd();
					for (uint8_t & byte : data) byte = rnd();
					encoder.append_without_mini(expected_blocks, seq.data(), no_mini_size, mini_pos, data.data());
					sm.write_compacted_sequence_without_mini(seq.data(), no_mini_size, mini_pos, data.data());
				}
				sm.close();
				file.close();
			}
			ifstream is(filename, ios::binary);
			vector<uint8_t> bytes((istreambuf_iterator<char>(is)), istreambuf_iterator<char>());
			const size_t blocks_start = section_start + 1 + encoder.minimizer_bytes + 8;
			EXPECT( vector<uint8_t>(bytes.begin() + blocks_start, bytes.end() - 3) == expected_blocks );
		}
		remove(filename.c_str());
	},
};

extern lest::tests & specification();

MODULE( specification(), module )
