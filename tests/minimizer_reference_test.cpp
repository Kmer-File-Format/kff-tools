// MinimizerSearcher::get_skmers() compared with the version of kff-tools fc989d2 (copied below,
// with its kmer_minimizer_compute()), on random and low complexity sequences.

#include <random>
#include <string>
#include <vector>

#include "lest.hpp"
#include "sequences.hpp"
#include "skmers.hpp"

using namespace std;

// Defined in sequences.cpp (used by kmer_minimizer_compute())
uint8_t * remove_mini(const uint8_t * seq, uint64_t seq_size, uint64_t start_idx, uint64_t k, uint64_t mini_pos, uint64_t m);

namespace {

int64_t ref_kmer_minimizer_compute(MinimizerSearcher & ms, const uint8_t * seq, uint64_t seq_size, uint64_t start_idx, uint64_t & mini_value, bool & multiple_mini, uint64_t & leftmost_mini)
{
	SlidingWindow window(ms.m, seq, seq_size, start_idx, ms.encoding);
	uint64_t fwd_mini = window.fwd;
	int64_t fwd_idx = 0;
	uint64_t rev_mini = window.rev;
	int64_t rev_idx = 0;

	// Compute all the candidates
	for (int64_t i=1 ; i<=ms.k-ms.m ; i++)
	{
		window.next_char();
		
		// --- Forward testing ---
		// New minimizer
		if (window.fwd < fwd_mini)
		{
			fwd_mini = window.fwd;
			fwd_idx = i;
		}
		// Equivalent value
		else if (window.fwd <= rev_mini and window.fwd == fwd_mini)
		{
			uint64_t window_interleaved = interleaved_size(ms.k, ms.m, i);
			uint64_t fwd_interleaved = interleaved_size(ms.k, ms.m, fwd_idx);

			// Prefer the longest interleaved (the more central minimizer)
			if (window_interleaved > fwd_interleaved)
			{
				fwd_mini = window.fwd;
				fwd_idx = i;
			}
		}

		// --- Reverse testing ---
		// New minimizer
		if (window.rev < rev_mini)
		{
			rev_mini = window.rev;
			rev_idx = i;
		}
		// Equivalent value
		else if (window.rev <= fwd_mini and window.rev == rev_mini)
		{
			uint64_t window_interleaved = interleaved_size(ms.k, ms.m, ms.k-ms.m-i);
			uint64_t rev_interleaved = interleaved_size(ms.k, ms.m, ms.k-ms.m-rev_idx);

			// Prefer the longest interleaved (the more central minimizer)
			if (window_interleaved > rev_interleaved)
			{
				rev_mini = window.rev;
				rev_idx = i;
			}
		}
	}

	// Final selection fwd/rev
	if (fwd_mini < rev_mini) {
		mini_value = fwd_mini;
		multiple_mini = false;
		return fwd_idx + start_idx;
	}
	else if (rev_mini < fwd_mini) {
		mini_value = rev_mini;
		multiple_mini = false;
		return -1 - rev_idx - start_idx;
	}

	// Equality hell
	uint64_t fwd_interleaved = interleaved_size(ms.k, ms.m, fwd_idx);
	uint64_t rev_interleaved = interleaved_size(ms.k, ms.m, ms.k-ms.m-rev_idx);

	multiple_mini = true;
	mini_value = fwd_mini;
	leftmost_mini = min(fwd_idx, rev_idx);

	if (fwd_interleaved > rev_interleaved)
		return fwd_idx + start_idx;
	else if (rev_interleaved > fwd_interleaved)
		return -1 - rev_idx - start_idx;

	// Minimizer + position equality hell
	// remove_mini(uint8_t * seq, uint64_t seq_size, uint64_t start_idx, uint64_t ms.k, uint64_t mini_pos, uint64_t ms.m)
	uint8_t * fwd_seq = remove_mini(seq, seq_size, start_idx, ms.k, fwd_idx, ms.m);
	uint8_t * rev_seq = remove_mini(seq, seq_size, start_idx, ms.k, rev_idx, ms.m);
	ms.rc.rev_comp(rev_seq, ms.k-ms.m);

	// Compare sequences
	int cmp_res = sequence_compare(
						fwd_seq, ms.k-ms.m, 0, ms.k-ms.m-1,
						rev_seq, ms.k-ms.m, 0, ms.k-ms.m-1
	);

	delete[] fwd_seq;
	delete[] rev_seq;

	if (cmp_res <= 0)
		return fwd_idx + start_idx;
	else
		return -1 - rev_idx - start_idx;
}


vector<skmer> ref_get_skmers(MinimizerSearcher & ms, const uint8_t * seq, const uint seq_size) {
	uint64_t unused_mini = 0;
	bool unused_multi = false;
	uint64_t unused_left = 0;
	// Init the sliding window
	SlidingWindow window(ms.m, seq, seq_size, 0, ms.encoding);

	// Init the first minimizer value
	uint64_t minimizer = min(window.fwd, window.rev);
	int64_t mini_pos = window.fwd < window.rev ? 0 : -1;
	bool multi_mini = window.fwd == window.rev;
	uint64_t left_multi_idx = 0;

	// --- Init the first kmer minimizer ---
	for (uint64_t idx=1 ; idx<=ms.k-ms.m ; idx++)
	{
		window.next_char();

		// Get the last minimizer candidate values
		uint64_t current_mini = min(window.fwd, window.rev);
		bool current_rev = window.rev < window.fwd;
		int64_t current_mini_pos = idx;
		bool current_double_mini = window.rev == window.fwd;

		// New minimizer
		if (current_mini < minimizer)
		{
			mini_pos = current_rev ? -1 -current_mini_pos : current_mini_pos;
			minimizer = current_mini;
			multi_mini = current_double_mini;
			left_multi_idx = idx;
		}
		// Multiple minimizer in the same kmer
		else if (current_mini == minimizer)
		{
			multi_mini = true;
			left_multi_idx = mini_pos;
		}
	}

	// Skmer variables init
	vector<skmer> skmers;
	uint64_t skmer_start = 0;

	// Go through all the kmers to split super kmers
	for(int64_t idx=ms.k-ms.m+1 ; idx<=seq_size-ms.m ; idx++)
	{
		window.next_char();

		// Get the last minimizer candidate values
		uint64_t current_mini = window.fwd;
		int64_t current_mini_pos = idx;
		bool current_double_mini = false;
		if (window.rev <= window.fwd)
		{
			current_mini = window.rev;
			current_mini_pos = -1 - idx;

			if (window.rev == window.fwd)
				current_double_mini = true;
		}

		// New minimizer
		if (current_mini < minimizer)
		{
			// Add the previous skmer to the list
			skmers.push_back({skmer_start, static_cast<uint64_t>(idx+ms.m-2), mini_pos, minimizer});
			// Update the values
			skmer_start = idx + ms.m - ms.k;
			minimizer = current_mini;
			mini_pos = current_mini_pos;
			multi_mini = current_double_mini;
			left_multi_idx = idx;
		}
		// Minimizer out of window
		else if ( ((mini_pos >= 0) and (idx >  mini_pos   + ms.k - ms.m))
			 or   ((mini_pos <  0) and (idx > -mini_pos-1 + ms.k - ms.m)) )
		{
			// Save the previous skmer
			skmers.push_back({skmer_start, static_cast<uint64_t>(idx+ms.m-2), mini_pos, minimizer});
			skmer_start = idx + ms.m - ms.k;

			// Compute new minimizer values
			mini_pos = ref_kmer_minimizer_compute(ms, seq, seq_size, idx+ms.m-ms.k, minimizer, multi_mini, left_multi_idx);
		}
		// Same minimizer than previously
		else if (current_mini == minimizer)
		{
			// No change on minimizer value
			multi_mini = true;
			left_multi_idx = current_double_mini ? idx : mini_pos;
			
			current_mini_pos = ref_kmer_minimizer_compute(ms, seq, seq_size, idx+ms.m-ms.k, unused_mini, unused_multi, unused_left);
			// Verify previous skmer continuity
			if (current_mini_pos != mini_pos)
			{
				// Add the previous skmer to the list
				skmers.push_back({skmer_start, static_cast<uint64_t>(idx+ms.m-2), mini_pos, minimizer});
				// Update the values
				skmer_start = idx + ms.m - ms.k;
				mini_pos = current_mini_pos;
			}
		}
		// Cursed case. Recompute everything on each loop (small probability)
		else if (multi_mini)
		{
			// No minimizer value change
			current_mini_pos = ref_kmer_minimizer_compute(ms, seq, seq_size, idx+ms.m-ms.k, unused_mini, unused_multi, unused_left);
			// Verify previous skmer continuity
			if (current_mini_pos != mini_pos)
			{
				// Add the previous skmer to the list
				skmers.push_back({skmer_start, static_cast<uint64_t>(idx+ms.m-2), mini_pos, minimizer});
				// Update the values
				skmer_start = idx + ms.m - ms.k;
				mini_pos = current_mini_pos;
			}
		}
	}

	// Add last skmer
	skmers.push_back({skmer_start, static_cast<uint64_t>(seq_size-1), mini_pos, minimizer});

	return skmers;
}


} // namespace

const lest::test module[] = {

	CASE("MinimizerSearcher::get_skmers gives the superkmers of the original version") {
		mt19937_64 rnd(60);
		uint8_t encodings[2][4] = {{0, 1, 3, 2}, {3, 2, 0, 1}};
		for (int test = 0; test < 3000; test++) {
			uint8_t * encoding = encodings[test % 2];
			const uint m = 1 + rnd() % 31;
			const uint k = m + 1 + rnd() % (test % 5 == 0 ? 40 : 15);
			const size_t n = k + rnd() % 300;
			// Random sequence, low complexity in some tests (many equal minimizer candidates)
			vector<uint8_t> nucl(n);
			const uint alphabet = test % 3 == 0 ? 1 + rnd() % 2 : 4;
			for (uint8_t & x : nucl) x = rnd() % alphabet;
			if (test % 7 == 0)
				for (size_t i = 0; i < n; i++)
					nucl[i] = (i / (1 + test % 5)) % 2;
			// 2-bit encoding, aligned on the right, random padding bits
			const size_t nb_bytes = (n + 3) / 4;
			const size_t offset = (4 - n % 4) % 4;
			vector<uint8_t> seq(nb_bytes + 8);
			for (uint8_t & b : seq) b = rnd();
			for (size_t i = 0; i < n; i++) {
				const size_t pos = offset + i;
				seq[pos / 4] &= ~(3 << (2 * (3 - pos % 4)));
				seq[pos / 4] |= encoding[nucl[i]] << (2 * (3 - pos % 4));
			}

			MinimizerSearcher reference(k, m, encoding);
			MinimizerSearcher searcher(k, m, encoding);
			const vector<skmer> expected = ref_get_skmers(reference, seq.data(), n);
			const vector<skmer> skmers = searcher.get_skmers(seq.data(), n);
			bool same = skmers.size() == expected.size();
			for (size_t i = 0; i < skmers.size() and same; i++)
				same = skmers[i].start_position == expected[i].start_position
				   and skmers[i].stop_position == expected[i].stop_position
				   and skmers[i].minimizer_position == expected[i].minimizer_position
				   and skmers[i].minimizer == expected[i].minimizer;
			EXPECT( same );

			// Minimizer of each kmer
			bool same_kmers = true;
			for (size_t start = 0; start + k <= n and same_kmers; start++) {
				uint64_t mini_ref = 0, left_ref = 0, mini = 0, left = 0;
				bool multi_ref = false, multi = false;
				const int64_t pos_ref = ref_kmer_minimizer_compute(reference, seq.data(), n, start, mini_ref, multi_ref, left_ref);
				const int64_t pos = searcher.kmer_minimizer_compute(seq.data(), n, start, mini, multi, left);
				same_kmers = pos == pos_ref and mini == mini_ref and multi == multi_ref and (not multi or left == left_ref);
			}
			EXPECT( same_kmers );
		}
	},
};

extern lest::tests & specification();

MODULE( specification(), module )
