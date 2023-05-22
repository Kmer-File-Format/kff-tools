#include <cstring>
#include <cassert>
#include <algorithm>

#include "sequences.hpp"
#include "skmers.hpp"


using namespace std;


/* Bitshift to the left all the bits in the array with a maximum of 7 bits.
 * Overflow on the left will be set into the previous cell.
 */
void leftshift8(uint8_t * bitarray, size_t length, size_t bitshift) {
	assert(bitshift < 8);

	for (uint64_t i=0 ; i<length-1 ; i++) {
		bitarray[i] = (bitarray[i] << bitshift) | (bitarray[i+1] >> (8-bitshift));
	}
	bitarray[length-1] <<= bitshift;
}

/* Similar to the previous function but on the right */
void rightshift8(uint8_t * bitarray, size_t length, size_t bitshift) {
	assert(bitshift < 8);

	for (uint64_t i=length-1 ; i>0 ; i--) {
		bitarray[i] = (bitarray[i-1] << (8-bitshift)) | (bitarray[i] >> bitshift);
	}
	bitarray[0] >>= bitshift;
}

/* Fusion to bytes into one.
 * The merge_index higher bits are from left_bits the others from right_bits
 */
uint8_t fusion8(uint8_t left_bits, uint8_t right_bits, size_t merge_index) {
	uint8_t mask = 0xFF << (8-merge_index);
	return (left_bits & mask) | (right_bits & ~mask);
}

int KffSeqStream::next_sequence(uint8_t * & seq, uint max_seq_size, uint8_t * & data, uint max_data_size) {
// int KffSeqStream::next_sequence(uint8_t * & seq, uint max_seq_size, uint8_t * & data, uint max_data_size) {
	if (this->reader.has_next()) {
		uint max_seq = this->reader.k + this->reader.max - 1;
		uint max_data = this->reader.max * this->reader.data_size;
		if (max_seq_size < max_seq or max_data_size < max_data)
			return -1;
		else
			return this->reader.next_block(seq, data);
	}

	return 0;
}




// #include <iostream>
void subsequence(const uint8_t * sequence, const uint seq_size, uint8_t * extracted, const uint begin_nucl, const uint end_nucl) {
	// Extract the correct slice
	uint seq_left_offset = (4 - seq_size % 4) % 4;
	uint extract_start_byte = (seq_left_offset + begin_nucl) / 4;
	uint extract_stop_byte = (seq_left_offset + end_nucl) / 4;

	memcpy(extracted, sequence + extract_start_byte, extract_stop_byte - extract_start_byte + 1);

	// Align the bits
	uint extract_left_offset = (seq_left_offset + begin_nucl) % 4;
	uint extract_right_offset = (seq_size - end_nucl - 1) % 4;

	if (extract_right_offset < 4 - extract_left_offset) {
		rightshift8(extracted, extract_stop_byte - extract_start_byte + 1, extract_right_offset * 2);
	} else {
		leftshift8(extracted, extract_stop_byte - extract_start_byte + 1, (4 - extract_right_offset) * 2);
	}
}


void subsequence_bis(const uint8_t * sequence, const uint seq_size, uint8_t * extracted, const uint begin_nucl, const uint end_nucl) {
	// Extract the correct slice
	uint seq_left_offset = (4 - seq_size % 4) % 4;
	uint extract_start_byte = (seq_left_offset + begin_nucl) / 4;

    uint64_t reads = end_nucl - begin_nucl + 1; // number of nucls to read
    uint64_t position = (seq_left_offset + begin_nucl) % 4; // the position of the nucls to read
    uint64_t mask = 0b11 << ((3 - position) * 2); // the mask to read the wanted bits
    uint64_t section = 0;

    for (uint i = 0; i < reads; i++) {
        extracted[section] <<= 2;
        extracted[section] |= ((sequence[extract_start_byte]) & mask) >> ((3 - position) * 2);
        if (position == 3) { // we arrive at the end of the section of the vector
            extract_start_byte++;
            position = 0;
            mask = 0b11000000;
        } else { // we continue to read in the same section of the vector
            mask >>= 2;
            position++;
        }
        if (i % 4 == 0 && i > 0) { // end of the section of extracted
            section++;
        }
    }
}

int sequence_compare(const uint8_t * seq1, const uint seq1_size,
											const uint seq1_start, const uint seq1_stop,
											const uint8_t * seq2, const uint seq2_size,
											const uint seq2_start, const uint seq2_stop) {

	// If inequal size
	if (seq1_stop - seq1_start != seq2_stop - seq2_start)
		return (seq1_stop - seq1_start) < (seq2_stop - seq2_start) ? -1 : 1;

	// Extraction of subsequences
	uint subseq_size = seq1_stop - seq1_start + 1;
	uint subseq_bytes = 1 + (subseq_size + 3) / 4;


    uint8_t * subseq1 = new uint8_t[subseq_bytes];
	memset(subseq1, 0, subseq_bytes);
	subsequence_bis(seq1, seq1_size, subseq1, seq1_start, seq1_stop);
	uint8_t * subseq2 = new uint8_t[subseq_bytes];
	memset(subseq2, 0, subseq_bytes);
	subsequence_bis(seq2, seq2_size, subseq2, seq2_start, seq2_stop);

	// comparison of the subsequences (same size)
	uint return_val = 0;
	// First byte
	uint offset = (4 - (subseq_size % 4)) % 4;
	uint mask = (1 << (2 * (4 - offset))) - 1;
	if ((subseq1[0] & mask) != (subseq2[0] & mask))
		return_val = (subseq1[0] & mask) < (subseq2[0] & mask) ? -1 : 1;


	// All following bytes
	for (uint b=1 ; b<subseq_bytes and return_val==0 ; b++) {
		if (subseq1[b] != subseq2[b]) {
			return_val = subseq1[b] < subseq2[b] ? -1 : 1;
		}
	}

	delete[] subseq1;
	delete[] subseq2;
	return return_val;
}


SlidingWindow::SlidingWindow(uint64_t window_size, const uint8_t * seq, uint64_t seq_size, uint64_t start_idx, uint8_t encoding[4])
							: fwd(0), rev(0), seq(seq), window_size(window_size)
{
	this->seq_offset = (4 - (seq_size % 4)) % 4;
	this->idx = this->seq_offset + start_idx;

	// Init a mask of the window size
	if (window_size == 32)
		this->mask = 0xFFFFFFFFFFFFFFFFul;
	else
		this->mask = (1ul << (2 * window_size)) - 1;

	// Prepare the complement table
	this->complement[encoding[0]] = encoding[3];
	this->complement[encoding[1]] = encoding[2];
	this->complement[encoding[2]] = encoding[1];
	this->complement[encoding[3]] = encoding[0];

	// Init the first w-mer
	for (uint64_t i=0 ; i<window_size ; i++)
		this->next_char();
};


/** Compute the next forward and reverse sequence. You have to be sure that there are remaining nucleotides in the sequence.
**/
void SlidingWindow::next_char()
{
	// Extract nucleotide
	uint64_t byte_idx = this->idx / 4;
	uint64_t bit_shift = 2 * (3 - (this->idx % 4));
	uint64_t nucl = (this->seq[byte_idx] >> bit_shift) & 0b11;

	// Update candidate
	this->fwd <<= 2;
	this->fwd  += nucl;
	this->fwd  &= this->mask;

	nucl = this->complement[nucl];
	this->rev >>= 2;
	this->rev += nucl << (2 * (window_size - 1));

	this->idx += 1;
};


/** Copy and return the input kmer without its minimizer.
 * @param The sequence containing the kmer.
 * @param seq_size Sequence size. Allows to know the number of useless bits at the begining.
 * @param start_idx Start (in nucleotides) of the kmer inside of the sequence.
 * @param k kmer size
 * @param mini_pos Minimizer position inside of the sequence
 * @param m Minimizer size
 * @return A memory allocated sequence of size k-m
 **/
uint8_t * remove_mini(const uint8_t * seq, uint64_t seq_size, uint64_t start_idx, uint64_t k, uint64_t mini_pos, uint64_t m)
{
	uint64_t seq_offset = (4 - (seq_size % 4)) % 4;
	uint64_t subseq_offset = (4 - ((k - m) % 4)) % 4;

	// init the sub sequence
	uint8_t * subseq = new uint8_t[(k - m + 3) / 4];
	memset(subseq, 0, (k - m + 3) / 4);

	// Prefix
	for (uint64_t nucl_idx=0 ; nucl_idx<mini_pos ; nucl_idx++)
	{
		uint64_t seq_idx = start_idx + nucl_idx;
		uint64_t abs_seq_idx = seq_idx + seq_offset;

		// Extract the nucleotide
		uint64_t seq_byte = abs_seq_idx / 4;
		uint64_t seq_bit = 2 * (3 - (abs_seq_idx % 4));
		uint64_t nucl = (seq[seq_byte] >> seq_bit) & 0b11;

		uint64_t subseq_idx = nucl_idx;
		uint64_t abs_subseq_idx = subseq_idx + subseq_offset;

		// Insert the nucleotide
		uint64_t subseq_byte = abs_subseq_idx / 4;
		uint64_t subseq_bit = 2 * (3 - (abs_subseq_idx % 4));
		subseq[subseq_byte] |= nucl << subseq_bit;
	}

	// Suffix
	for (uint64_t nucl_idx=mini_pos+m ; nucl_idx<k ; nucl_idx++)
	{
		uint64_t seq_idx = start_idx + nucl_idx;
		uint64_t abs_seq_idx = seq_idx + seq_offset;

		// Extract the nucleotide
		uint64_t seq_byte = abs_seq_idx / 4;
		uint64_t seq_bit = 2 * (3 - (abs_seq_idx % 4));
		uint64_t nucl = (seq[seq_byte] >> seq_bit) & 0b11;

		uint64_t subseq_idx = nucl_idx-m;
		uint64_t abs_subseq_idx = subseq_idx + subseq_offset;

		// Insert the nucleotide
		uint64_t subseq_byte = abs_subseq_idx / 4;
		uint64_t subseq_bit = 2 * (3 - (abs_subseq_idx % 4));
		subseq[subseq_byte] |= nucl << subseq_bit;
	}

	return subseq;
};


/** Compute the position of the minimizer inside of the kmer starting at position start_idx
 * @param seq A 2-bits encoded sequence of nucleotides
 * @param seq_size Size of the sequence in nucleotides. If not multiple of 4, the first seq byte
 * contains non used bits.
 * @param start_idx First nucleotide to look at inside the sequence. Non used bits are 
 * automatically skipped.
 * @return index of the minimizer in the sequence. If the minimizer is on the reverse side, the 
 * return value is (- 1 - index)
 **/
int64_t MinimizerSearcher::kmer_minimizer_compute(const uint8_t * seq, uint64_t seq_size, uint64_t start_idx)
{
	uint64_t mini = 0;
	bool multi_mini = false;
	uint64_t leftmost_mini = 0;

	return this->kmer_minimizer_compute(seq, seq_size, start_idx, mini, multi_mini, leftmost_mini);
}

/** Compute the position of the minimizer inside of the kmer starting at position start_idx
 * @param seq A 2-bits encoded sequence of nucleotides
 * @param seq_size Size of the sequence in nucleotides. If not multiple of 4, the first seq byte
 * contains non used bits.
 * @param start_idx First nucleotide to look at inside the sequence. Non used bits are 
 * automatically skipped.
 * @param mini_value Minimizer value
 * @param multiple_mini True if the minimizer value is present at multiple spots
 * @return index of the minimizer in the sequence. If the minimizer is on the reverse side, the 
 * return value is (- 1 - index)
 **/
int64_t MinimizerSearcher::kmer_minimizer_compute(const uint8_t * seq, uint64_t seq_size, uint64_t start_idx, uint64_t & mini_value, bool & multiple_mini, uint64_t & leftmost_mini)
{
	SlidingWindow window(m, seq, seq_size, start_idx, this->encoding);
	uint64_t fwd_mini = window.fwd;
	int64_t fwd_idx = 0;
	uint64_t rev_mini = window.rev;
	int64_t rev_idx = 0;

	// Compute all the candidates
	for (int64_t i=1 ; i<=k-m ; i++)
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
			uint64_t window_interleaved = interleaved_size(k, m, i);
			uint64_t fwd_interleaved = interleaved_size(k, m, fwd_idx);

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
			uint64_t window_interleaved = interleaved_size(k, m, k-m-i);
			uint64_t rev_interleaved = interleaved_size(k, m, k-m-rev_idx);

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
	uint64_t fwd_interleaved = interleaved_size(k, m, fwd_idx);
	uint64_t rev_interleaved = interleaved_size(k, m, k-m-rev_idx);

	multiple_mini = true;
	mini_value = fwd_mini;
	leftmost_mini = min(fwd_idx, rev_idx);

	if (fwd_interleaved > rev_interleaved)
		return fwd_idx + start_idx;
	else if (rev_interleaved > fwd_interleaved)
		return -1 - rev_idx - start_idx;

	// Minimizer + position equality hell
	// remove_mini(uint8_t * seq, uint64_t seq_size, uint64_t start_idx, uint64_t k, uint64_t mini_pos, uint64_t m)
	uint8_t * fwd_seq = remove_mini(seq, seq_size, start_idx, k, fwd_idx, m);
	uint8_t * rev_seq = remove_mini(seq, seq_size, start_idx, k, rev_idx, m);
	this->rc.rev_comp(rev_seq, k-m);

	// Compare sequences
	int cmp_res = sequence_compare(
						fwd_seq, k-m, 0, k-m-1,
						rev_seq, k-m, 0, k-m-1
	);

	delete[] fwd_seq;
	delete[] rev_seq;

	if (cmp_res <= 0)
		return fwd_idx + start_idx;
	else
		return -1 - rev_idx - start_idx;
}


vector<skmer> MinimizerSearcher::get_skmers(const uint8_t * seq, const uint seq_size) {
	// Init the sliding window
	SlidingWindow window(this->m, seq, seq_size, 0, this->encoding);

	// Init the first minimizer value
	uint64_t minimizer = min(window.fwd, window.rev);
	int64_t mini_pos = window.fwd < window.rev ? 0 : -1;
	bool multi_mini = window.fwd == window.rev;
	uint64_t left_multi_idx = 0;

	// --- Init the first kmer minimizer ---
	for (uint64_t idx=1 ; idx<=k-m ; idx++)
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
	for(int64_t idx=k-m+1 ; idx<=seq_size-m ; idx++)
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
			skmers.push_back({skmer_start, static_cast<uint64_t>(idx+m-2), mini_pos, minimizer});
			// Update the values
			skmer_start = idx + m - k;
			minimizer = current_mini;
			mini_pos = current_mini_pos;
			multi_mini = current_double_mini;
			left_multi_idx = idx;
		}
		// Minimizer out of window
		else if ( ((mini_pos >= 0) and (idx >  mini_pos   + k - m))
			 or   ((mini_pos <  0) and (idx > -mini_pos-1 + k - m)) )
		{
			// Save the previous skmer
			skmers.push_back({skmer_start, static_cast<uint64_t>(idx+m-2), mini_pos, minimizer});
			skmer_start = idx + m - k;

			// Compute new minimizer values
			mini_pos = this->kmer_minimizer_compute(seq, seq_size, idx+m-k, minimizer, multi_mini, left_multi_idx);
		}
		// Same minimizer than previously
		else if (current_mini == minimizer)
		{
			// No change on minimizer value
			multi_mini = true;
			left_multi_idx = current_double_mini ? idx : mini_pos;
			
			current_mini_pos = this->kmer_minimizer_compute(seq, seq_size, idx+m-k);
			// Verify previous skmer continuity
			if (current_mini_pos != mini_pos)
			{
				// Add the previous skmer to the list
				skmers.push_back({skmer_start, static_cast<uint64_t>(idx+m-2), mini_pos, minimizer});
				// Update the values
				skmer_start = idx + m - k;
				mini_pos = current_mini_pos;
			}
		}
		// Cursed case. Recompute everything on each loop (small probability)
		else if (multi_mini)
		{
			// No minimizer value change
			current_mini_pos = this->kmer_minimizer_compute(seq, seq_size, idx+m-k);
			// Verify previous skmer continuity
			if (current_mini_pos != mini_pos)
			{
				// Add the previous skmer to the list
				skmers.push_back({skmer_start, static_cast<uint64_t>(idx+m-2), mini_pos, minimizer});
				// Update the values
				skmer_start = idx + m - k;
				mini_pos = current_mini_pos;
			}
		}
	}

	// Add last skmer
	skmers.push_back({skmer_start, static_cast<uint64_t>(seq_size-1), mini_pos, minimizer});

	return skmers;
}


vector<skmer> MinimizerSearcher::get_skmers_old(const uint8_t * seq, const uint seq_size) {
	this->compute_candidates(seq, seq_size);
	uint nb_kmers = seq_size - k + 1;
	this->compute_minimizers(nb_kmers);
	this->compute_skmers(nb_kmers);

	vector<skmer> skmers(this->skmers.size(), {0,0,0,0});

	for (uint sk_idx=0 ; sk_idx<this->skmers.size() ; sk_idx++) {
		uint64_t start = this->skmers[sk_idx].first;
		skmers[sk_idx].start_position = start;
		skmers[sk_idx].stop_position = this->skmers[sk_idx].second;
		int64_t mini_pos = this->mini_pos[start];
		skmers[sk_idx].minimizer_position = mini_pos;
		if (mini_pos >= 0)
			skmers[sk_idx].minimizer = this->mini_buffer[mini_pos];
		else {
			mini_pos = - mini_pos - 1;
			skmers[sk_idx].minimizer = this->mini_buffer[this->mini_buffer.size() / 2 + mini_pos];
		}
	}

	return skmers;
}


void MinimizerSearcher::compute_candidates(const uint8_t * seq, const uint seq_size) {
	if (seq_size > this->max_seq_size) {
		this->max_seq_size = seq_size;
		this->mini_buffer.resize((max_seq_size - m + 1) * 2, 0);
		this->minis.resize(this->max_seq_size - k + 1, 0);
		this->mini_pos.resize(this->max_seq_size - k + 1, 0);
	}

	uint offset = (4 - (seq_size % 4)) % 4;
	uint64_t current_value = 0;
	uint64_t current_rev_value = 0;
	// Compute prefix of first candidate
	for (uint i=0 ; i<this->m-1 ; i++) {
		uint idx = offset + i;
		uint byte_idx = idx/4;
		uint nucl_shift = 3 - (idx % 4);

		uint nucl = (seq[byte_idx] >> (nucl_shift * 2)) & 0b11;
		current_value = (current_value << 2) + nucl;
		current_rev_value = (current_rev_value >> 2) + (this->rc.reverse[nucl] << (2 * (this->m - 1)));
	}
	
	// Compute minimizer candidates
	uint64_t m_mask = (1 << (this->m*2)) - 1;
	for (uint i=this->m-1, kmer_idx=0 ; i<seq_size ; i++, kmer_idx++) {
		uint idx = offset + i;
		uint byte_idx = idx/4;
		uint nucl_shift = 3 - (idx % 4);

		uint nucl = (seq[byte_idx] >> (nucl_shift * 2)) & 0b11;
		current_value = ((current_value << 2) + nucl) & m_mask;
		current_rev_value = (current_rev_value >> 2) + this->nucl_rev[idx%4][seq[byte_idx]];//(this->rc.reverse[nucl] << (2 * (m - 1)));
		this->mini_buffer[kmer_idx] = current_value;
		this->mini_buffer[this->mini_buffer.size()/2 + kmer_idx] = current_rev_value;
    }
}


void MinimizerSearcher::compute_minimizers(const uint nb_kmers) {
    // Compute the minimizer of each sliding window of size k - m
    uint max_nb_candidates = this->mini_buffer.size()/2;
    for (uint i=0 ; i<nb_kmers ; i++) {
        auto smallest_fwd = min_element(
			this->mini_buffer.begin() + i,
			this->mini_buffer.begin()+i+(k-m)+1
		);
        auto smallest_rev = smallest_fwd;
        if (not this->single_side) {
            smallest_rev = min_element(
				this->mini_buffer.begin()+max_nb_candidates+i,
				this->mini_buffer.begin()+max_nb_candidates+i+(k-m)+1
			);
        }

        int64_t mini_pos = INT64_MAX;


		if (this->single_side or (*smallest_fwd) < (*smallest_rev)) {
			mini_pos = smallest_fwd - this->mini_buffer.begin();
		}

		else if ((*smallest_fwd) == (*smallest_rev)) {
			if (smallest_fwd < smallest_rev) {
				mini_pos = smallest_fwd - this->mini_buffer.begin();
			} else if (smallest_fwd > smallest_rev) {
				mini_pos = - (smallest_rev - (this->mini_buffer.begin() + this->mini_buffer.size() / 2)) - 1;
			} else { // smallest_fwd == smallest_rev, taking the smallest one
                for (uint64_t j = 0; j < k - m + 1; j++) {
                    if (mini_buffer[i + j] < mini_buffer[i + max_nb_candidates + j]) {
                        mini_pos = smallest_fwd - this->mini_buffer.begin();
                        break;
                    } else if (mini_buffer[i + j] > mini_buffer[i + max_nb_candidates + j]) {
                        mini_pos = - (smallest_rev - (this->mini_buffer.begin() + this->mini_buffer.size() / 2)) - 1;
                        break;
                    }
                }
            }
		}

        else {
            mini_pos = - (smallest_rev - (this->mini_buffer.begin() + this->mini_buffer.size() / 2)) - 1;
        }

        this->mini_pos[i] = mini_pos;
    }
}

void MinimizerSearcher::compute_skmers(const uint nb_kmers) {
	this->skmers.clear();

	uint last_mini_start = 0;
	for (uint idx=1 ; idx<nb_kmers ; idx++) {
		if (this->mini_pos[idx] != this->mini_pos[last_mini_start]) {
			this->skmers.emplace_back(last_mini_start, idx - 1 + this->k - 1);
			last_mini_start = idx;
		}
	}
	this->skmers.emplace_back(last_mini_start, nb_kmers - 1 + this->k - 1);
}


uint64_t seq_to_uint(const uint8_t * seq, uint seq_size) {
	uint nucl_to_extract = seq_size;
	if (nucl_to_extract > 32)
		nucl_to_extract = 32;

	uint seq_offset = (4 - (seq_size % 4)) % 4;
	uint seq_bytes = (seq_size + 3) / 4;
	uint useless_seq_nucl = seq_size - nucl_to_extract;

	uint suff_offset = (4 - (nucl_to_extract % 4)) % 4;
	uint mask = (1 << (2 * (4 - suff_offset))) - 1;
	uint suff_first_byte = (seq_offset + useless_seq_nucl) / 4;

	uint64_t val = seq[suff_first_byte] & mask;
	for (uint idx=suff_first_byte+1 ; idx<seq_bytes ; idx++) {
		val <<= 8;
		val += seq[idx];
	}

	return val;
}


uint64_t subseq_to_uint(const uint8_t * seq, uint seq_size, uint start_nucl, uint end_nucl) {
	// Trunkate too long sequences
	if (end_nucl - start_nucl + 1 > 32)
		start_nucl = end_nucl - 31;

	// Determine boundary bytes
	uint seq_offset = (4 - (seq_size % 4)) % 4;
	uint first_sub_byte = (seq_offset + start_nucl) / 4;
	uint last_sub_byte = (seq_offset + end_nucl) / 4;

	if (first_sub_byte == last_sub_byte) {
		uint64_t val = seq[first_sub_byte];
		uint shift = 2 * (3 - seq_offset - end_nucl);
		uint mask = (1 << ((end_nucl - start_nucl + 1) * 2)) - 1;
		return (val >> shift) & mask;
	}

	// First byte integration
	uint mask = (1 << (2 * (4 - ((seq_offset + start_nucl) % 4)))) - 1;
	uint64_t sub_val = seq[first_sub_byte] & mask;

	// Middle bytes
	for (uint b=first_sub_byte+1 ; b<last_sub_byte ; b++) {
		sub_val <<= 8;
		sub_val += seq[b];
	}

	// End byte
	uint last_shift = (seq_size - end_nucl - 1) % 4;
	uint8_t end_byte = seq[last_sub_byte] >> (2 * last_shift);
	sub_val <<= 2 * (4 - last_shift);
	sub_val += end_byte;

	return sub_val;
}


void uint_to_seq(uint seq, uint8_t * bin_seq, uint size) {
	uint seq_bytes = (size + 3) / 4;

	for (int idx=seq_bytes-1 ; idx>=0 ; idx--) {
		bin_seq[idx] = seq & 0b11111111;
		seq >>= 8;
	}
}