#include <cstring>
#include <cassert>
#include <algorithm>

#include "sequences.hpp"



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


uint KffSeqStream::next_sequence(uint8_t * & seq, uint8_t * & data) {
	if (this->reader.has_next()) {
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


int sequence_compare(const uint8_t * seq1, const uint seq1_size,
											const uint seq1_start, const uint seq1_stop,
											const uint8_t * seq2, const uint seq2_size,
											const uint seq2_start, const uint seq2_stop) {
	// If inequal size
	if (seq1_stop - seq1_start != seq2_stop - seq2_start)
		return (seq1_stop - seq1_start) < (seq2_stop - seq2_start) ? -1 : 1;

	// Extraction of subsequences
	uint subseq_size = seq1_stop - seq1_start + 1;
	uint subseq_bytes = (subseq_size + 3) / 4;
	uint8_t * subseq1 = new uint8_t[subseq_bytes];
	subsequence(seq1, seq1_size, subseq1, seq1_start, seq1_stop);
	uint8_t * subseq2 = new uint8_t[subseq_bytes];
	subsequence(seq2, seq2_size, subseq2, seq2_start, seq2_stop);

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


Minimizer_Creator::Minimizer_Creator(const uint64_t k, const uint64_t m, RevComp & rc)
										: k(k)
										, m(m)
										, r(rc)
{
	this->candidate_vector_size = 1;
	this->candidates_fwd = new uint64_t[1];
  this->candidates_rev = new uint64_t[1];
};


Minimizer_Creator::~Minimizer_Creator() {
	delete[] this->candidates_fwd;
	delete[] this->candidates_rev;
}


void Minimizer_Creator::shift4_compute_mini_candidates(const uint8_t * seq, const uint size) {
	uint64_t tmp_mini_candidates[128];

	uint8_t encoding[] = {0, 1, 3, 2};
	Stringifyer strif(encoding);
	cout << strif.translate(seq, size) << " " << endl;

	uint8_t * shift_buffer;
	uint64_t seq_bytes = (size + 3) / 4;
	uint64_t left_missing_bytes = 8 - (m + 3) / 4;
	uint64_t buffer_size = seq_bytes + left_missing_bytes;
	cout << "buffer size " << buffer_size << endl;
	shift_buffer = new uint8_t[buffer_size];
	memset(shift_buffer, 0, buffer_size);

	// Copy the sequence into the buffer
	memcpy(shift_buffer + left_missing_bytes, seq, seq_bytes);

	uint64_t nb_mini_candidates = size - this->m + 1;
	uint64_t mini_mask = (1 << (2 * this->m)) - 1;
	cout << "mask " << mini_mask << endl;
	// Read initialy aligned candidates
	for (int64_t m_idx= nb_mini_candidates-1 ; m_idx>=0 ; m_idx -= 4) {
		cout << "m_idx " << m_idx << " nb_mini_candidates " << nb_mini_candidates << endl;
		uint64_t right_m_idx = nb_mini_candidates - 1 - m_idx;
		cout << "right idx " << right_m_idx << endl;
		cout << "position " << buffer_size - 8 - right_m_idx/4 << endl;
		uint8_t * buffer_position = shift_buffer + buffer_size - 8 - right_m_idx/4;
		cout << (uint64_t)*(buffer_position+7) << endl;
		uint64_t candidate = (*((uint64_t *)buffer_position)) & mini_mask;
		cout << m_idx << ": " << candidate << " " << this->candidates_fwd[m_idx] << endl;
	}

	for (uint shift=2 ; shift<8 ; shift += 2) {
		// Shift the array
		// Read aligned candidates
	}

	exit(0);
}

void Minimizer_Creator::naive_compute_mini_candidates(const uint8_t * seq, const uint size) {
	if (size > this->candidate_vector_size) {
		delete[] this->candidates_fwd;
		delete[] this->candidates_rev;

		this->candidate_vector_size = size;
		this->candidates_fwd = new uint64_t[size];
		this->candidates_rev = new uint64_t[size];
	}

	uint offset = (4 - (size % 4)) % 4;
	uint64_t current_value = 0;
	uint64_t current_rev_value = 0;
	// Compute prefix of first candidate
	for (uint i=0 ; i<m-1 ; i++) {
		uint idx = offset + i;
		uint byte_idx = idx/4;
		uint nucl_shift = 3 - (idx % 4);

		uint nucl = (seq[byte_idx] >> (nucl_shift * 2)) & 0b11;
		current_value = (current_value << 2) + nucl;
		// cout << nucl << "->" << (uint)r.reverse[nucl] << endl;
		current_rev_value = (current_rev_value >> 2) + (r.reverse[nucl] << (2 * (m - 1)));
	}

	// Compute minimizer candidates
	uint64_t m_mask = (1 << (m*2)) - 1;
	for (uint i=m-1 ; i<size ; i++) {
		uint idx = offset + i;
		uint byte_idx = idx/4;
		uint nucl_shift = 3 - (idx % 4);

		uint nucl = (seq[byte_idx] >> (nucl_shift * 2)) & 0b11;
		current_value = ((current_value << 2) + nucl) & m_mask;
		current_rev_value = (current_rev_value >> 2) + (r.reverse[nucl] << (2 * (m - 1)));
		this->candidates_fwd[i - m + 1] = current_value;
		this->candidates_rev[i - m + 1] = current_rev_value;
	}

	//this->shift4_compute_mini_candidates(seq, size);
}

vector<pair<int, uint64_t> > Minimizer_Creator::compute_minizers(const uint8_t * seq, const uint size, const bool single_side) {
	vector<pair<int, uint64_t> > minimizers;
	// Get all the candidates
	this->naive_compute_mini_candidates(seq, size);

	int prev_pos = k + 2;
	// Compute the minimizer of each sliding window of size k - m
	for (uint i=0 ; i<=size-k ; i++) {
		auto smallest_fwd = min_element(
			this->candidates_fwd + i,
			this->candidates_fwd + i + (k - m) + 1
		);
		auto smallest_rev = smallest_fwd;
		if (not single_side)
			smallest_rev = min_element(
				this->candidates_rev + i,
				this->candidates_rev + i + (k - m) + 1
			);

		int mini_pos;
		uint64_t minimizer;
		if (single_side or *smallest_fwd < *smallest_rev) {
			mini_pos = smallest_fwd - this->candidates_fwd;
			minimizer = *smallest_fwd;
		} 
		else if (*smallest_fwd == *smallest_rev) {
			if (smallest_fwd - smallest_rev <= 0) {
				mini_pos = smallest_fwd - this->candidates_fwd;
				minimizer = *smallest_fwd;
			} else {
				mini_pos = - (this->candidates_rev - smallest_rev);
				minimizer = *smallest_rev;
			}
		}
		else {
			mini_pos = - (this->candidates_rev - smallest_rev);
			minimizer = *smallest_rev;
		}

		// New minimizer ?
		if (mini_pos != prev_pos) {
			prev_pos = mini_pos;
			minimizers.emplace_back(mini_pos, minimizer);
		}
	}

	return minimizers;
}

std::vector<pair<int, int> > Minimizer_Creator::compute_skmers(const uint seq_size, std::vector<std::pair<int, uint64_t> > & minimizers) {
	// Superkmer list
	vector<pair<int, int> > skmers;

	uint current_begin = 0;
	for (uint i=0 ; i<minimizers.size()-1 ; i++) {
		pair<int, uint64_t> & current_mini = minimizers[i];
		pair<int, uint64_t> & next_mini = minimizers[i+1];

		// cout << "current_mini " << current_mini.first << " " << current_mini.second << endl;
		// cout << "next_mini " << next_mini.first << " " << next_mini.second << endl;


		uint end = 0;
		uint new_begin = 0;
		bool is_fwd = true;
		// Minimizer position
		int mini_pos = current_mini.first;
		// Reverse
		if (current_mini.first < 0) {
			is_fwd = false;
			mini_pos = seq_size + current_mini.first - m + 1;
			// cout << "absolute mini_pos " << mini_pos << endl;
		}

		// First minimizer is dominant
		if (current_mini.second <= next_mini.second) {
			// cout << "dominant" << endl;
			end = mini_pos + k - 1;
			new_begin = mini_pos + 1;
			// cout << "end " << end << " new_begin " << new_begin << endl;
		}
		// Second minimizer is dominant
		else {
			int next_mini_pos = next_mini.first;
			if (next_mini_pos < 0)
				next_mini_pos = seq_size + next_mini_pos - m + 1;
			new_begin = next_mini_pos + m - k;
			end = new_begin + k - 2;
		}

		// Add the superkmer and save position
		if (is_fwd)
			skmers.emplace_back(current_begin, end);
		else
			skmers.emplace_back(-(seq_size - end), -(seq_size - current_begin));
		current_begin = new_begin;
	}

	// Save the last skmer
	pair<int, uint64_t> & last_mini = *(minimizers.end()-1);
	if (last_mini.first >= 0)
		skmers.emplace_back(current_begin, seq_size-1);
	else
		skmers.emplace_back(-1, -(seq_size - current_begin));

	return skmers;
}

std::vector<pair<int, int> > Minimizer_Creator::compute_skmers(const uint8_t * seq, const uint size,  const bool single_side) {
	std::vector<std::pair<int, uint64_t> > minimizers = compute_minizers(seq, size, single_side);

	return compute_skmers(size, minimizers);
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
