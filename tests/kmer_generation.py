import random


def generate_sequences(nb_seq, size_min, size_max=0):
    if size_max == 0:
        size_max = size_min

    nucleotides = "ACGT"
    for _ in range(nb_seq):
        yield ''.join(random.choice(nucleotides) for __ in range(random.randint(size_min, size_max)))


def generate_counts(n, max_val):
    return [str(random.randint(1, max_val)) for _ in range(n)]


def generate_sequences_file(filename, nb_seq, size_min, size_max=0, max_count=0):
    with open(filename, "w") as fp:
        for seq in generate_sequences(nb_seq, size_min, size_max):
            if max_count == 0:
                fp.write(f"{seq}\n")
            else:
                fp.write(f"{seq} {','.join(generate_counts(len(seq) - size_min + 1, max_count))}\n")


def generate_random_kmers_file(filename, nb_kmers, k, max_count=0, overlapping=False):
    with open(filename, "w") as fp:
        generated = 0

        while (generated < nb_kmers):
            seq = None
            seq_kmers = 0

            if overlapping:
                seq = generate_sequences(1, k, 4*k).__next__()
                seq_kmers = min(len(seq) - k + 1, nb_kmers - generated)
            else:
                seq = generate_sequences(1, k, k).__next__()
                seq_kmers = 1

            for i in range(seq_kmers):
                fp.write(seq[i:i+k])

                if max_count > 0:
                    fp.write(f" {random.randint(1, max_count)}")
                fp.write("\n")

            generated += seq_kmers




if __name__ == "__main__":
    generate_sequences_file("mytest.txt", 10, 32, 39, 255);
