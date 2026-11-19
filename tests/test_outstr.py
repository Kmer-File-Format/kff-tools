import os
import subprocess
import tempfile
import unittest


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


class TestOutstrCanonical(unittest.TestCase):
    """outstr -c prints the smallest of a kmer and its reverse complement in the
    encoding order of the file (A=0 C=1 T=2 G=3 for instr), for any k."""

    def test_canonical_kmers(self):
        order = {"A": "0", "C": "1", "T": "2", "G": "3"}
        key = lambda s: "".join(order[c] for c in s)
        kmers = ["GGGGAAAA", "TTTTCCCC", "ACGTTGCA", "GCATGCAT", "CCCCAAAA", "AGTCAGTC"]
        for k in (4, 7, 8, 9, 12):
            seqs = [(kmer * 2)[:k] for kmer in kmers]
            with tempfile.TemporaryDirectory() as tmp:
                txt = os.path.join(tmp, "in.txt")
                with open(txt, "w") as f:
                    f.write("\n".join(seqs) + "\n")
                kff = os.path.join(tmp, "in.kff")
                subprocess.run(["./bin/kff-tools", "instr", "-i", txt, "-o", kff, "-k", str(k), "-m", "1"], check=True)
                out = subprocess.run(["./bin/kff-tools", "outstr", "-c", "-i", kff], check=True,
                                     capture_output=True, text=True).stdout.split()
            expected = [min(s, revcomp(s), key=key) for s in seqs]
            self.assertEqual(out, expected, f"k = {k}")


if __name__ == "__main__":
    unittest.main()
