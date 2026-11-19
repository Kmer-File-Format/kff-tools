import os
import subprocess
import tempfile
import unittest


class TestInstrSplitLines(unittest.TestCase):
    """A line with more than max kmers is written in several blocks: each block
    must keep the data of its own kmers."""

    def run_instr(self, lines, k, max_kmers, data_size):
        with tempfile.TemporaryDirectory() as tmp:
            txt = os.path.join(tmp, "in.txt")
            kff = os.path.join(tmp, "out.kff")
            with open(txt, "w") as f:
                f.write("".join(line + "\n" for line in lines))
            subprocess.run(["./bin/kff-tools", "instr", "-i", txt, "-o", kff, "-k", str(k),
                            "-m", str(max_kmers), "-d", str(data_size)], check=True)
            out = subprocess.run(["./bin/kff-tools", "outstr", "-i", kff], check=True,
                                 capture_output=True, text=True).stdout
        return out.splitlines()

    def test_counts_of_split_lines(self):
        seq = "ACGTACGTACCA"   # 8 kmers of size 5
        counts = [11, 12, 13, 14, 15, 16, 17, 18]
        lines = [seq + " " + ",".join(map(str, counts)), "TTGCA 99"]
        for max_kmers in (1, 2, 3, 8, 20):
            out = self.run_instr(lines, 5, max_kmers, 1)
            expected = [f"{seq[i:i + 5]} {c}" for i, c in enumerate(counts)] + ["TTGCA 99"]
            self.assertEqual(out, expected, f"max = {max_kmers}")

    def test_lines_shorter_than_k_are_skipped(self):
        out = self.run_instr(["ACG 1", "ACGTA 7"], 5, 4, 2)
        self.assertEqual(out, ["ACGTA 7"])


if __name__ == "__main__":
    unittest.main()
