import os
import subprocess
import tempfile
import unittest


class TestShuffle(unittest.TestCase):
    """shuffle keeps the kmers and their data, with or without data."""

    def check(self, lines, data_size):
        with tempfile.TemporaryDirectory() as tmp:
            txt = os.path.join(tmp, "in.txt")
            with open(txt, "w") as f:
                f.write("\n".join(lines) + "\n")
            kff = os.path.join(tmp, "in.kff")
            out = os.path.join(tmp, "out.kff")
            subprocess.run(["./bin/kff-tools", "instr", "-i", txt, "-o", kff, "-k", "7",
                            "-d", str(data_size), "-m", "10"], check=True)
            subprocess.run(["./bin/kff-tools", "shuffle", "-i", kff, "-o", out], check=True,
                           stderr=subprocess.DEVNULL)
            kmers = lambda path: sorted(subprocess.run(["./bin/kff-tools", "outstr", "-i", path],
                                        check=True, capture_output=True, text=True).stdout.splitlines())
            self.assertEqual(kmers(kff), kmers(out))

    def test_without_data(self):
        self.check(["ACGTACGTAC", "TTTGGGCCAA", "ACGTTGCA"], 0)

    def test_with_data(self):
        self.check(["ACGTACGTAC 1,2,3,4", "TTTGGGCCAA 5,6,7,8", "ACGTTGCA 9,10"], 1)


if __name__ == "__main__":
    unittest.main()
