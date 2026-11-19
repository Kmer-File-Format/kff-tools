import os
import random
import subprocess
import tempfile
import unittest


class TestCompactKeepsKmers(unittest.TestCase):
    """compact (greedy and sorted) keeps exactly the kmers of its input, including
    large k (k - m >= 128, minimizer positions on 2 bytes in the output)."""

    def kmers(self, path):
        out = subprocess.run(["./bin/kff-tools", "outstr", "-i", path], check=True,
                             capture_output=True, text=True).stdout
        return sorted(line.split()[0] for line in out.splitlines())

    def check(self, k, m, length):
        rnd = random.Random(k)
        genome = "".join(rnd.choice("ACGT") for _ in range(length))
        with tempfile.TemporaryDirectory() as tmp:
            txt = os.path.join(tmp, "kmers.txt")
            with open(txt, "w") as f:
                f.write("".join(genome[i:i + k] + "\n" for i in range(len(genome) - k + 1)))
            kff = os.path.join(tmp, "kmers.kff")
            bucket = os.path.join(tmp, "bucket.kff")
            subprocess.run(["./bin/kff-tools", "instr", "-i", txt, "-o", kff, "-k", str(k), "-m", "1"], check=True)
            subprocess.run(["./bin/kff-tools", "bucket", "-i", kff, "-o", bucket, "-m", str(m)], check=True)
            expected = self.kmers(bucket)
            for options in ([], ["-s"]):
                out = os.path.join(tmp, "compacted.kff")
                subprocess.run(["./bin/kff-tools", "compact", "-i", bucket, "-o", out] + options, check=True,
                               stderr=subprocess.DEVNULL)
                self.assertEqual(self.kmers(out), expected, f"k={k} m={m} {options}")

    def test_small_k(self):
        self.check(15, 10, 3000)

    def test_large_k(self):
        self.check(141, 11, 2000)


if __name__ == "__main__":
    unittest.main()
