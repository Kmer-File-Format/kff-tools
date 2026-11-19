import os
import random
import subprocess
import tempfile
import unittest


class TestSortReproducible(unittest.TestCase):
    """Sorting the same file twice gives the same file, even with equal
    sequences whose size is a multiple of 4."""

    def test_same_output_twice(self):
        rnd = random.Random(3)
        seqs = ["".join(rnd.choice("ACGT") for _ in range(rnd.choice([12, 16, 13])))
                for _ in range(30)]
        lines = []
        for _ in range(3000):
            s = rnd.choice(seqs)
            lines.append(s + " " + ",".join(str(rnd.randint(0, 255)) for _ in range(len(s) - 10)))
        with tempfile.TemporaryDirectory() as tmp:
            txt = os.path.join(tmp, "in.txt")
            with open(txt, "w") as f:
                f.write("\n".join(lines) + "\n")
            kff = os.path.join(tmp, "in.kff")
            subprocess.run(["./bin/kff-tools", "instr", "-i", txt, "-o", kff, "-k", "11",
                            "-d", "1", "-m", "10"], check=True)
            outputs = []
            for i in range(3):
                out = os.path.join(tmp, f"sorted{i}.kff")
                subprocess.run(["./bin/kff-tools", "sort", "-i", kff, "-o", out], check=True,
                               stderr=subprocess.DEVNULL)
                with open(out, "rb") as f:
                    outputs.append(f.read())
            self.assertEqual(outputs[0], outputs[1])
            self.assertEqual(outputs[0], outputs[2])
            # Same kmers and data
            kmers = lambda path: sorted(subprocess.run(["./bin/kff-tools", "outstr", "-i", path],
                                        check=True, capture_output=True, text=True).stdout.splitlines())
            self.assertEqual(kmers(kff), kmers(os.path.join(tmp, "sorted0.kff")))


if __name__ == "__main__":
    unittest.main()
