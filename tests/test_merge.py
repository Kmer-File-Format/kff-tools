import unittest
import os
import subprocess

import kmer_generation as kg


class TestMerge(unittest.TestCase):

    def test_merge_index(self):
        print(f"\n-- TestMergeSplit test_merge_index")
        print("  init - generate 3 random kmer files")
        txt_file_1 = f"txt1_test.txt"
        kff_file_1 = f"txt1_test.kff"
        txt_file_2 = f"txt2_test.txt"
        kff_file_2 = f"txt2_test.kff"
        kg.generate_sequences_file(txt_file_1, 10, 32, size_max=42)
        kg.generate_random_kmers_file(txt_file_2, 10, 17)

        # Prepare a file with sections
        print(f"  1/3 Generate kff files from txts.")
        self.assertEqual(0, os.system(f"./bin/kff-tools instr -i {txt_file_1} -o {kff_file_1} -k 32 -m 11"))
        # Read 1 kmer per line
        self.assertEqual(0, os.system(f"./bin/kff-tools instr -i {txt_file_2} -o {kff_file_2} -k 17"))

        print(f"  2/3 Merge files and split it again")
        merged = "merged_raw_test.kff"
        self.assertEqual(0, os.system(f"./bin/kff-tools merge -i {kff_file_1} {kff_file_2} -o {merged}"))
        index_output = []
        try:
            # Analyse the outfile
            output = subprocess.check_output(f"./bin/kff-tools validate --infile {merged} -v", stderr=subprocess.STDOUT, shell=True, text=True)

            # extract index output
            is_index_section = False
            for line in output.split("\n"):
                if line.startswith("=== Section i ==="):
                    is_index_section = True
                elif line.startswith("=== "):
                    is_index_section = False

                if is_index_section:
                    index_output += [line]
        except CalledProcessError:
            self.fail("Error raised on merged file validation test")

        # Analyse the index of the outfile
        self.assertTrue(index_output[3].startswith('v'), f"found: {index_output[3]}")
        self.assertTrue(index_output[4].startswith('r'), f"found: {index_output[4]}")
        self.assertTrue(index_output[5].startswith('v'), f"found: {index_output[5]}")
        self.assertTrue(index_output[6].startswith('r'), f"found: {index_output[6]}")

        print("  clean the test area")
        os.system(f"rm -r {merged} {txt_file_1} {kff_file_1} {txt_file_2} {kff_file_2}")



if __name__ == '__main__':
  unittest.main()
