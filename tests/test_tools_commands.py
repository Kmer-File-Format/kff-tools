import unittest
import os

import kmer_generation as kg


# file = "kmers"
file = "ecoli_count_dsk"


class TestInOut(unittest.TestCase):
    def test_raw_sections(self):
        print("\n-- TestInOut test_raw_section")
        for k in range(12, 16):
            # Generate random kmer
            print(f"Generate a random {k}-mers file.")
            txt_file = f"inout_raw_k{k}_test.txt"
            kg.generate_random_kmers_file(txt_file, 100, k, max_count=511, overlapping=True)

            # Generate a kff file from a textual kmer count
            print("  1/3 Generate a kff file from a textual kmer count")
            kff_file = f"inout_raw_k{k}_test.kff"
            self.assertEqual(0, os.system(f"./bin/kff-tools instr --kmer-size {k} --data-size 2 --infile {txt_file} --outfile {kff_file}"))

            # Regenerate a textual file from the kff
            print("  2/3 Regenerate a txt file from the kff")
            txt_out_file = f"inout_final_k{k}_test.kff"
            self.assertEqual(0, os.system(f"./bin/kff-tools outstr --infile {kff_file} > {txt_out_file}"))

            # Compate the original file to the final translated file
            print("  3/3 Compare initial and final txt files")
            stream = os.popen(f"diff {txt_file} {txt_out_file}")
            self.assertEqual(stream.read(), "")
            stream.close()

            print("  Clean the directory")
            self.assertEqual(0, os.system(f"rm {txt_file} {kff_file} {txt_out_file}"))

    def test_data_sequences(self):
        print("\n-- TestInOut test_data_sequences")
        # Create a test file
        print("Generate a data sequence test file")
        seqfilename = "seqfile_test.txt";
        with open(seqfilename, "w") as seqfile:
            seqfile.write("AGTTCT 12,3\n")
            seqfile.write("GAGCT 4\n")
            seqfile.write("TCTTACC 1,2,7\n")

        compfilename = "compfile_test.txt"
        with open(compfilename, "w") as compfile:
            compfile.write("AGTTC 12\nGTTCT 3\n")
            compfile.write("GAGCT 4\n")
            compfile.write("TCTTA 1\nCTTAC 2\nTTACC 7\n")

        # Convert file into kff
        print("1/3 Generate a kff file from the txt")
        kff_file = seqfilename[:-4] + ".kff"
        print(kff_file)
        self.assertEqual(0, os.system(f"./bin/kff-tools instr --kmer-size 5 --data-size 1 --infile {seqfilename} --outfile {kff_file}"))

        print("2/3 Output the kff file as kmer list")
        outfile = seqfilename + ".out"
        self.assertEqual(0, os.system(f"./bin/kff-tools outstr --infile {kff_file} > {outfile}"))

        print("3/3 Check the kmers")
        stream = os.popen(f"diff {compfilename} {outfile}")
        diff = stream.read()
        stream.close()
        self.assertEqual(diff, "")

        print("  Clean the directory")
        self.assertEqual(0, os.system(f"rm {seqfilename} {kff_file} {outfile} {compfilename}"))


class TestSplitMerge(unittest.TestCase):

    def test_merge_split(self):
        print(f"\n-- TestMergeSplit test_merge_split")
        print("  init - generate 3 random kmer files")
        txt_file_1 = f"txt1_test.txt"
        kff_file_1 = f"txt1_test.kff"
        txt_file_2 = f"txt2_test.txt"
        kff_file_2 = f"txt2_test.kff"
        txt_file_3 = f"txt3_test.txt"
        kff_file_3 = f"txt3_test.kff"
        kg.generate_sequences_file(txt_file_1, 100, 32, size_max=42)
        kg.generate_random_kmers_file(txt_file_2, 100, 17)
        kg.generate_random_kmers_file(txt_file_3, 100, 33, max_count=511, overlapping=False)

        # Prepare a file with sections
        print(f"  1/3 Generate kff files from txts.")
        self.assertEqual(0, os.system(f"./bin/kff-tools instr -i {txt_file_1} -o {kff_file_1} -k 32 -m 11"))
        # Read 1 kmer per line
        self.assertEqual(0, os.system(f"./bin/kff-tools instr -i {txt_file_2} -o {kff_file_2} -k 17"))
        # Read 1 kmer per line with its counts (up to 255)
        self.assertEqual(0, os.system(f"./bin/kff-tools instr -i {txt_file_3} -o {kff_file_3} -k 33 -d 2"))


        print(f"  2/3 Merge files and split it again")
        merged = "merged_raw_test.kff"
        split_dir = "split_test/"
        self.assertEqual(0, os.system(f"rm -r {split_dir} ; mkdir {split_dir}"))
        self.assertEqual(0, os.system(f"./bin/kff-tools merge -i {kff_file_1} {kff_file_2} {kff_file_3} -o {merged}"))
        self.assertEqual(0, os.system(f"./bin/kff-tools split --infile {merged} --outdir {split_dir}"))
        

        print(f"  3/3 Compare outputs")
        self.assertEqual(0, os.system(f"./bin/kff-tools outstr -i {kff_file_1} > {split_dir}/{txt_file_1}_original"))
        self.assertEqual(0, os.system(f"./bin/kff-tools outstr -i {split_dir}/r_0.kff > {split_dir}/{txt_file_1}_split"))
        stream = os.popen(f"diff {split_dir}/{txt_file_1}_original {split_dir}/{txt_file_1}_split")
        stream_val = stream.read()
        stream.close()
        self.assertEqual(stream_val, "")

        self.assertEqual(0, os.system(f"./bin/kff-tools outstr -i {kff_file_2} > {split_dir}/{txt_file_2}_original"))
        self.assertEqual(0, os.system(f"./bin/kff-tools outstr -i {split_dir}/r_1.kff > {split_dir}/{txt_file_2}_split"))
        stream = os.popen(f"diff {split_dir}/{txt_file_2}_original {split_dir}/{txt_file_2}_split")
        stream_val = stream.read()
        stream.close()
        self.assertEqual(stream_val, "")

        self.assertEqual(0, os.system(f"./bin/kff-tools outstr -i {kff_file_3} > {split_dir}/{txt_file_3}_original"))
        self.assertEqual(0, os.system(f"./bin/kff-tools outstr -i {split_dir}/r_2.kff > {split_dir}/{txt_file_3}_split"))
        stream = os.popen(f"diff {split_dir}/{txt_file_3}_original {split_dir}/{txt_file_3}_split")
        stream_val = stream.read()
        stream.close()
        self.assertEqual(stream_val, "")


        print("  clean the test area")
        os.system(f"rm -r {merged} {split_dir} {txt_file_1} {kff_file_1} {txt_file_2} {kff_file_2} {txt_file_3} {kff_file_3}")


    def test_merge_split_file(self):
        print(f"\n-- TestMergeSplit test_merge_split")
        print("  init - generate 3 random kmer files")
        txt_file_1 = f"txt1_test.txt"
        kff_file_1 = f"txt1_test.kff"
        txt_file_2 = f"txt2_test.txt"
        kff_file_2 = f"txt2_test.kff"
        txt_file_3 = f"txt3_test.txt"
        kff_file_3 = f"txt3_test.kff"
        kg.generate_sequences_file(txt_file_1, 100, 32, size_max=42)
        kg.generate_random_kmers_file(txt_file_2, 100, 17)
        kg.generate_random_kmers_file(txt_file_3, 100, 33, max_count=511, overlapping=False)

        # Prepare a file with sections
        print(f"  1/3 Generate kff files from txts.")
        self.assertEqual(0, os.system(f"./bin/kff-tools instr -i {txt_file_1} -o {kff_file_1} -k 32 -m 11"))
        # Read 1 kmer per line
        self.assertEqual(0, os.system(f"./bin/kff-tools instr -i {txt_file_2} -o {kff_file_2} -k 17"))
        # Read 1 kmer per line with its counts (up to 255)
        self.assertEqual(0, os.system(f"./bin/kff-tools instr -i {txt_file_3} -o {kff_file_3} -k 33 -d 2"))

        filelist = "filelist_test.txt"
        with open(filelist, "w") as fp:
            print(kff_file_1, file=fp)
            print(kff_file_2, file=fp)
            print(kff_file_3, file=fp)


        print(f"  2/3 Merge files and split it again")
        merged = "merged_raw_test.kff"
        split_dir = "split_test/"
        self.assertEqual(0, os.system(f"rm -r {split_dir} ; mkdir {split_dir}"))
        self.assertEqual(0, os.system(f"valgrind ./bin/kff-tools merge -f {filelist} -o {merged}"))
        self.assertEqual(0, os.system(f"./bin/kff-tools split --infile {merged} --outdir {split_dir}"))
        

        print(f"  3/3 Compare outputs")
        self.assertEqual(0, os.system(f"./bin/kff-tools outstr -i {kff_file_1} > {split_dir}/{txt_file_1}_original"))
        self.assertEqual(0, os.system(f"./bin/kff-tools outstr -i {split_dir}/r_0.kff > {split_dir}/{txt_file_1}_split"))
        stream = os.popen(f"diff {split_dir}/{txt_file_1}_original {split_dir}/{txt_file_1}_split")
        stream_val = stream.read()
        stream.close()
        self.assertEqual(stream_val, "")

        self.assertEqual(0, os.system(f"./bin/kff-tools outstr -i {kff_file_2} > {split_dir}/{txt_file_2}_original"))
        self.assertEqual(0, os.system(f"./bin/kff-tools outstr -i {split_dir}/r_1.kff > {split_dir}/{txt_file_2}_split"))
        stream = os.popen(f"diff {split_dir}/{txt_file_2}_original {split_dir}/{txt_file_2}_split")
        stream_val = stream.read()
        stream.close()
        self.assertEqual(stream_val, "")

        self.assertEqual(0, os.system(f"./bin/kff-tools outstr -i {kff_file_3} > {split_dir}/{txt_file_3}_original"))
        self.assertEqual(0, os.system(f"./bin/kff-tools outstr -i {split_dir}/r_2.kff > {split_dir}/{txt_file_3}_split"))
        stream = os.popen(f"diff {split_dir}/{txt_file_3}_original {split_dir}/{txt_file_3}_split")
        stream_val = stream.read()
        stream.close()
        self.assertEqual(stream_val, "")


        print("  clean the test area")
        os.system(f"rm -r {merged} {split_dir} {txt_file_1} {kff_file_1} {txt_file_2} {kff_file_2} {txt_file_3} {kff_file_3} {filelist}")


class TestBucketting(unittest.TestCase):

    def test_basic_bucketting(self):
        print(f"\n-- TestBucketting - basic bucketting")
        print("  init - generate a random sequence file")
        txt = f"txt_test.txt"
        kff_raw = f"kff_raw_test.kff"
        kff_bucket = f"kff_bucket_test.kff"
        kg.generate_sequences_file(txt, 100, 32, size_max=42)
        # kg.generate_sequences_file(txt, 1, 32, size_max=32)

        # Prepare a file with sections
        print(f"  1/3 Generate the kff raw file.")
        self.assertEqual(0, os.system(f"./bin/kff-tools instr -i {txt} -o {kff_raw} -k 32 -m 11"))


        print(f"  2/3 bucket the file")
        print(f"./bin/kff-tools bucket -i {kff_raw} -o {kff_bucket} -m 11")
        self.assertEqual(0, os.system(f"./bin/kff-tools bucket -i {kff_raw} -o {kff_bucket} -m 11"))
        

        print(f"  3/3 Compare outputs")
        self.assertEqual(0, os.system(f"./bin/kff-tools outstr -c -i {kff_raw} | sort > {kff_raw}_sorted.txt"))
        self.assertEqual(0, os.system(f"./bin/kff-tools outstr -c -i {kff_bucket} | sort > {kff_bucket}_sorted.txt"))
        stream = os.popen(f"diff {kff_raw}_sorted.txt {kff_bucket}_sorted.txt")
        stream_val = stream.read()
        stream.close()
        self.assertEqual(stream_val, "")


        print("  clean the test area")
        os.system(f"rm -r {txt} {kff_raw}* {kff_bucket}*")



class TestCompaction(unittest.TestCase):

    def test_greedy_compaction(self):
        print(f"\n-- TestCompaction - greedy compaction")
        print("  init - generate a random sequence file")
        txt = "test.txt"
        kff_raw = "raw_test.kff"
        kff_disjoin = "disjoin_test.kff"
        kff_bucket = "bucket_test.kff"
        kff_compacted = "compact_test.kff"
        kg.generate_sequences_file(txt, 100, 32, size_max=42)

        # Prepare a file with sections
        print(f"  1/4 Generate the kff raw file.")
        self.assertEqual(0, os.system(f"./bin/kff-tools instr -i {txt} -o {kff_raw} -k 32 -m 11"))
        self.assertEqual(0, os.system(f"./bin/kff-tools disjoin -i {kff_raw} -o {kff_disjoin}"))

        print(f"  2/4 Bucket the file")
        self.assertEqual(0, os.system(f"./bin/kff-tools bucket -i {kff_disjoin} -o {kff_bucket} -m 11"))

        print(f"  3/4 Compact kmers")
        self.assertEqual(0, os.system(f"./bin/kff-tools compact -i {kff_bucket} -o {kff_compacted}"))
        

        print(f"  4/4 Compare outputs")
        self.assertEqual(0, os.system(f"./bin/kff-tools outstr -c -i {kff_raw} | sort > {kff_raw}_sorted.txt"))
        self.assertEqual(0, os.system(f"./bin/kff-tools outstr -c -i {kff_compacted} | sort > {kff_compacted}_sorted.txt"))
        stream = os.popen(f"diff {kff_raw}_sorted.txt {kff_compacted}_sorted.txt")
        stream_val = stream.read()
        stream.close()
        self.assertEqual(stream_val, "")


        print("  clean the test area")
        os.system(f"rm -r {txt} {kff_raw}* {kff_bucket}* {kff_compacted}* {kff_disjoin}")


if __name__ == '__main__':
  unittest.main()
