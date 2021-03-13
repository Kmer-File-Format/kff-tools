import unittest
import os


# file = "kmers"
file = "ecoli_count_dsk"


class TestInOut(unittest.TestCase):
  def test_raw_section(self):
    print("\n-- TestInOut test_raw_section")
    # Generate a kff file from a textual kmer count
    print("1/3 - Generate a kff file from a textual kmer count")
    self.assertEqual(0, os.system(f"./kff-tools instr --data-size 2 --infile data/{file}_sorted.txt --outprefix data/{file}_instr"))
    # Regenerate a textual file from the kff
    print("2/3 - Regenerate a txt file from the kff")
    self.assertEqual(0, os.system(f"./kff-tools outstr --infile data/{file}_instr.kff > data/{file}_instr_outstr.txt"))
    self.assertEqual(0, os.system(f"rm data/{file}_instr.kff"))
    # Compate the original file to the final translated file
    print("3/3 - Compare initial and final txt files")
    stream = os.popen(f"diff data/{file}_sorted.txt data/{file}_instr_outstr.txt")
    self.assertEqual(0, os.system(f"rm data/{file}_instr_outstr.txt"))
    self.assertEqual(stream.read(), "")
    stream.close()

  def test_bucketized(self):
    bucket_size = [4]
    for size in bucket_size:
      print(f"\n-- TestInOut test_bucketized_{size}")
      # Generate a kff file from a textual kmer count bucketing the kmers using a minimizer of size 4
      print(f"1/4 - Generate a kff file from a textual kmer count. Bucket size : {size}")
      self.assertEqual(0, os.system(f"./kff-tools instr --data-size 2 --mini-size {size} --infile data/{file}_sorted.txt --outprefix data/{file}_instr_b{size}"))
      # Regenerate a textual file from the kff
      print("2/4 - Regenerate a txt file from the kff")
      self.assertEqual(0, os.system(f"./kff-tools outstr --infile data/{file}_instr_b{size}.kff > data/{file}_instr_b{size}_outstr.txt"))
      self.assertEqual(0, os.system(f"rm data/{file}_instr_b4.kff"))
      # Sort the regenerated file to match the original order
      print("3/4 - Sort final output to match original file")
      self.assertEqual(0, os.system(f"sort data/{file}_instr_b{size}_outstr.txt > data/{file}_instr_b{size}_outstr_sorted.txt"))
      self.assertEqual(0, os.system(f"rm data/{file}_instr_b{size}_outstr.txt"))
      # Compare initial and final files
      print("4/4 - Compare initial and final txt files")
      stream = os.popen(f"diff data/{file}_sorted.txt data/{file}_instr_b{size}_outstr_sorted.txt")
      os.system(f"rm data/{file}_instr_b{size}_outstr_sorted.txt")
      content = stream.read()
      stream.close()
      self.assertEqual(content, "")



class TestSplitMerge(unittest.TestCase):
  def test_basic_split(self):
    print(f"\n-- TestSplitMerge test_basic_split")
    # Prepare a file with sections
    print(f"1/4 - Generate a kff file from a textual kmer count. Bucket size : 4")
    self.assertEqual(0, os.system(f"./kff-tools instr --data-size 2 --mini-size 4 --infile data/{file}_sorted.txt --outprefix data/{file}_instr_b4"))
    # Split the file into one file per section
    self.assertEqual(0, os.system(f"rm -rf data/{file}_instr_b4/ && mkdir data/{file}_instr_b4/"))
    self.assertEqual(0, os.system(f"./kff-tools split --infile data/{file}_instr_b4.kff --outdir data/{file}_instr_b4/"))
    # Merge outstrs
    print(f"3/4 - outstr each file")
    sections = list(os.listdir(f"data/{file}_instr_b4/"))
    self.assertEqual(0, os.system(f"touch data/{file}_instr_b4/concat.txt"))
    for section_file in sections:
      self.assertEqual(0, os.system(f"./kff-tools outstr -i data/{file}_instr_b4/{section_file} >> data/{file}_instr_b4/concat.txt"))
    # Sort
    print(f"4/4 - Sort and compare output")
    self.assertEqual(0, os.system(f"sort data/{file}_instr_b4/concat.txt > data/{file}_instr_b4/concat_sorted.txt"))
    self.assertEqual(0, os.system(f"rm data/{file}_instr_b4/concat.txt"))
    stream = os.popen(f"diff data/{file}_sorted.txt data/{file}_instr_b4/concat_sorted.txt")
    os.system(f"rm -rf data/{file}_instr_b4/")
    self.assertEqual(stream.read(), "")
    stream.close()

  def test_split_merge(self):
    print(f"\n-- TestSplitMerge test_split_merge")
    # Prepare a file with sections
    print(f"1/3 - Generate a kff file from a textual kmer count. Bucket size : 4")
    self.assertEqual(0, os.system(f"./kff-tools instr --data-size 2 --mini-size 4 --infile data/{file}_sorted.txt --outprefix data/{file}_instr_b4"))
    # Split the file into one file per section
    self.assertEqual(0, os.system(f"rm -rf data/{file}_instr_b4/ && mkdir data/{file}_instr_b4/"))
    self.assertEqual(0, os.system(f"./kff-tools split --infile data/{file}_instr_b4.kff --outdir data/{file}_instr_b4/"))
    # Merge outstrs
    print(f"2/3 - Merge splited files")
    sections = list(os.listdir(f"data/{file}_instr_b4/"))
    sections = [f"data/{file}_instr_b4/{section_file}" for section_file in sections]
    self.assertEqual(0, os.system(f"./kff-tools merge -i {' '.join(sections)} -o data/{file}_instr_b4_split_merged.kff"))
    self.assertEqual(0, os.system(f"rm -rf data/{file}_instr_b4/"))
    # Sort
    print(f"3/3 - Sort and compare output")
    self.assertEqual(0, os.system(f"./kff-tools outstr -i data/{file}_instr_b4_split_merged.kff > data/{file}_instr_b4_split_merged.txt"))
    self.assertEqual(0, os.system(f"rm -rf data/{file}_instr_b4_split_merged.kff"))
    self.assertEqual(0, os.system(f"sort data/{file}_instr_b4_split_merged.txt > data/{file}_instr_b4_split_merged_sorted.txt"))
    self.assertEqual(0, os.system(f"rm -rf data/{file}_instr_b4_split_merged.txt"))
    stream = os.popen(f"diff data/{file}_sorted.txt data/{file}_instr_b4_split_merged_sorted.txt")
    os.system(f"rm -rf data/{file}_instr_b4/concat_sorted.txt")
    self.assertEqual(stream.read(), "")
    stream.close()


if __name__ == '__main__':
  unittest.main()
