# kff-tools

This repository contains a list of tools that are used to manipulate kff files.
kff file format is [described here](https://github.com/yoann-dufresne/kmer_file_format).

kff-tools is a program containing a set of small programs allowing kff files manipulations.
Each following part describes one of these tools.


## kff-tools instr

Convert a text kmer file or a text sequence file into a kff file. ACTG encoding is used.
Kmer files must contain 1 kmer per line.
Example (k=3):
```
  ACC
  CCT
  GTA
```
Sequence files must contains 1 sequence of compacted kmers per line.
Example (k=3):
```
  ACCT   # contains the 2 kmers ACC and CCT
  TGC
```
The counts flag can be added for kmer files.
It allows counts at the end of each line after any 1 char separator.
If the c flag is activated, do not forget to set the data size sufficiently large.
Example:
```
  ACC 12
  CCT 1
  GTA 42
```

Parameters:
* **-i &lt;input.kff&gt;** \[required\]: File to translate.
* **-o &lt;output.kff&gt;** \[required\]: Output kff file.
* **-k kmer_size** \[required\]: kmer size.

* **-d nb_bytes**: Data size in Bytes (Default 0, max 8).
* **-c**: Read counts at the end of the lines in a kmer txt file.
* **-m max_size**: Set the maximum number of kmer per sequence in the kff file. Set to 0 when 1 kmer per line. Split the sequences that are too long when reading a sequence file.



Usage:
```bash
  # Read 1 kmer per line
  kff-tools instr -i kmers.txt -o kmers.kff -k 12
  # Read 1 kmer per line with its counts (up to 255)
  kff-tools instr -i counts.txt -o counts.kff -k 12 -c -d 1
  # Read sequences and split them if the contains more than 256 kmers
  kff-tools instr -i sequences.txt -o sequences.kff -k 12 -m 256
```

## kff-tools outstr

Read a kff file and print to stdout the kmers and data as strings (one kmer per line)

Parameters:
* **-i &lt;input.kff&gt;** \[required\]: File to print.
* **-c**: Print the encoding lexicagraphic minimal string between each kmer and its reverse complement (can affact the computation time).

Usage:
```bash
  kff-tools outstr -i file.kff
```

## kff-tools validate

Read a kff file and exit raising an error if a file corruption is detected.
Print details of the file on verbose mode.

Parameters:
* **-i &lt;input.kff&gt;** \[required\]: File to read.
* **-v** : Verbose mode.

Usage:
```bash
  kff-tools validate -i file.kff -v
```


## kff-tools split

Split a kff file into one kff file per section.

Parameters:
* **-i &lt;input.kff&gt;** \[required\]: Input file to split.
* **-o &lt;path&gt;**: Directory where the split output files are written (Default ./).

Usage:
```bash
  kff-tools split -i to_split.kff -o split_dir/
```

## kff-tools merge

Merge a list of kff files into only one.
The order of the input file will be preserved in the merged output.

Parameters:
* **-i &lt;input1.kff&gt; &lt;input2.kff&gt; ...** \[required\]: Input file list to merge.
All the files must share the same encoding.
If not, please first translate them (you can use the translate subprogram of kff-tools).
* **-o &lt;output.kff&gt;** \[required\]: Name of the merged kff file.

Usage:
```bash
  kff-tools merge -i to_merge_1.kff to_merge_2.kff to_merge_3.kff -o merged.kff
```


## kff-tools compact

TODO

## kff-tools disjoin

This tool is the opposite of the compact tool.
Each block in each section is split into one block per kmer.
The number of kmers inside of each section is preserved.
The number of blocks per section is increased to 1 block per kmer.

Parameters:
* **-i &lt;input.kff&gt;** \[required\]: File to disjoin.
* **-o &lt;disjoin.kff&gt;** \[required\]: Disjoin output file.

Usage:
```bash
  kff-tools disjoin -i file.kff -o disjoin.kff
```



## kff-tools translate

Read and rewrite a kff file changing the nucleotide encoding.

Parameters:
* **-i &lt;input.kff&gt;** \[required\]: File to translate.
* **-o &lt;output.kff&gt;** \[required\]: Translated file.
* **-e &lt;encoding&gt;** \[required\]: 4 chars encoding. All the letters A, C, G and T must be present in the encoding order.
For example, AGTC represent the encoding A=0, G=1, T=2, C=3.

Usage:
```bash
  kff-tools translate -i to_encode.kff -o encoded.kff -e AGTC
```

## kff-tools data-rm

Read a kff file and write the same one with a data size of 0.
It means that all the data are removed and the file only preserve sequences.

Parameters:
* **-i &lt;input.kff&gt;** \[required\]: File to copy.
* **-o &lt;output.kff&gt;** \[required\]: Copied file without data.

Usage:
```bash
  kff-tools data-rm -i file.kff -o file_nodata.kff
```




# Test the code

Run functional tests from the root of the project

```bash
  python3 -m unittest discover -s tests/
```
