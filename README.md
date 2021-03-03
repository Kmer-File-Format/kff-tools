# kff-tools

This repository contains a list of tools that are used to manipulate kff files.
kff file format is [described here](https://github.com/yoann-dufresne/kmer_file_format).

kff-tools is a program containing a set of small programs allowing kff files manipulations.
Each following part describes one of these tools.

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

## kff-tools outstr

Read a file and print to stdout the kmers and data as strings (one kmer per line)

Parameters:
* **-i &lt;input.kff&gt;** \[required\]: File to print.

Usage:
```bash
  kff-tools outstr -i to_encode.kff
```

## kff-tools instr

Translate a textual kmer file (tsv format) into one or multiple kff file(s).
Instr supposes that the data are integers.

Parameters:
* **-i &lt;input.kff&gt;** \[required\]: File to translate.
* **-o &lt;outprefix&gt;** \[required\]: prefix of the outputed kff file(s). For unique file the output will be &lt;prefix&gt;.kff. &lt;prefix&gt;\_&lt;minimizer&gt;.kff for multiple files.
* **-s**: Split flag. If set, the kmers will be distributed into one file per minimizer.
* **-m mini_size**: The minimizer size to use if the s flag is set.
* **-d nb_bytes**: Number of Bytes to use to store the counts (Default 0).



Usage:
```bash
  kff-tools instr -i to_encode.txt -o encoded.kff -e AGTC
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

## kff-tools validate

Read a kff file and crash if a bad format is detected.
Print details of the file on verbose mode.

Parameters:
* **-i &lt;input.kff&gt;** \[required\]: File to read.
* **-v** : Verbose mode.

Usage:
```bash
  kff-tools validate -i file.kff -v
```

## kff-tools disjoin

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
