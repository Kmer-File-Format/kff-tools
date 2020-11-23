# kff-tools

This repository contains a list of tools that are used to manipulate kff files.
kff file format is [described here](https://github.com/yoann-dufresne/kmer_file_format).

kff-tools is a program containing a set of small programs allowing kff files manipulations.
Each following part describes one of these tools.

## kff-tools split

Split a kff file into one kff file per section.

Parameters:
* **-i <input.kff>** \[required\]: Input file to split.
* **-o <path>**: Directory where the split output files are writen (Default ./).

Usage:
```bash
  kff-tools split -i to_split.kff -o split_dir/
```

Remaining work: The outdir must be created if not exists.

## kff-tools merge

Merge a list of kff files into only one.

code status: TODO

## kff-tools enlarge

Each block in each section is split into one block per kmer.

code status: TODO

## kff-tools outstr

Read a file and print out the kmers and data as strings

code status: TODO

## kff-tools translate

Read and rewrite a kff file changing the nucleotide encoding.

code status: TODO

## kff-tools data-rm

Read a kff file and write the same one with a data size of 0.
i.e. all the data are removed to only keep kmer sequences

code status: TODO

## kff-tools diff

Read two kff files and print the differences between them

code status: TODO
