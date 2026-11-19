#!/usr/bin/env bash
# Download the public datasets of the functional tests and derive the text
# inputs of `kff-tools instr`.
#
# usage: tests/functional/prepare_data.sh DATA_DIR
#
# Needs: curl, python3, kmc and kmc_tools (KMC 3) in the PATH.
# Downloads ~420 MB (E. coli, yeast and C. elegans genomes, one E. coli
# Illumina run) and writes ~1.2 GB of text files in DATA_DIR.
set -euo pipefail

if [ $# -ne 1 ]; then
	echo "usage: $0 DATA_DIR" >&2
	exit 1
fi

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
mkdir -p "$1"
cd "$1"

for tool in curl python3 kmc kmc_tools; do
	command -v $tool > /dev/null || { echo "missing: $tool" >&2; exit 1; }
done

get() {
	[ -s "$2" ] && return
	curl -sSfL --retry 3 -o "$2.part" "$1"
	mv "$2.part" "$2"
}

# --- Raw datasets ---
# E. coli K-12 MG1655, S. cerevisiae S288C and C. elegans reference genomes (NCBI)
get https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz ecoli.fna.gz
get https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz yeast.fna.gz
get https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.fna.gz celegans.fna.gz
# E. coli REL606 Illumina reads, SRR2584863 (ENA), first mate only
get https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_1.fastq.gz reads_ecoli_1.fastq.gz

# --- One sequence per line (split on non ACGT characters) ---
F2L="python3 $HERE/fasta2lines.py"
[ -s sars.seq.txt ]     || $F2L "$REPO/data/sars.fa" > sars.seq.txt
[ -s ecoli.seq.txt ]    || $F2L ecoli.fna.gz > ecoli.seq.txt
[ -s ecoli.w80.txt ]    || $F2L ecoli.fna.gz --wrap 80 --min 31 > ecoli.w80.txt
[ -s yeast.seq.txt ]    || $F2L yeast.fna.gz > yeast.seq.txt
[ -s celegans.seq.txt ] || $F2L celegans.fna.gz > celegans.seq.txt
[ -s reads.seq.txt ]    || $F2L reads_ecoli_1.fastq.gz --min 31 > reads.seq.txt

# --- k-mer sets and counts (KMC, canonical k-mers, sorted dumps) ---
mkdir -p kmc_tmp
if [ ! -s reads_k31.counts.txt ]; then
	kmc -k31 -ci2 -cs65535 -t8 -fq reads_ecoli_1.fastq.gz reads_k31 kmc_tmp > /dev/null 2>&1
	kmc_tools transform reads_k31 dump -s reads_k31.dump > /dev/null 2>&1
	tr '\t' ' ' < reads_k31.dump > reads_k31.counts.txt
	rm reads_k31.dump
fi
for k in 21 31 63; do
	if [ ! -s ecoli_k$k.kmers.txt ]; then
		kmc -k$k -ci1 -t8 -fm ecoli.fna.gz ecoli_k$k kmc_tmp > /dev/null 2>&1
		kmc_tools transform ecoli_k$k dump -s ecoli_k$k.dump > /dev/null 2>&1
		cut -d$'\t' -f1 ecoli_k$k.dump > ecoli_k$k.kmers.txt
		rm ecoli_k$k.dump
	fi
done
rm -rf kmc_tmp ./*.kmc_pre ./*.kmc_suf

# The reference hashes are only valid for these exact inputs
sha256sum -c "$HERE/inputs.sha256"
