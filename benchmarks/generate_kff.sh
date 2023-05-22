#!/bin/bash

# Generate a kff file from a dataset and a k value
# Script command line: bash generate_kff sequences.fasta k output.kff

set -e

if [ "$#" -eq 4 ]; then
	echo "cmd: bash generate_kff sequences.fasta k output.kff"
	exit 1
fi

fa=$1
k=$2
out=$3

mkdir -p ${out}_tmpdir

kmc -k$k -fm -ci1 $fa ${fa}_kmc ${out}_tmpdir
kmc_dump -ci1 ${fa}_kmc /dev/stdout | cut -f1 -d$'\t' > ${fa}.txt
kff-tools  instr -i ${fa}.txt -o $out -k $k -m 1

rm -rf ${out}_tmpdir ${fa}.txt ${fa}_kmc_pre ${fa}_kmc_suf
