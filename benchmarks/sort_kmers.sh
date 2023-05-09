#!/bin/bash

# Produce a compacted and sorted kff file from any input kff file
# Comand line: bash sort_kmers.sh input.kff m output_sorted.kff output_compacted.kff

set -e

if [ "$#" -eq 5 ]; then
	echo "cmd: bash sort_kmers.sh input.kff m output_sorted.kff output_compacted.kff"
	exit 1
fi

in=$1
m=$2
sort=$3
compact=$4

# 1/4 of the time
kff-tools bucket -i $in -o ${out}_bucket.kff -m $m
# 3/4 of the time
kff-tools compact -i ${out}_bucket.kff -s -o $sort
kff-tools compact -i ${out}_bucket.kff -o $compact

rm -rf ${out}_bucket.kff
