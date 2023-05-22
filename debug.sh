#!/bin/bash

set -e

cmake . && make -j

name=sars
k=63
m=10
# name=small
# k=5
# m=3

kmc -k$k -fm -ci1 data/${name}.fa data/kmc_${name} /tmp
kmc_dump -ci1 data/kmc_${name} /dev/stdout | cut -f1 -d$'\t' > data/${name}.txt

./bin/kff-tools instr -i data/${name}.txt -o data/${name}.kff -k $k -m 10
./bin/kff-tools outstr -i data/${name}.kff -c | sort > data/${name}_comp.txt
./bin/kff-tools bucket -i data/${name}.kff -o data/${name}_m$m.kff -m $m
# ./bin/kff-tools validate -v -i data/${name}_m$m.kff 
# valgrind --gen-suppressions=yes 
valgrind ./bin/kff-tools compact -i data/${name}_m$m.kff -o data/${name}_m${m}_sorted.kff -s
./bin/kff-tools validate -i data/${name}_m${m}_sorted.kff
./bin/kff-tools outstr -i data/${name}_m${m}_sorted.kff -c | sort > data/${name}_m${m}_sorted.txt
diff data/${name}_comp.txt data/${name}_m${m}_sorted.txt