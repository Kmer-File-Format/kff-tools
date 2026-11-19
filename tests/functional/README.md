# Functional tests on real data

These tests run every `kff-tools` command on real datasets (bacterial,
yeast and nematode genomes, a real Illumina run, k-mer sets and counts for
k = 21, 31 and 63) and record the SHA-256 of every output file and of the
standard output. Comparing these hashes with `reference.tsv` proves that a
modified version still produces **bit-identical** outputs.

`reference.tsv` was produced by kff-tools `9b90270` (main) with kff-cpp-api
`49dd236`. The times and memory peaks it contains were measured on an Intel
Core Ultra 7 165H and are only indicative.

## Prepare the data (once)

Needs `curl`, `python3` and KMC 3 (`kmc`, `kmc_tools`).

```bash
tests/functional/prepare_data.sh /path/to/kff-data   # ~420 MB downloaded
```

## Run and compare

```bash
# quick set (E. coli, ~2 min) or full set (+ reads, yeast, C. elegans)
python3 tests/functional/golden.py run --bin bin --data /path/to/kff-data \
        --work /tmp/kff-work --out new.tsv --set quick
python3 tests/functional/golden.py diff tests/functional/reference.tsv new.tsv
```

`diff` prints the time and memory of each command in both runs and ends with
`ALL IDENTICAL` when every output is bit-identical (exit code 0).

`shuffle` is random by design: only its sorted k-mer content is compared.
