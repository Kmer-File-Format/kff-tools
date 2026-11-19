# Functional tests on real data

These tests run every `kff-tools` command on real datasets (bacterial,
yeast and nematode genomes, a real Illumina run, k-mer sets and counts for
k = 21, 31 and 63) and record the SHA-256 of every output file and of the
standard output. Comparing these hashes with `reference.tsv` proves that a
modified version still produces **bit-identical** outputs.

`reference.tsv` was produced by kff-tools `83fff4d` (optimisations and
fixes) with kff-cpp-api `b733a09`. Its times and memory peaks were measured
on an Intel Core Ultra 7 165H on AC power, with no other load, and are only
indicative.

`reference_9b90270.tsv` holds the outputs of the version before the
optimisations and fixes (kff-tools `9b90270`, kff-cpp-api `49dd236`), with
their k-mer content. The fixes changed the bytes of 12 outputs (choice of the
bucket of some k-mers, and the outputs computed from these files) but their
k-mer content is the same, except for `compact_s_kmers_k63`: the previous
version lost 82 k-mers and created 82 k-mers that were not in its input.
Every compacted file of the reference has the k-mer content of its input.

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

### Comparing the k-mer content

`diff --mode content` compares the k-mer content of the outputs instead of
their bytes, to validate a change that modifies the layout of the files but
must keep the same k-mers (for example the choice of the bucket of a k-mer):
- kff files and directories of kff files: canonical k-mers and their data
  (decoded with `outstr -c`), sorted, so duplicates count;
- standard output of outstr/shuffle: sorted lines;
- validate: exact text; validate -v: final status line only (its dump
  describes the blocks).

The content digests are computed by `run`; for a run made before (or with
another version of the harness), add them from the outputs left in its work
directory:

```bash
python3 tests/functional/golden.py content --manifest run.tsv --work /tmp/kff-work \
        --decoder bin/kff-tools --out run_with_content.tsv
python3 tests/functional/golden.py diff --mode content tests/functional/reference.tsv run_with_content.tsv
```

`shuffle` is random by design: only its sorted k-mer content is compared.
