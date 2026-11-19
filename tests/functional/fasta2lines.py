#!/usr/bin/env python3
"""Convert a FASTA/FASTQ (optionally gzipped) into one uppercase ACGT sequence per line.
Sequences are split on any non-ACGT character; fragments shorter than --min are dropped."""
import gzip, re, sys, argparse
p = argparse.ArgumentParser()
p.add_argument("input"); p.add_argument("--min", type=int, default=1)
p.add_argument("--wrap", type=int, default=0, help="cut fragments into lines of this size")
p.add_argument("--max-records", type=int, default=0)
a = p.parse_args()
op = gzip.open if a.input.endswith(".gz") else open
split = re.compile(r"[^ACGT]+")
out = sys.stdout
def emit(seq):
    for frag in split.split(seq.upper()):
        if len(frag) < a.min: continue
        if a.wrap:
            for i in range(0, len(frag), a.wrap):
                piece = frag[i:i + a.wrap]
                if len(piece) >= a.min: out.write(piece + "\n")
        else:
            out.write(frag + "\n")
with op(a.input, "rt") as f:
    first = f.read(1); f.seek(0) if not hasattr(f, "peek") else None
with op(a.input, "rt") as f:
    n = 0
    if first == "@":
        for i, line in enumerate(f):
            if i % 4 == 1:
                emit(line.strip()); n += 1
                if a.max_records and n >= a.max_records: break
    else:
        buf = []
        for line in f:
            if line.startswith(">"):
                if buf: emit("".join(buf)); buf = []
            else: buf.append(line.strip())
        if buf: emit("".join(buf))
