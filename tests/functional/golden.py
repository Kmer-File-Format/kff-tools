#!/usr/bin/env python3
"""Functional harness for kff-tools on real data.

Runs a fixed list of kff-tools commands on real datasets and records, for each
command, the exit code, the wall time, the peak RSS and two digests of every
output:
- bytes:   SHA-256 of the output file (or of stdout);
- content: SHA-256 of its k-mer content, independent of the layout of the file:
           the canonical k-mers and their data (decoded with `outstr -c`),
           sorted, for kff files and directories of kff files; the sorted
           lines for the standard output of outstr/shuffle; the exact text for
           validate; only the final status line for validate -v (its dump
           describes the blocks).

Two runs are compared bit for bit (default) or on the k-mer content.

usage (see tests/functional/README.md):
  golden.py run     --bin BIN --data DATA --work WORK --out run.tsv [--only REGEX] [--set quick|full]
  golden.py content --manifest run.tsv --work WORK --decoder BIN --out run_with_content.tsv
  golden.py diff    [--mode bytes|content] ref.tsv new.tsv
"""
import argparse, hashlib, os, re, shlex, shutil, subprocess, sys, time

# (name, set, command template, outputs)
# {T} = kff-tools binary, {D} = data dir, {W} = work dir.
# Outputs are paths relative to {W}; "stdout" hashes the standard output.
# Commands are run in order: later commands may use earlier outputs.
TESTS = [
    # --- instr: text -> kff ---------------------------------------------------
    ("instr_sars_k21",       "quick", "{T} instr -i {D}/sars.seq.txt -o sars21.kff -k 21 -m 64", ["sars21.kff"]),
    ("instr_genome_k31",     "quick", "{T} instr -i {D}/ecoli.seq.txt -o G31.kff -k 31", ["G31.kff"]),
    ("instr_genome_w80_k31", "quick", "{T} instr -i {D}/ecoli.w80.txt -o W31.kff -k 31 -m 50", ["W31.kff"]),
    ("instr_genome_k63",     "quick", "{T} instr -i {D}/ecoli.seq.txt -o G63.kff -k 63 -m 128", ["G63.kff"]),
    ("instr_kmers_k21",      "quick", "{T} instr -i {D}/ecoli_k21.kmers.txt -o K21.kff -k 21 -m 1", ["K21.kff"]),
    ("instr_kmers_k31",      "quick", "{T} instr -i {D}/ecoli_k31.kmers.txt -o K31.kff -k 31 -m 1", ["K31.kff"]),
    ("instr_kmers_k63",      "quick", "{T} instr -i {D}/ecoli_k63.kmers.txt -o K63.kff -k 63 -m 1", ["K63.kff"]),
    ("instr_counts_k31",     "quick", "{T} instr -i {D}/reads_k31.counts.txt -o C31.kff -k 31 -d 2 -m 1", ["C31.kff"]),
    ("instr_reads_k31",      "full",  "{T} instr -i {D}/reads.seq.txt -o R31.kff -k 31 -m 200", ["R31.kff"]),
    ("instr_yeast_k31",      "full",  "{T} instr -i {D}/yeast.seq.txt -o Y31.kff -k 31", ["Y31.kff"]),
    ("instr_celegans_k31",   "full",  "{T} instr -i {D}/celegans.seq.txt -o E31.kff -k 31", ["E31.kff"]),
    # --- outstr: kff -> text --------------------------------------------------
    ("outstr_sars",          "quick", "{T} outstr -i sars21.kff", ["stdout"]),
    ("outstr_genome_k31",    "quick", "{T} outstr -i G31.kff", ["stdout"]),
    ("outstr_genome_k63_rc", "quick", "{T} outstr -c -i G63.kff", ["stdout"]),
    ("outstr_counts_k31",    "quick", "{T} outstr -i C31.kff", ["stdout"]),
    ("outstr_counts_k31_rc", "quick", "{T} outstr -c -i C31.kff", ["stdout"]),
    ("outstr_reads_k31",     "full",  "{T} outstr -i R31.kff", ["stdout"]),
    # --- validate -------------------------------------------------------------
    ("validate_genome_k31",  "quick", "{T} validate -i G31.kff", ["stdout"]),
    ("validate_v_genome_k31","quick", "{T} validate -v -i G31.kff", ["stdout"]),
    ("validate_v_counts_k31","quick", "{T} validate -v -i C31.kff", ["stdout"]),
    # --- bucket: r -> m sections ----------------------------------------------
    ("bucket_sars_m7",       "quick", "{T} bucket -i sars21.kff -o Bsars.kff -m 7", ["Bsars.kff"]),
    ("bucket_genome_k31_m10","quick", "{T} bucket -i G31.kff -o BG31.kff -m 10", ["BG31.kff"]),
    ("bucket_kmers_k31_m10", "quick", "{T} bucket -i K31.kff -o BK31.kff -m 10", ["BK31.kff"]),
    ("bucket_kmers_k21_m8",  "quick", "{T} bucket -i K21.kff -o BK21.kff -m 8", ["BK21.kff"]),
    ("bucket_kmers_k63_m13", "quick", "{T} bucket -i K63.kff -o BK63.kff -m 13", ["BK63.kff"]),
    ("bucket_counts_k31_m9", "quick", "{T} bucket -i C31.kff -o BC31.kff -m 9", ["BC31.kff"]),
    ("bucket_reads_k31_m11", "full",  "{T} bucket -i R31.kff -o BR31.kff -m 11", ["BR31.kff"]),
    ("validate_v_bucket",    "quick", "{T} validate -v -i Bsars.kff", ["stdout"]),
    ("validate_bucket_k31",  "quick", "{T} validate -i BC31.kff", ["stdout"]),
    ("outstr_bucket_k31",    "quick", "{T} outstr -i BC31.kff", ["stdout"]),
    # --- compact: m sections -> super-kmers -----------------------------------
    ("compact_sars",         "quick", "{T} compact -i Bsars.kff -o CBsars.kff", ["CBsars.kff"]),
    ("compact_s_sars",       "quick", "{T} compact -s -i Bsars.kff -o CSBsars.kff", ["CSBsars.kff"]),
    ("compact_kmers_k31",    "quick", "{T} compact -i BK31.kff -o CBK31.kff", ["CBK31.kff"]),
    ("compact_s_kmers_k31",  "quick", "{T} compact -s -i BK31.kff -o CSBK31.kff", ["CSBK31.kff"]),
    ("compact_kmers_k21",    "quick", "{T} compact -i BK21.kff -o CBK21.kff", ["CBK21.kff"]),
    ("compact_s_kmers_k63",  "quick", "{T} compact -s -i BK63.kff -o CSBK63.kff", ["CSBK63.kff"]),
    ("compact_counts_k31",   "quick", "{T} compact -i BC31.kff -o CBC31.kff", ["CBC31.kff"]),
    ("compact_s_counts_k31", "quick", "{T} compact -s -i BC31.kff -o CSBC31.kff", ["CSBC31.kff"]),
    ("outstr_compact_k31",   "quick", "{T} outstr -i CSBC31.kff", ["stdout"]),
    # --- disjoin --------------------------------------------------------------
    ("disjoin_genome_k31",   "quick", "{T} disjoin -i G31.kff -o DG31.kff", ["DG31.kff"]),
    ("disjoin_bucket_k31",   "quick", "{T} disjoin -i BG31.kff -o DBG31.kff", ["DBG31.kff"]),
    ("disjoin_compact_k31",  "quick", "{T} disjoin -i CBC31.kff -o DCBC31.kff", ["DCBC31.kff"]),
    # --- data-rm --------------------------------------------------------------
    ("datarm_counts_k31",    "quick", "{T} data-rm -i C31.kff -o RC31.kff", ["RC31.kff"]),
    ("datarm_bucket_k31",    "quick", "{T} data-rm -i BC31.kff -o RBC31.kff", ["RBC31.kff"]),
    # --- translate ------------------------------------------------------------
    ("translate_genome_k31", "quick", "{T} translate -e AGTC -i G31.kff -o TG31.kff", ["TG31.kff"]),
    ("translate_bucket_k31", "quick", "{T} translate -e TGCA -i BC31.kff -o TBC31.kff", ["TBC31.kff"]),
    ("outstr_translated",    "quick", "{T} outstr -i TBC31.kff", ["stdout"]),
    # --- index ----------------------------------------------------------------
    ("index_genome_k31",     "quick", "{T} index -i G31.kff -o IG31.kff", ["IG31.kff"]),
    ("index_bucket_k31",     "quick", "{T} index -i BK31.kff -o IBK31.kff", ["IBK31.kff"]),
    # --- merge / split --------------------------------------------------------
    ("merge_mixed",          "quick", "{T} merge -i G31.kff C31.kff K21.kff BC31.kff -o M.kff", ["M.kff"]),
    ("outstr_merge",         "quick", "{T} outstr -i M.kff", ["stdout"]),
    ("split_mixed",          "quick", "mkdir -p split_M && {T} split -i M.kff -o split_M", ["split_M"]),
    ("split_bucket_sars",    "quick", "mkdir -p split_B && {T} split -i Bsars.kff -o split_B", ["split_B"]),
    # --- sort / shuffle -------------------------------------------------------
    ("sort_counts_k31",      "quick", "{T} sort -i C31.kff -o SC31.kff", ["SC31.kff"]),
    ("sort_genome_k31",      "quick", "{T} sort -i W31.kff -o SW31.kff", ["SW31.kff"]),
    # shuffle is random by design: only the sorted k-mer content is compared.
    ("shuffle_counts_k31",   "quick", "{T} shuffle -i C31.kff -o XC31.kff && {T} outstr -i XC31.kff | LC_ALL=C sort", ["stdout"]),
    # --- larger inputs (performance) ------------------------------------------
    ("outstr_celegans",      "full",  "{T} outstr -i E31.kff", ["stdout"]),
    ("bucket_celegans_m11",  "full",  "{T} bucket -i E31.kff -o BE31.kff -m 11", ["BE31.kff"]),
    ("compact_celegans",     "full",  "{T} compact -i BE31.kff -o CBE31.kff", ["CBE31.kff"]),
    ("compact_s_celegans",   "full",  "{T} compact -s -i BE31.kff -o CSBE31.kff", ["CSBE31.kff"]),
    ("disjoin_celegans",     "full",  "{T} disjoin -i E31.kff -o DE31.kff", ["DE31.kff"]),
    ("translate_celegans",   "full",  "{T} translate -e AGTC -i E31.kff -o TE31.kff", ["TE31.kff"]),
    ("validate_celegans",    "full",  "{T} validate -i BE31.kff", ["stdout"]),
]


def sha256_file(path, h=None):
    h = h or hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h


def digest(path):
    """Hash a file, or a directory (sorted file names + contents)."""
    if os.path.isdir(path):
        h = hashlib.sha256()
        for name in sorted(os.listdir(path)):
            h.update(name.encode() + b"\0")
            sha256_file(os.path.join(path, name), h)
        return h.hexdigest()
    return sha256_file(path).hexdigest()


def output_kind(name, output, path):
    """How the content of an output is compared"""
    if output == "stdout":
        if name.startswith("validate_v"):
            return "status"
        if name.startswith("validate"):
            return "text"
        return "lines"
    return "kmers"


def content_digest(name, output, path, decoder, work):
    """Digest of the k-mer content of an output (see the module docstring)"""
    if not os.path.exists(path):
        return "MISSING"
    kind = output_kind(name, output, path)
    if kind == "text":
        return sha256_file(path).hexdigest()
    if kind == "status":
        lines = [line for line in open(path, errors="replace").read().splitlines() if line.strip()]
        return hashlib.sha256((lines[-1] if lines else "").encode()).hexdigest()
    tmp = os.path.join(work, ".sort_tmp")
    os.makedirs(tmp, exist_ok=True)
    sort = f"LC_ALL=C sort -S 2G --parallel=8 -T {shlex.quote(tmp)}"
    if kind == "lines":
        cmd = f"{sort} {shlex.quote(path)} | sha256sum"
    else:
        if os.path.isdir(path):
            decode = f'for f in {shlex.quote(path)}/*.kff; do {shlex.quote(decoder)} outstr -c -i "$f"; done'
        else:
            decode = f"{shlex.quote(decoder)} outstr -c -i {shlex.quote(path)}"
        cmd = f"( {decode} ) | {sort} | sha256sum"
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
    return out.split()[0]


def run(args):
    tool = os.path.abspath(os.path.join(args.bin, "kff-tools"))
    data = os.path.abspath(args.data)
    work = os.path.abspath(args.work)
    os.makedirs(work, exist_ok=True)
    only = re.compile(args.only) if args.only else None
    sets = {"quick"} if args.set == "quick" else {"quick", "full"}
    rows = []
    for name, tset, cmd, outputs in TESTS:
        if tset not in sets or (only and not only.search(name)):
            continue
        cmd = cmd.format(T=tool, D=data, W=work)
        for o in outputs:
            p = os.path.join(work, o)
            if o == "stdout" or not os.path.lexists(p):
                continue
            if os.path.isdir(p):
                shutil.rmtree(p)
            else:
                os.remove(p)
        stdout_path = os.path.join(work, name + ".stdout")
        timing = os.path.join(work, name + ".time")
        full = f"/usr/bin/time -f '%e %M' -o {timing} bash -c \"{cmd}\" > {stdout_path} 2> {os.path.join(work, name + '.stderr')}"
        t0 = time.time()
        code = subprocess.call(full, shell=True, cwd=work)
        wall = time.time() - t0
        try:
            secs, rss = open(timing).read().split()[-2:]
        except Exception:
            secs, rss = f"{wall:.2f}", "0"
        hashes, contents = [], []
        decoder = os.path.abspath(args.decoder) if args.decoder else tool
        for o in outputs:
            p = stdout_path if o == "stdout" else os.path.join(work, o)
            hashes.append(f"{o}={digest(p) if os.path.exists(p) else 'MISSING'}")
            contents.append(f"{o}={content_digest(name, o, p, decoder, work)}")
        row = (name, str(code), secs, rss, " ".join(hashes), " ".join(contents))
        rows.append(row)
        print("\t".join(row[:4]) + "\t" + row[4][:90], flush=True)
    write_manifest(args.out, rows)


def write_manifest(path, rows):
    with open(path, "w") as f:
        f.write("name\texit\tseconds\tmax_rss_kb\toutputs\tcontent\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


def content(args):
    """Adds the content digests to a manifest, from the outputs left in its work directory"""
    work = os.path.abspath(args.work)
    decoder = os.path.abspath(args.decoder)
    rows = []
    with open(args.manifest) as f:
        next(f)
        for line in f:
            fields = line.rstrip("\n").split("\t")
            name = fields[0]
            outputs = [o.split("=")[0] for o in fields[4].split()]
            contents = []
            for o in outputs:
                p = os.path.join(work, name + ".stdout") if o == "stdout" else os.path.join(work, o)
                contents.append(f"{o}={content_digest(name, o, p, decoder, work)}")
            rows.append(tuple(fields[:5]) + (" ".join(contents),))
            print(name, flush=True)
    write_manifest(args.out, rows)


def load(path):
    rows = {}
    with open(path) as f:
        next(f)
        for line in f:
            fields = line.rstrip("\n").split("\t")
            name, code, secs, rss, outs = fields[:5]
            rows[name] = (code, float(secs), int(rss), outs, fields[5] if len(fields) > 5 else None)
    return rows


def diff(args):
    ref, new = load(args.ref), load(args.new)
    column = 4 if args.mode == "content" else 3
    bad = 0
    print(f"comparison: {'k-mer content' if args.mode == 'content' else 'bytes'}")
    print(f"{'test':26s} {'ref s':>8s} {'new s':>8s} {'speedup':>8s} {'ref MB':>8s} {'new MB':>8s}  status")
    tot_r = tot_n = 0.0
    for name in ref:
        if name not in new:
            continue
        rc, rs, rr = ref[name][:3]
        nc, ns, nr = new[name][:3]
        if ref[name][column] is None or new[name][column] is None:
            sys.exit(f"no content digests in {args.ref if ref[name][column] is None else args.new}: "
                     f"add them with the content command")
        same = (rc == nc and ref[name][column] == new[name][column])
        bad += not same
        tot_r += rs; tot_n += ns
        sp = rs / ns if ns > 0 else float("inf")
        print(f"{name:26s} {rs:8.2f} {ns:8.2f} {sp:7.2f}x {rr/1024:8.1f} {nr/1024:8.1f}  {'OK' if same else 'DIFFERENT'}")
    print(f"{'TOTAL':26s} {tot_r:8.2f} {tot_n:8.2f} {tot_r/max(tot_n,1e-9):7.2f}x")
    print("ALL IDENTICAL" if bad == 0 else f"{bad} DIFFERENCES")
    return 1 if bad else 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--bin", required=True); r.add_argument("--data", required=True)
    r.add_argument("--work", required=True); r.add_argument("--out", required=True)
    r.add_argument("--only"); r.add_argument("--set", default="quick", choices=["quick", "full"])
    r.add_argument("--decoder", help="kff-tools used to decode the outputs (default: the tested one)")
    c = sub.add_parser("content")
    c.add_argument("--manifest", required=True); c.add_argument("--work", required=True)
    c.add_argument("--decoder", required=True); c.add_argument("--out", required=True)
    d = sub.add_parser("diff")
    d.add_argument("--mode", default="bytes", choices=["bytes", "content"])
    d.add_argument("ref"); d.add_argument("new")
    a = ap.parse_args()
    sys.exit({"run": run, "content": content, "diff": diff}[a.cmd](a))
