#!/usr/bin/env python3

import argparse
import csv
import gzip
import json
import os
import re
from collections import Counter


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize demultiplex assignment and unknown barcode diagnostics"
    )
    parser.add_argument("--assigned-fastqs", nargs="+", required=True)
    parser.add_argument("--unknown-fastq", required=True)
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--out-prefix", default="demux_unknown_barcode_qc")
    parser.add_argument("--top-n", type=int, default=10)
    return parser.parse_args()


def open_text_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def count_fastq_reads(path):
    lines = 0
    with open_text_maybe_gz(path) as handle:
        for _ in handle:
            lines += 1
    return lines // 4


def extract_header_token(header_line):
    # Supports Illumina-like headers and fqtk outputs.
    # Takes the last whitespace token and strips lane fields if present.
    tok = header_line.strip().split()[-1]
    if ":" in tok:
        tok = tok.rsplit(":", 1)[-1]
    return tok


def barcode_from_token(token):
    token = token.strip()
    if token in {"", "0", "N", "NNNN", "unknown", "UNKNOWN"}:
        return "UNKNOWN"
    return token


def read_unknown_barcode_counts(path):
    counts = Counter()
    total = 0
    with open_text_maybe_gz(path) as handle:
        for i, line in enumerate(handle):
            if i % 4 != 0:
                continue
            token = extract_header_token(line)
            barcode = barcode_from_token(token)
            counts[barcode] += 1
            total += 1
    return counts, total


def parse_samplesheet(path):
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.read().splitlines()

    data_start = 0
    for i, line in enumerate(lines):
        if line.strip().lower() == "[data]":
            data_start = i + 1
            break

    data_lines = [line for line in lines[data_start:] if line.strip()]
    reader = csv.DictReader(data_lines)

    known_pairs = set()
    known_tokens = set()
    known_i1 = set()
    known_i2 = set()

    for row in reader:
        i1 = (row.get("index") or row.get("Index") or "").strip()
        i2 = (row.get("index2") or row.get("Index2") or "").strip()
        if i1:
            known_i1.add(i1)
        if i2:
            known_i2.add(i2)
        if i1 and i2:
            known_pairs.add((i1, i2))
            known_tokens.add(f"{i1}+{i2}")

    return known_pairs, known_tokens, known_i1, known_i2


def hamming_if_same_len(a, b):
    if len(a) != len(b):
        return None
    return sum(1 for x, y in zip(a, b) if x != y)


def classify_against_samplesheet(barcode, known_pairs, known_tokens, known_i1, known_i2):
    if barcode == "UNKNOWN":
        return "random"

    if barcode in known_tokens:
        return "exact"

    parts = barcode.split("+")
    if len(parts) >= 2:
        i1 = parts[0]
        i2 = parts[1]
        if (i1, i2) in known_pairs:
            return "exact"

        near_i1 = any((hamming_if_same_len(i1, k) is not None and hamming_if_same_len(i1, k) <= 1) for k in known_i1)
        near_i2 = any((hamming_if_same_len(i2, k) is not None and hamming_if_same_len(i2, k) <= 1) for k in known_i2)
        if near_i1 and near_i2:
            return "near"

    return "random"


def longest_homopolymer_run(seq):
    if not seq:
        return 0
    best = 1
    cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 1
    return best


def low_diversity_flags(barcode):
    if barcode == "UNKNOWN":
        return False, False
    parts = barcode.split("+")
    joined = "".join(parts)
    if not joined:
        return False, False

    run = longest_homopolymer_run(joined)
    base_counts = Counter(joined)
    max_base_frac = max(base_counts.values()) / float(len(joined))

    has_low_div = run >= 8 or max_base_frac >= 0.8
    has_g_homopolymer = "GGGGGGGG" in joined
    return has_low_div, has_g_homopolymer


def summarize_indicators(unknown_pct, exact_pct, near_pct, random_pct, low_div_pct, g_hpoly_pct, top10_pct):
    notes = []

    if unknown_pct > 20.0:
        notes.append("High unknown fraction (>20%): demultiplexing performance is suboptimal.")

    if exact_pct + near_pct >= 30.0:
        notes.append("Unknown barcodes substantially overlap expected barcodes (exact+near >=30%): possible index/sample-sheet mismatch or permissive mismatch behavior.")

    if low_div_pct >= 10.0:
        notes.append("Low-diversity unknown barcode signal is high (>=10%): possible sequencing/index diversity issue.")

    if random_pct >= 90.0 and top10_pct < 25.0:
        notes.append("Unknown barcodes are mostly random and weakly concentrated: pattern is more consistent with background noise/contamination than with a single barcode swap.")

    if random_pct >= 95.0 and low_div_pct < 5.0 and g_hpoly_pct < 2.0:
        notes.append("Very high random unknown fraction with low low-diversity signal: compatible with PhiX-like/background noise patterns.")

    if not notes:
        notes.append("No single dominant failure mode detected by heuristics.")

    return notes


def main():
    args = parse_args()

    assigned_r1 = [
        p for p in args.assigned_fastqs
        if re.search(r"(_R1_.*\.f(ast)?q(\.gz)?$)|(\.R1\.f(ast)?q(\.gz)?$)", os.path.basename(p))
    ]
    assigned_reads = sum(count_fastq_reads(p) for p in assigned_r1)

    unknown_counts, unknown_reads = read_unknown_barcode_counts(args.unknown_fastq)
    total_reads = assigned_reads + unknown_reads

    if total_reads == 0:
        raise SystemExit("No reads found in assigned or unknown FASTQ inputs")

    known_pairs, known_tokens, known_i1, known_i2 = parse_samplesheet(args.samplesheet)

    cls_counter = Counter()
    low_div_reads = 0
    g_hpoly_reads = 0

    for bc, ct in unknown_counts.items():
        cls = classify_against_samplesheet(bc, known_pairs, known_tokens, known_i1, known_i2)
        cls_counter[cls] += ct

        low_div, g_hpoly = low_diversity_flags(bc)
        if low_div:
            low_div_reads += ct
        if g_hpoly:
            g_hpoly_reads += ct

    top = unknown_counts.most_common(args.top_n)
    top10_total = sum(c for _, c in top)

    assigned_pct = (100.0 * assigned_reads) / total_reads
    unknown_pct = (100.0 * unknown_reads) / total_reads

    exact_pct = (100.0 * cls_counter["exact"]) / unknown_reads if unknown_reads else 0.0
    near_pct = (100.0 * cls_counter["near"]) / unknown_reads if unknown_reads else 0.0
    random_pct = (100.0 * cls_counter["random"]) / unknown_reads if unknown_reads else 0.0

    low_div_pct = (100.0 * low_div_reads) / unknown_reads if unknown_reads else 0.0
    g_hpoly_pct = (100.0 * g_hpoly_reads) / unknown_reads if unknown_reads else 0.0
    top10_pct = (100.0 * top10_total) / unknown_reads if unknown_reads else 0.0

    indicators = summarize_indicators(
        unknown_pct=unknown_pct,
        exact_pct=exact_pct,
        near_pct=near_pct,
        random_pct=random_pct,
        low_div_pct=low_div_pct,
        g_hpoly_pct=g_hpoly_pct,
        top10_pct=top10_pct,
    )

    metrics = {
        "assigned_reads": assigned_reads,
        "unknown_reads": unknown_reads,
        "total_reads": total_reads,
        "assigned_pct": assigned_pct,
        "unknown_pct": unknown_pct,
        "unknown_exact_match_pct": exact_pct,
        "unknown_near_match_pct": near_pct,
        "unknown_random_pct": random_pct,
        "unknown_low_diversity_pct": low_div_pct,
        "unknown_g_homopolymer_pct": g_hpoly_pct,
        "unknown_top10_share_pct": top10_pct,
    }

    with open(f"{args.out_prefix}.metrics.tsv", "w", encoding="utf-8") as out:
        out.write("metric\tvalue\n")
        for k, v in metrics.items():
            out.write(f"{k}\t{v}\n")

    with open(f"{args.out_prefix}.top_unknowns.tsv", "w", encoding="utf-8") as out:
        out.write("rank\tbarcode\tcount\tunknown_pct\n")
        for i, (bc, ct) in enumerate(top, start=1):
            pct = (100.0 * ct) / unknown_reads if unknown_reads else 0.0
            out.write(f"{i}\t{bc}\t{ct}\t{pct:.4f}\n")

    with open(f"{args.out_prefix}.indicators.txt", "w", encoding="utf-8") as out:
        out.write("Demux Unknown Barcode QC Summary\n")
        out.write(f"Assigned reads: {assigned_reads} ({assigned_pct:.2f}%)\n")
        out.write(f"Unknown reads: {unknown_reads} ({unknown_pct:.2f}%)\n")
        out.write(f"Unknown exact barcode matches to samplesheet: {exact_pct:.2f}%\n")
        out.write(f"Unknown near barcode matches to samplesheet: {near_pct:.2f}%\n")
        out.write(f"Unknown random/novel barcodes: {random_pct:.2f}%\n")
        out.write(f"Unknown low-diversity barcodes: {low_div_pct:.2f}%\n")
        out.write(f"Unknown G-homopolymer (>=8 G) barcodes: {g_hpoly_pct:.2f}%\n")
        out.write(f"Top {args.top_n} unknown barcode share: {top10_pct:.2f}%\n")
        out.write("\nPotential issue indicators\n")
        for note in indicators:
            out.write(f"- {note}\n")

    with open(f"{args.out_prefix}.metrics.json", "w", encoding="utf-8") as out:
        json.dump(
            {
                "metrics": metrics,
                "top_unknowns": [
                    {
                        "rank": i,
                        "barcode": bc,
                        "count": ct,
                        "unknown_pct": (100.0 * ct) / unknown_reads if unknown_reads else 0.0,
                    }
                    for i, (bc, ct) in enumerate(top, start=1)
                ],
                "indicators": indicators,
            },
            out,
            indent=2,
        )


if __name__ == "__main__":
    main()
