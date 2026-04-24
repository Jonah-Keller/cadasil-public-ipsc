#!/usr/bin/env python3
"""Compare FLT1 and FN1 isoform expression between CADASIL and controls."""

import csv, math, os

EXPORT = "/sc/arion/work/kellej10/cadasil_public_ipsc/export"

def load_tx(fname):
    with open(os.path.join(EXPORT, "variants", fname)) as f:
        lines = list(csv.reader(f))
    return lines[1], lines[2:]  # header, data (skip comment)

def to_cpm(counts, lib_sizes):
    return [(c / ls) * 1e6 if ls > 0 else 0 for c, ls in zip(counts, lib_sizes)]

def mean(vals):
    return sum(vals) / len(vals) if vals else 0

def log2fc(a, b):
    if b < 0.01:
        return float('inf') if a > 0.01 else 0
    return math.log2((a + 0.5) / (b + 0.5))

# Load library sizes
with open(os.path.join(EXPORT, "qc", "library_sizes.csv")) as f:
    reader = csv.DictReader(f)
    lib_sizes = {row['sample']: int(row['filtered_total_counts']) for row in reader}

# Groups (corrected metadata)
cadasil = ['WZ16_CADASIL_C224Y', 'XZ5_CADASIL_R153C']
controls = ['WZ12_CADASIL_Control', 'WZ4_CADASIL_Control',
            'WZ8_CADASIL_Control', 'XZ7_CADASIL_C224Y_isoCtrl']
isogenic = ['XZ6_CADASIL_R153C_isoCtrl']

def analyze_gene(fname, gene_name):
    header, rows = load_tx(fname)
    samples = header[2:]
    lib_arr = [lib_sizes[s] for s in samples]

    cad_idx = [samples.index(s) for s in cadasil]
    ctrl_idx = [samples.index(s) for s in controls]
    iso_idx = [samples.index(s) for s in isogenic]

    print(f"\n{'=' * 105}")
    print(f"{gene_name} ISOFORM ANALYSIS: CADASIL (n=2) vs Controls (n=4) vs Isogenic (n=1)")
    print(f"{'=' * 105}")
    print(f"{'Transcript':<45} {'CADASIL CPM':>12} {'Control CPM':>12} {'Isogenic CPM':>13} {'log2FC':>8}  Direction")
    print("-" * 105)

    results = []
    for row in rows:
        annot = row[1]
        counts = [int(x) for x in row[2:]]
        cpm = to_cpm(counts, lib_arr)

        cad_mean = mean([cpm[i] for i in cad_idx])
        ctrl_mean = mean([cpm[i] for i in ctrl_idx])
        iso_mean = mean([cpm[i] for i in iso_idx])
        fc = log2fc(cad_mean, ctrl_mean)

        direction = ""
        if abs(fc) > 0.5:
            direction = "UP in CADASIL" if fc > 0 else "DOWN in CADASIL"

        print(f"  {annot:<43} {cad_mean:>12.2f} {ctrl_mean:>12.2f} {iso_mean:>13.2f} {fc:>8.2f}  {direction}")
        results.append((annot, cpm, cad_mean, ctrl_mean, iso_mean, fc))

    return samples, lib_arr, results, cad_idx, ctrl_idx, iso_idx


# === FLT1 ===
samples, lib_arr, flt1_results, cad_idx, ctrl_idx, iso_idx = analyze_gene(
    "flt1_isoform_counts.csv", "FLT1")

# sFlt1 / mFlt1 ratio
print(f"\n  {'sFlt1/mFlt1 ratio per sample:'}")
mflt1 = sflt1_total = None
for annot, cpm, _, _, _, _ in flt1_results:
    if "mFlt1_full_length" in annot:
        mflt1 = cpm
    if "sFlt1" in annot:
        if sflt1_total is None:
            sflt1_total = [0] * len(cpm)
        for j, c in enumerate(cpm):
            sflt1_total[j] += c

if mflt1 and sflt1_total:
    for j, s in enumerate(samples):
        ratio = sflt1_total[j] / mflt1[j] if mflt1[j] > 0.01 else float('inf')
        group = "CADASIL" if j in cad_idx else ("ISO" if j in iso_idx else "CTRL")
        print(f"    {s:<40} {group:<8}  sFlt1/mFlt1 = {ratio:.1f}")

# === FN1 ===
samples, lib_arr, fn1_results, cad_idx, ctrl_idx, iso_idx = analyze_gene(
    "fn1_isoform_counts.csv", "FN1")

# EDA+/EDA- ratio
print(f"\n  {'EDA+EDB+ / EDA-EDB- ratio per sample:'}")
eda_plus = eda_minus = None
for annot, cpm, _, _, _, _ in fn1_results:
    if "EDA_EDB_plus" in annot:
        eda_plus = cpm
    elif "EDA_EDB_minus" in annot:
        eda_minus = cpm

if eda_plus and eda_minus:
    for j, s in enumerate(samples):
        ratio = eda_plus[j] / eda_minus[j] if eda_minus[j] > 0.01 else float('inf')
        group = "CADASIL" if j in cad_idx else ("ISO" if j in iso_idx else "CTRL")
        print(f"    {s:<40} {group:<8}  EDA+/EDA- = {ratio:.4f}")

# === NOTCH3 ===
analyze_gene("notch3_isoform_counts.csv", "NOTCH3")

print("\n")
