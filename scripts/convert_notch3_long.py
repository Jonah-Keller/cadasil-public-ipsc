#!/usr/bin/env python3
"""Emit long-format notch3_variant_summary.csv for step 09 consumption.

Bridges the legacy wide CSV (one row per sample, columns per variant) to the
iEC-style long format (one row per sample x variant, with chrom/pos/ref/alt
for cluster-aware joining). Run once after 03_build_matrices.py.

Columns: sample,variant,cdna,expected_donor,chrom,pos,ref_expected,alt_expected,depth,ref_count,alt_count,vaf
"""
import csv
from pathlib import Path

PROJ = Path("/sc/arion/work/kellej10/cadasil_public_ipsc")
WIDE = PROJ / "notch3_variant_summary.csv"
OUT = PROJ / "results/variants/notch3_variant_summary.csv"
OUT.parent.mkdir(parents=True, exist_ok=True)

# Authoritative variant coords (NOTCH3 is minus-strand; pileup reports + strand)
VARIANTS = {
    "p.R153C": {"chrom": "chr19", "pos": 15192182, "ref": "G", "alt": "A",
                "cdna": "c.457C>T", "expected_donor": ""},
    "p.C224Y": {"chrom": "chr19", "pos": 15191968, "ref": "C", "alt": "T",
                "cdna": "c.671G>A", "expected_donor": ""},
}

with open(WIDE) as f, open(OUT, "w", newline="") as g:
    r = csv.DictReader(f)
    w = csv.writer(g)
    w.writerow(["sample", "variant", "cdna", "expected_donor", "chrom", "pos",
                "ref_expected", "alt_expected", "depth", "ref_count",
                "alt_count", "vaf"])
    for row in r:
        s = row["sample_name"]
        for var, vi in VARIANTS.items():
            cd = vi["cdna"]
            depth = row.get(f"{var}_{cd}_depth", "0")
            ref = row.get(f"{var}_{cd}_ref_count", "0")
            alt = row.get(f"{var}_{cd}_alt_count", "0")
            vaf = row.get(f"{var}_{cd}_VAF", "0")
            w.writerow([s, var, cd, vi["expected_donor"], vi["chrom"],
                        vi["pos"], vi["ref"], vi["alt"],
                        depth, ref, alt, vaf])

print(f"wrote {OUT}")
