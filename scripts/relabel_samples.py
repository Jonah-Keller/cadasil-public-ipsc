#!/usr/bin/env python3
"""Add display_label / reported_status / actual_status columns to
sample_metadata_reassigned.csv so every downstream report uses the same
unambiguous (reported / actual) naming.

Convention:
  Token : WT  (no NOTCH3 mutation, never carried one)
        : CAD-<MUT>   (CADASIL mutant, carries <MUT>)
        : iCAD-<MUT>  (CADASIL-corrected isogenic control, was <MUT>)
  Label : "<sample> (reported/actual)"  e.g. XZ5 (CAD-R153C/CAD-C224Y)
"""
import csv
from pathlib import Path

CSV = Path("/sc/arion/work/kellej10/cadasil_public_ipsc/data/sample_metadata_reassigned.csv")

# ENA-submitted reported status, derived from sample_type + notch3_mutation_reported
REPORTED = {
    "WZ12_CADASIL_Control":      "WT",
    "WZ16_CADASIL_C224Y":        "CAD-C224Y",
    "WZ4_CADASIL_Control":       "WT",
    "WZ8_CADASIL_Control":       "WT",
    "XZ5_CADASIL_R153C":         "CAD-R153C",
    "XZ6_CADASIL_R153C_isoCtrl": "iCAD-R153C",
    "XZ7_CADASIL_C224Y_isoCtrl": "iCAD-C224Y",
}

# Empirical actual status from NGSCheckMate cluster + NOTCH3 pileup
ACTUAL = {
    "WZ12_CADASIL_Control":      "WT",
    "WZ16_CADASIL_C224Y":        "CAD-C224Y",
    "WZ4_CADASIL_Control":       "WT",
    "WZ8_CADASIL_Control":       "WT",
    "XZ5_CADASIL_R153C":         "CAD-C224Y",   # donor swap: clusters with C224Y
    "XZ6_CADASIL_R153C_isoCtrl": "iCAD-C224Y",  # corrected, but of C224Y donor
    "XZ7_CADASIL_C224Y_isoCtrl": "WT",          # not from C224Y donor at all
}

SHORT = {
    "WZ12_CADASIL_Control":      "WZ12",
    "WZ16_CADASIL_C224Y":        "WZ16",
    "WZ4_CADASIL_Control":       "WZ4",
    "WZ8_CADASIL_Control":       "WZ8",
    "XZ5_CADASIL_R153C":         "XZ5",
    "XZ6_CADASIL_R153C_isoCtrl": "XZ6",
    "XZ7_CADASIL_C224Y_isoCtrl": "XZ7",
}


def label(sid):
    rep, act = REPORTED[sid], ACTUAL[sid]
    return f"{SHORT[sid]} ({rep}/{act})" if rep != act else f"{SHORT[sid]} ({rep})"


def main():
    with open(CSV) as f:
        rows = list(csv.DictReader(f))
    fields = list(rows[0].keys())
    for col in ("short_id", "reported_status", "actual_status", "display_label"):
        if col not in fields:
            fields.append(col)
    for r in rows:
        sid = r["sample_id"]
        r["short_id"]        = SHORT[sid]
        r["reported_status"] = REPORTED[sid]
        r["actual_status"]   = ACTUAL[sid]
        r["display_label"]   = label(sid)
    with open(CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"updated {CSV}")
    for r in rows:
        print(f"  {r['display_label']}")


if __name__ == "__main__":
    main()
