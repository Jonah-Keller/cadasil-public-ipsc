#!/usr/bin/env python3
# 24_annotate_isoforms.py -- add biological annotation to isoform DE tables.
#
# Input  : results/de_reassigned/isoform/tx_de_<contrast>_panel.csv
#          (produced by 16_isoform_de.R, covers FLT1/FLT4/FN1)
# Output : results/de_reassigned/isoform/{FLT1,FLT4,FN1}_isoform_stats_annotated.csv
#
# Annotations come from _isoform_annotations.ANNOT (single source of truth
# shared with the iEC repo).
import csv
import sys
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _isoform_annotations import ANNOT, strip_version  # noqa: E402

ROOT = Path("/sc/arion/work/kellej10/cadasil_public_ipsc")
ISO_DIR = ROOT / "results/de_reassigned/isoform"

FIELDS = ["transcript_id", "ensembl_name", "isoform_label", "aa_length",
          "n_exons", "biotype", "notes", "contrast",
          "AveExpr", "logFC", "t", "P.Value", "adj.P.Val"]


def annotate(rows_by_gene, gene_label):
    rows = rows_by_gene.get(gene_label, [])
    out_path = ISO_DIR / f"{gene_label}_isoform_stats_annotated.csv"
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        n_missing = 0
        for r in rows:
            tid_full = r["transcript_id"]
            a = ANNOT.get(strip_version(tid_full))
            if a is None:
                n_missing += 1
                a = {"ensembl_name": "", "isoform_label": "(unannotated)",
                     "aa_length": None, "n_exons": None,
                     "biotype": r.get("transcript_type", ""), "notes": ""}
            w.writerow({
                "transcript_id":  tid_full,
                "ensembl_name":   a["ensembl_name"],
                "isoform_label":  a["isoform_label"],
                "aa_length":      a["aa_length"] if a["aa_length"] is not None else "NA",
                "n_exons":        a["n_exons"]   if a["n_exons"]   is not None else "NA",
                "biotype":        a["biotype"],
                "notes":          a["notes"],
                "contrast":       r["contrast"],
                "AveExpr":        r["AveExpr"],
                "logFC":          r["logFC"],
                "t":              r["t"],
                "P.Value":        r["P.Value"],
                "adj.P.Val":      r["adj.P.Val"],
            })
        if n_missing:
            print(f"[warn] {gene_label}: {n_missing} rows had no annotation entry")
    print(f"  wrote {out_path} ({len(rows)} rows)")


def main():
    print("=== annotating isoform DE tables ===")
    rows_by_gene = defaultdict(list)
    for panel_csv in sorted(ISO_DIR.glob("tx_de_*_panel.csv")):
        with open(panel_csv) as f:
            for r in csv.DictReader(f):
                g = r.get("gene_name", "")
                if g in ("FLT1", "FLT4", "FN1"):
                    rows_by_gene[g].append(r)
    for gene in ("FLT1", "FLT4", "FN1"):
        annotate(rows_by_gene, gene)
    print("=== done ===")


if __name__ == "__main__":
    main()
