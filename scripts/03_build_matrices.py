#!/usr/bin/env python3
"""
Build expression matrix + metadata from featureCounts output.
Parses NOTCH3, FLT1, and FN1 transcript variants.

Usage:
    python 03_build_matrices.py

Outputs (to /sc/arion/work/kellej10/cadasil_public_ipsc/):
    - gene_expression_matrix.csv      (genes x samples, raw counts)
    - transcript_variant_counts.csv   (FLT1/FN1/NOTCH3 transcripts x samples)
    - sample_metadata.csv             (sample annotations)
    - notch3_variant_summary.csv      (variant allele evidence per sample)
"""

import os
import csv
import re
from collections import defaultdict

# Paths
SCRATCH = "/sc/arion/scratch/kellej10/cadasil_public_ipsc"
COUNTS_DIR = os.path.join(SCRATCH, "counts")
REF_DIR = os.path.join(SCRATCH, "reference")
OUTPUT_DIR = "/sc/arion/work/kellej10/cadasil_public_ipsc"

# Sample metadata from E-MTAB-16303 SDRF
SAMPLE_META = {
    "WZ12_CADASIL_Control": {
        "err_id": "ERR15965983",
        "genotype": "wild_type",
        "disease": "normal",
        "notch3_mutation": "none",
        "sample_type": "control",
    },
    "WZ4_CADASIL_Control": {
        "err_id": "ERR15965985",
        "genotype": "wild_type",
        "disease": "normal",
        "notch3_mutation": "none",
        "sample_type": "control",
    },
    "WZ8_CADASIL_Control": {
        "err_id": "ERR15965986",
        "genotype": "wild_type",
        "disease": "normal",
        "notch3_mutation": "none",
        "sample_type": "control",
    },
    "WZ16_CADASIL_C224Y": {
        "err_id": "ERR15965984",
        "genotype": "NOTCH3_p.C224Y",
        "disease": "CADASIL",
        "notch3_mutation": "p.C224Y",
        "sample_type": "CADASIL",
    },
    "XZ5_CADASIL_R153C": {
        "err_id": "ERR15965987",
        "genotype": "NOTCH3_p.R153C",
        "disease": "CADASIL",
        "notch3_mutation": "p.R153C",
        "notch3_mutation_observed": "p.C224Y",  # RNA-seq pileup shows C224Y (VAF=0.42), NOT R153C
        "sample_type": "CADASIL",
    },
    "XZ6_CADASIL_R153C_isoCtrl": {
        "err_id": "ERR15965988",
        "genotype": "mutation_corrected",
        "disease": "normal",
        "notch3_mutation": "corrected_R153C",
        "notch3_mutation_observed": "none",  # confirmed corrected
        "sample_type": "isogenic_control",
    },
    "XZ7_CADASIL_C224Y_isoCtrl": {
        "err_id": "ERR15965989",
        "genotype": "mutation_corrected",
        "disease": "normal",
        "notch3_mutation": "corrected_C224Y",
        "sample_type": "isogenic_control",
    },
}

# FLT1 transcript isoforms (Gencode v44, verified from GTF)
# mFlt1 = 30 exons (full-length, membrane-bound receptor with kinase domain)
# sFlt1 = 13-15 exons (soluble, lacks transmembrane + kinase domains)
FLT1_TRANSCRIPTS = {
    "ENST00000282397": "FLT1-201_mFlt1_full_length",    # 30 exons, 7123bp — membrane-bound, MANE_Select
    "ENST00000615840": "FLT1-207_sFlt1_i13",            # 13 exons, 6502bp — soluble, intronic polyA intron 13
    "ENST00000639477": "FLT1-209_sFlt1_i14",            # 14 exons, 6337bp — soluble, intronic polyA intron 14
    "ENST00000541932": "FLT1-204_sFlt1_e15a",           # 15 exons, 2969bp — soluble, alternative exon 15a
    "ENST00000615840": "FLT1-207_sFlt1_i13",            # 13 exons — key soluble isoform
    "ENST00000617835": "FLT1-208_short",                 # 3 exons, 575bp — minor
    "ENST00000543394": "FLT1-205_partial",               # 10 exons, 1275bp
    "ENST00000540678": "FLT1-203_partial",               # 17 exons, 1927bp
    "ENST00000706527": "FLT1-210_partial",               # 16 exons, 1865bp
    "ENST00000539099": "FLT1-202_NMD",                   # 14 exons, 1763bp — NMD
}

# FN1 transcript isoforms (Gencode v44, verified from GTF)
# EDA/EDB inclusion tracked by exon count: 47=EDA+EDB+, 46=one, 45=other, 44=neither
FN1_TRANSCRIPTS = {
    "ENST00000323926": "FN1-201_EDA_EDB_plus",           # 47 exons, 8708bp — both extra domains
    "ENST00000354785": "FN1-203_canonical_MANE",         # 46 exons, 8390bp — MANE_Select, appris_principal
    "ENST00000336916": "FN1-202_46exon",                 # 46 exons, 8435bp
    "ENST00000359671": "FN1-206_45exon",                 # 45 exons, 8524bp
    "ENST00000421182": "FN1-207_45exon",                 # 45 exons, 8103bp
    "ENST00000357867": "FN1-205_EDA_EDB_minus",          # 44 exons, 7898bp — neither extra domain
    "ENST00000356005": "FN1-204_44exon",                 # 44 exons, 7846bp
    "ENST00000446046": "FN1-212_46exon",                 # 46 exons, 7952bp
    "ENST00000443816": "FN1-211_45exon",                 # 45 exons, 7762bp
    "ENST00000432072": "FN1-209_45exon",                 # 45 exons, 7759bp
    "ENST00000426059": "FN1-208_short",                  # 13 exons, 2388bp
}

# NOTCH3 transcripts (Gencode v44, verified from GTF)
NOTCH3_TRANSCRIPTS = {
    "ENST00000263388": "NOTCH3-201_full_length_MANE",    # 33 exons, 8680bp — canonical, MANE_Select
    "ENST00000601011": "NOTCH3-206_truncated",           # 23 exons, 3858bp
    "ENST00000597756": "NOTCH3-204_truncated",           # 3 exons, 559bp
    "ENST00000595514": "NOTCH3-203_NMD",                 # 5 exons, 564bp — NMD
}

# NOTCH3 variant positions (GRCh38, from Ensembl VEP)
# NOTCH3 is on the MINUS strand of chr19, so pileup ref/alt are reverse-complemented:
#   cDNA C>T → genomic G>A (minus strand)
#   cDNA G>A → genomic C>T (minus strand)
NOTCH3_VARIANTS = {
    "p.R153C": {"chrom": "chr19", "pos": 15192182, "ref": "G", "alt": "A", "cdna": "c.457C>T"},
    "p.C224Y": {"chrom": "chr19", "pos": 15191968, "ref": "C", "alt": "T", "cdna": "c.671G>A"},
}

# All transcripts of interest (gene-level tracking)
GENES_OF_INTEREST = ["FLT1", "FN1", "NOTCH3", "VEGFA", "VEGFB", "KDR", "ACTA2",
                     "MYH11", "CNN1", "TAGLN", "MYOCD"]


def parse_featurecounts(filepath):
    """Parse featureCounts output file into a dict of {feature: {sample: count}}."""
    counts = {}
    samples = []
    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if parts[0] == "Geneid":
                # Header row — columns 6+ are sample BAM paths
                samples = [os.path.basename(p).replace("_Aligned.sortedByCoord.out.bam", "")
                           for p in parts[6:]]
                continue
            feature_id = parts[0]
            count_values = [int(x) for x in parts[6:]]
            counts[feature_id] = dict(zip(samples, count_values))
    return counts, samples


def strip_version(ensembl_id):
    """Remove version suffix from Ensembl ID (e.g., ENSG00000001.5 -> ENSG00000001)."""
    return ensembl_id.split(".")[0]


def build_gene_matrix():
    """Build gene-level expression matrix from featureCounts output."""
    gene_counts_file = os.path.join(COUNTS_DIR, "gene_counts.txt")
    if not os.path.exists(gene_counts_file):
        print(f"WARNING: {gene_counts_file} not found. Run pipeline first.")
        return None, None

    counts, samples = parse_featurecounts(gene_counts_file)
    print(f"Gene-level matrix: {len(counts)} genes x {len(samples)} samples")
    return counts, samples


def build_transcript_variants(samples):
    """Extract FLT1, FN1, NOTCH3 transcript-level counts."""
    tx_counts_file = os.path.join(COUNTS_DIR, "transcript_counts.txt")
    if not os.path.exists(tx_counts_file):
        print(f"WARNING: {tx_counts_file} not found. Skipping transcript variants.")
        return None

    counts, tx_samples = parse_featurecounts(tx_counts_file)

    # Collect transcripts of interest
    all_tx = {}
    all_tx.update(FLT1_TRANSCRIPTS)
    all_tx.update(FN1_TRANSCRIPTS)
    all_tx.update(NOTCH3_TRANSCRIPTS)

    results = {}
    for tx_id_versioned, tx_counts in counts.items():
        tx_id = strip_version(tx_id_versioned)
        if tx_id in all_tx:
            label = all_tx[tx_id]
            results[f"{tx_id}|{label}"] = tx_counts

    # Also search by gene name in the featureCounts summary
    # (featureCounts gene_id column may have gene names or Ensembl IDs)

    if results:
        print(f"Transcript variants found: {len(results)}")
    else:
        print("WARNING: No transcript variants matched. IDs may need version suffixes.")
        # Fallback: try matching with versions
        for tx_id_versioned, tx_counts in counts.items():
            tx_base = strip_version(tx_id_versioned)
            if tx_base in all_tx:
                label = all_tx[tx_base]
                results[f"{tx_id_versioned}|{label}"] = tx_counts

    return results


def parse_notch3_pileup():
    """Parse NOTCH3 pileup files for targeted variant allele counts.

    Looks for exact variant positions and counts ref/alt alleles from
    mpileup output, mirroring the AVM project's targeted pileup approach.
    """
    summaries = {}
    for sample in SAMPLE_META:
        pileup_file = os.path.join(COUNTS_DIR, f"{sample}_NOTCH3_pileup.txt")
        if not os.path.exists(pileup_file):
            continue

        # Index pileup by position for fast lookup
        pileup_by_pos = {}
        with open(pileup_file) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 5:
                    chrom, pos, ref_base, depth, bases = parts[0], int(parts[1]), parts[2], int(parts[3]), parts[4]
                    pileup_by_pos[pos] = {
                        "chrom": chrom, "ref_base": ref_base,
                        "depth": depth, "bases": bases,
                    }

        sample_results = {
            "expected_mutation": SAMPLE_META[sample]["notch3_mutation"],
        }

        for var_name, var_info in NOTCH3_VARIANTS.items():
            pos = var_info["pos"]
            alt_base = var_info["alt"]
            ref_base = var_info["ref"]

            if pos in pileup_by_pos:
                p = pileup_by_pos[pos]
                depth = p["depth"]
                bases = p["bases"].upper()
                # Count matches (. and , = ref) and specific alt alleles
                ref_count = bases.count(".") + bases.count(",")
                alt_count = bases.count(alt_base)
                vaf = alt_count / depth if depth > 0 else 0.0
                sample_results[f"{var_name}_depth"] = depth
                sample_results[f"{var_name}_ref_count"] = ref_count
                sample_results[f"{var_name}_alt_count"] = alt_count
                sample_results[f"{var_name}_vaf"] = round(vaf, 4)
            else:
                sample_results[f"{var_name}_depth"] = 0
                sample_results[f"{var_name}_ref_count"] = 0
                sample_results[f"{var_name}_alt_count"] = 0
                sample_results[f"{var_name}_vaf"] = 0.0

        summaries[sample] = sample_results
    return summaries


def write_csv(filepath, header, rows):
    """Write a CSV file."""
    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"  Written: {filepath}")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 1. Sample metadata
    print("\n=== Building sample metadata ===")
    meta_header = ["sample_name", "err_id", "genotype", "disease",
                   "notch3_mutation_reported", "notch3_mutation_observed",
                   "sample_type", "organism", "cell_type",
                   "organism_part", "platform"]
    meta_rows = []
    for sample, meta in sorted(SAMPLE_META.items()):
        meta_rows.append([
            sample, meta["err_id"], meta["genotype"], meta["disease"],
            meta["notch3_mutation"],
            meta.get("notch3_mutation_observed", meta["notch3_mutation"]),
            meta["sample_type"],
            "Homo sapiens",
            "iPSC_derived_neural_crest_VSMC",
            "skin",
            "Illumina_HiSeq4000_PE76",
        ])
    write_csv(os.path.join(OUTPUT_DIR, "sample_metadata.csv"), meta_header, meta_rows)

    # 2. Gene expression matrix
    print("\n=== Building gene expression matrix ===")
    gene_counts, samples = build_gene_matrix()
    if gene_counts and samples:
        sorted_samples = sorted(samples)
        gene_header = ["gene_id"] + sorted_samples
        gene_rows = []
        for gene_id in sorted(gene_counts.keys()):
            row = [gene_id] + [gene_counts[gene_id].get(s, 0) for s in sorted_samples]
            gene_rows.append(row)
        write_csv(os.path.join(OUTPUT_DIR, "gene_expression_matrix.csv"),
                  gene_header, gene_rows)

        # Also write a filtered version with genes of interest
        # featureCounts uses Ensembl gene IDs — map to gene names
        print("\n=== Extracting genes of interest ===")
        # Ensembl gene IDs for genes of interest (Gencode v44 / GRCh38)
        GOI_ENSEMBL = {
            "ENSG00000102755": "FLT1",
            "ENSG00000115414": "FN1",
            "ENSG00000074181": "NOTCH3",
            "ENSG00000112715": "VEGFA",
            "ENSG00000173511": "VEGFB",
            "ENSG00000128052": "KDR",
            "ENSG00000107796": "ACTA2",
            "ENSG00000133392": "MYH11",
            "ENSG00000130598": "CNN1",
            "ENSG00000163017": "TAGLN",
            "ENSG00000141068": "MYOCD",
            "ENSG00000069869": "NEDD4",
            "ENSG00000196136": "SERPINA3",
            "ENSG00000163631": "ALB",
            "ENSG00000170345": "FOS",
            "ENSG00000177606": "JUN",
        }
        goi_rows = []
        for gene_id, counts_dict in gene_counts.items():
            gene_base = strip_version(gene_id)
            if gene_base in GOI_ENSEMBL:
                gene_name = GOI_ENSEMBL[gene_base]
                row = [f"{gene_id}|{gene_name}"] + [counts_dict.get(s, 0) for s in sorted_samples]
                goi_rows.append(row)

        if goi_rows:
            write_csv(os.path.join(OUTPUT_DIR, "genes_of_interest_counts.csv"),
                      gene_header, goi_rows)
            for row in goi_rows:
                print(f"  {row[0]:<45} {row[1:]}")
        else:
            print("  No genes of interest matched.")

    # 3. Transcript variant counts
    print("\n=== Building transcript variant matrix ===")
    if samples:
        tx_variants = build_transcript_variants(samples)
        if tx_variants:
            sorted_samples = sorted(samples)
            tx_header = ["transcript_id|annotation"] + sorted_samples
            tx_rows = []
            for tx_label in sorted(tx_variants.keys()):
                row = [tx_label] + [tx_variants[tx_label].get(s, 0) for s in sorted_samples]
                tx_rows.append(row)
            write_csv(os.path.join(OUTPUT_DIR, "transcript_variant_counts.csv"),
                      tx_header, tx_rows)

    # 4. NOTCH3 variant summary — targeted allele counts at each mutation site
    print("\n=== NOTCH3 variant pileup summary ===")
    notch3_summary = parse_notch3_pileup()
    if notch3_summary:
        var_names = sorted(NOTCH3_VARIANTS.keys())
        n3_header = ["sample_name", "sample_type", "expected_mutation"]
        for vn in var_names:
            vi = NOTCH3_VARIANTS[vn]
            n3_header.extend([
                f"{vn}_{vi['cdna']}_depth",
                f"{vn}_{vi['cdna']}_ref_count",
                f"{vn}_{vi['cdna']}_alt_count",
                f"{vn}_{vi['cdna']}_VAF",
            ])
        n3_rows = []
        for sample in sorted(notch3_summary.keys()):
            s = notch3_summary[sample]
            row = [sample, SAMPLE_META[sample]["sample_type"], s["expected_mutation"]]
            for vn in var_names:
                row.extend([
                    s.get(f"{vn}_depth", 0),
                    s.get(f"{vn}_ref_count", 0),
                    s.get(f"{vn}_alt_count", 0),
                    s.get(f"{vn}_vaf", 0.0),
                ])
            n3_rows.append(row)
        write_csv(os.path.join(OUTPUT_DIR, "notch3_variant_summary.csv"),
                  n3_header, n3_rows)

        # Print a quick summary table to stdout
        print("\n  NOTCH3 Variant Detection Results:")
        print(f"  {'Sample':<35} {'Type':<18} {'Expected':<15} {'R153C VAF':>10} {'C224Y VAF':>10}")
        print("  " + "-" * 90)
        for sample in sorted(notch3_summary.keys()):
            s = notch3_summary[sample]
            print(f"  {sample:<35} {SAMPLE_META[sample]['sample_type']:<18} "
                  f"{s['expected_mutation']:<15} "
                  f"{s.get('p.R153C_vaf', 0.0):>10.4f} "
                  f"{s.get('p.C224Y_vaf', 0.0):>10.4f}")

    print("\n=== All outputs written to", OUTPUT_DIR, "===")
    print("Files:")
    for f in os.listdir(OUTPUT_DIR):
        if f.endswith(".csv"):
            fpath = os.path.join(OUTPUT_DIR, f)
            size = os.path.getsize(fpath)
            print(f"  {f} ({size:,} bytes)")


if __name__ == "__main__":
    main()
