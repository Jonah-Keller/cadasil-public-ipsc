#!/usr/bin/env python3
"""
Build a clean, QC'd export package for local CADASIL signature analysis.

Standard bulk RNA-seq preprocessing:
  1. Map Ensembl IDs to gene names via GTF
  2. Remove mitochondrial genes (chrM / MT-)
  3. Remove ribosomal protein genes (RPS*, RPL*)
  4. Remove non-coding / low-information biotypes (rRNA, snoRNA, snRNA, misc_RNA, pseudogenes)
  5. Filter lowly expressed genes (require >= 1 CPM in at least 2 samples)
  6. Compute CPM and log2(CPM + 1)
  7. Write QC summary with filtering stats

Outputs:
  export/
    expression/
      raw_counts_unfiltered.csv     — All 62,754 genes, Ensembl ID + gene name
      raw_counts_filtered.csv       — After QC filtering
      cpm_normalized.csv            — CPM on filtered genes
      log2cpm_normalized.csv        — log2(CPM + 1) on filtered genes
    metadata/
      sample_metadata.csv           — Corrected annotations
    variants/
      notch3_variant_calls.csv
      flt1_isoform_counts.csv
      fn1_isoform_counts.csv
      notch3_isoform_counts.csv
    qc/
      filtering_summary.csv         — Per-step gene removal counts
      library_sizes.csv             — Pre- and post-filter library sizes
      mitochondrial_fraction.csv    — MT read fraction per sample
      sample_correlation.csv        — Pairwise Pearson correlation (log2CPM)
    README.txt
"""

import os
import csv
import math
from collections import defaultdict

WORK = "/sc/arion/work/kellej10/cadasil_public_ipsc"
SCRATCH = "/sc/arion/scratch/kellej10/cadasil_public_ipsc"
EXPORT = os.path.join(WORK, "export")
GTF = os.path.join(SCRATCH, "reference", "gencode.v44.primary_assembly.annotation.gtf")


def build_gene_info(gtf_path):
    """Parse GTF to get gene_id -> (gene_name, gene_type, chromosome)."""
    info = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom = fields[0]
            attrs = fields[8]
            gid = gname = gtype = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    gid = attr.split('"')[1]
                elif attr.startswith('gene_name'):
                    gname = attr.split('"')[1]
                elif attr.startswith('gene_type'):
                    gtype = attr.split('"')[1]
            if gid:
                info[gid] = {"name": gname or "", "type": gtype or "", "chrom": chrom}
    return info


def load_counts(filepath):
    """Load CSV count matrix. Returns header, list of rows."""
    with open(filepath) as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)
    return header, rows


def pearson(x, y):
    """Pearson correlation between two lists."""
    n = len(x)
    mx = sum(x) / n
    my = sum(y) / n
    sx = math.sqrt(sum((xi - mx) ** 2 for xi in x) / n)
    sy = math.sqrt(sum((yi - my) ** 2 for yi in y) / n)
    if sx == 0 or sy == 0:
        return 0.0
    return sum((xi - mx) * (yi - my) for xi, yi in zip(x, y)) / (n * sx * sy)


def main():
    for subdir in ["expression", "metadata", "variants", "qc"]:
        os.makedirs(os.path.join(EXPORT, subdir), exist_ok=True)

    # =========================================================================
    # LOAD DATA
    # =========================================================================
    print("Parsing GTF for gene annotations...")
    gene_info = build_gene_info(GTF)
    print(f"  {len(gene_info)} genes annotated")

    print("Loading raw count matrix...")
    header, rows = load_counts(os.path.join(WORK, "gene_expression_matrix.csv"))
    samples = header[1:]
    n_samples = len(samples)
    print(f"  {len(rows)} genes x {n_samples} samples")

    # =========================================================================
    # STEP 0: Annotate all genes
    # =========================================================================
    annotated = []
    for row in rows:
        eid = row[0]
        gi = gene_info.get(eid, {})
        if not gi:
            base = eid.split(".")[0]
            for k, v in gene_info.items():
                if k.startswith(base):
                    gi = v
                    break
        counts = [int(x) for x in row[1:]]
        annotated.append({
            "ensembl_id": eid,
            "gene_name": gi.get("name", ""),
            "gene_type": gi.get("type", ""),
            "chrom": gi.get("chrom", ""),
            "counts": counts,
        })

    total_start = len(annotated)

    # =========================================================================
    # Write unfiltered raw counts (with gene names)
    # =========================================================================
    unfiltered_header = ["ensembl_id", "gene_name", "gene_type", "chromosome"] + samples
    with open(os.path.join(EXPORT, "expression", "raw_counts_unfiltered.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(unfiltered_header)
        for g in annotated:
            w.writerow([g["ensembl_id"], g["gene_name"], g["gene_type"], g["chrom"]] + g["counts"])
    print(f"\n  Written: expression/raw_counts_unfiltered.csv ({total_start} genes)")

    # =========================================================================
    # STEP 1: Compute mitochondrial fraction BEFORE filtering
    # =========================================================================
    print("\n=== QC: Mitochondrial fraction ===")
    mt_counts = [0] * n_samples
    total_counts = [0] * n_samples
    for g in annotated:
        for j, c in enumerate(g["counts"]):
            total_counts[j] += c
            if g["chrom"] == "chrM":
                mt_counts[j] += c

    mt_frac_header = ["sample", "total_counts", "mt_counts", "mt_fraction_pct"]
    mt_frac_rows = []
    for j, s in enumerate(samples):
        frac = 100 * mt_counts[j] / total_counts[j] if total_counts[j] > 0 else 0
        mt_frac_rows.append([s, total_counts[j], mt_counts[j], f"{frac:.2f}"])
        print(f"  {s}: {frac:.2f}% mitochondrial")

    with open(os.path.join(EXPORT, "qc", "mitochondrial_fraction.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(mt_frac_header)
        w.writerows(mt_frac_rows)

    # =========================================================================
    # STEP 2: Remove mitochondrial genes
    # =========================================================================
    print("\n=== Filtering ===")
    before = len(annotated)
    annotated = [g for g in annotated if g["chrom"] != "chrM"]
    n_mt_removed = before - len(annotated)
    print(f"  Step 1: Removed {n_mt_removed} mitochondrial genes (chrM)")

    # =========================================================================
    # STEP 3: Remove ribosomal protein genes
    # =========================================================================
    before = len(annotated)
    ribo_prefixes = ("RPS", "RPL", "MRPS", "MRPL")
    annotated = [g for g in annotated if not g["gene_name"].startswith(ribo_prefixes)]
    n_ribo_removed = before - len(annotated)
    print(f"  Step 2: Removed {n_ribo_removed} ribosomal protein genes (RPS/RPL/MRPS/MRPL)")

    # =========================================================================
    # STEP 4: Remove non-coding / low-information biotypes
    # =========================================================================
    before = len(annotated)
    remove_biotypes = {
        "rRNA", "rRNA_pseudogene",
        "snoRNA", "snRNA", "scaRNA", "scRNA",
        "misc_RNA", "ribozyme", "vault_RNA", "sRNA",
        "Mt_rRNA", "Mt_tRNA",
        "processed_pseudogene", "unprocessed_pseudogene",
        "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene",
        "translated_processed_pseudogene", "translated_unprocessed_pseudogene",
        "unitary_pseudogene", "polymorphic_pseudogene",
        "pseudogene", "IG_pseudogene", "TR_pseudogene",
        "IG_C_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene",
        "TR_J_pseudogene", "TR_V_pseudogene",
    }
    annotated = [g for g in annotated if g["gene_type"] not in remove_biotypes]
    n_biotype_removed = before - len(annotated)
    print(f"  Step 3: Removed {n_biotype_removed} genes with non-informative biotypes (pseudogenes, rRNA, sn/snoRNA)")

    # =========================================================================
    # STEP 5: Filter lowly expressed genes
    #         Require >= 1 CPM in at least 2 samples
    # =========================================================================
    before = len(annotated)

    # Recompute library sizes after MT removal
    lib_sizes_post = [0] * n_samples
    for g in annotated:
        for j, c in enumerate(g["counts"]):
            lib_sizes_post[j] += c

    filtered = []
    for g in annotated:
        n_above = 0
        for j, c in enumerate(g["counts"]):
            cpm = (c / lib_sizes_post[j]) * 1e6 if lib_sizes_post[j] > 0 else 0
            if cpm >= 1:
                n_above += 1
        if n_above >= 2:
            filtered.append(g)

    n_lowexpr_removed = before - len(filtered)
    annotated = filtered
    print(f"  Step 4: Removed {n_lowexpr_removed} lowly expressed genes (< 1 CPM in >= 2 samples)")
    print(f"\n  Final: {len(annotated)} genes retained (from {total_start} total)")

    # =========================================================================
    # Filtering summary
    # =========================================================================
    filt_summary = [
        ["Starting genes", total_start, ""],
        ["Mitochondrial (chrM)", n_mt_removed, f"{total_start - n_mt_removed} remaining"],
        ["Ribosomal proteins (RPS/RPL/MRPS/MRPL)", n_ribo_removed,
         f"{total_start - n_mt_removed - n_ribo_removed} remaining"],
        ["Non-informative biotypes", n_biotype_removed,
         f"{total_start - n_mt_removed - n_ribo_removed - n_biotype_removed} remaining"],
        ["Lowly expressed (< 1 CPM in 2+ samples)", n_lowexpr_removed,
         f"{len(annotated)} remaining"],
        ["Final filtered gene count", len(annotated), ""],
    ]
    with open(os.path.join(EXPORT, "qc", "filtering_summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["step", "genes_removed", "note"])
        w.writerows(filt_summary)
    print("  Written: qc/filtering_summary.csv")

    # =========================================================================
    # Recompute library sizes on filtered gene set
    # =========================================================================
    lib_sizes = [0] * n_samples
    for g in annotated:
        for j, c in enumerate(g["counts"]):
            lib_sizes[j] += c

    lib_header = ["sample", "raw_total_counts", "filtered_total_counts", "pct_retained"]
    lib_rows = []
    for j, s in enumerate(samples):
        pct = 100 * lib_sizes[j] / total_counts[j] if total_counts[j] > 0 else 0
        lib_rows.append([s, total_counts[j], lib_sizes[j], f"{pct:.1f}"])

    with open(os.path.join(EXPORT, "qc", "library_sizes.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(lib_header)
        w.writerows(lib_rows)
    print("  Written: qc/library_sizes.csv")

    print("\n  Library sizes (filtered):")
    for j, s in enumerate(samples):
        print(f"    {s}: {lib_sizes[j]:,}")

    # =========================================================================
    # WRITE FILTERED RAW COUNTS
    # =========================================================================
    filt_header = ["ensembl_id", "gene_name"] + samples
    with open(os.path.join(EXPORT, "expression", "raw_counts_filtered.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(filt_header)
        for g in annotated:
            w.writerow([g["ensembl_id"], g["gene_name"]] + g["counts"])
    print(f"\n  Written: expression/raw_counts_filtered.csv ({len(annotated)} genes)")

    # =========================================================================
    # CPM NORMALIZATION
    # =========================================================================
    cpm_matrix = []  # for correlation calculation later
    with open(os.path.join(EXPORT, "expression", "cpm_normalized.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(filt_header)
        for g in annotated:
            cpm_vals = [(c / lib_sizes[j]) * 1e6 if lib_sizes[j] > 0 else 0
                        for j, c in enumerate(g["counts"])]
            cpm_matrix.append(cpm_vals)
            w.writerow([g["ensembl_id"], g["gene_name"]] +
                       [f"{v:.4f}" for v in cpm_vals])
    print(f"  Written: expression/cpm_normalized.csv")

    # =========================================================================
    # LOG2(CPM + 1) NORMALIZATION
    # =========================================================================
    log2cpm_matrix = []
    with open(os.path.join(EXPORT, "expression", "log2cpm_normalized.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(filt_header)
        for i, g in enumerate(annotated):
            log_vals = [math.log2(cpm_matrix[i][j] + 1) for j in range(n_samples)]
            log2cpm_matrix.append(log_vals)
            w.writerow([g["ensembl_id"], g["gene_name"]] +
                       [f"{v:.4f}" for v in log_vals])
    print(f"  Written: expression/log2cpm_normalized.csv")

    # =========================================================================
    # SAMPLE CORRELATION MATRIX (Pearson on log2CPM)
    # =========================================================================
    print("\n=== Sample correlation (Pearson, log2CPM) ===")
    # Transpose: sample x gene
    sample_vecs = [[log2cpm_matrix[i][j] for i in range(len(annotated))]
                   for j in range(n_samples)]

    corr_header = ["sample"] + samples
    corr_rows = []
    for j1, s1 in enumerate(samples):
        row = [s1]
        for j2, s2 in enumerate(samples):
            r = pearson(sample_vecs[j1], sample_vecs[j2])
            row.append(f"{r:.4f}")
        corr_rows.append(row)
        print(f"  {s1}: min_r={min(float(x) for x in row[1:]):.3f}")

    with open(os.path.join(EXPORT, "qc", "sample_correlation.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(corr_header)
        w.writerows(corr_rows)
    print("  Written: qc/sample_correlation.csv")

    # =========================================================================
    # BIOTYPE COMPOSITION (top 10)
    # =========================================================================
    biotype_counts = defaultdict(int)
    for g in annotated:
        biotype_counts[g["gene_type"]] += 1
    print("\n  Gene biotype composition (filtered set):")
    for bt, n in sorted(biotype_counts.items(), key=lambda x: -x[1])[:10]:
        print(f"    {bt}: {n}")

    # =========================================================================
    # METADATA
    # =========================================================================
    print("\n=== Metadata ===")
    meta_header = [
        "sample_id", "ena_run", "reported_genotype", "observed_genotype",
        "disease_status", "sample_type", "patient_id",
        "cell_type", "organism", "platform", "notes"
    ]
    meta_rows = [
        ["WZ12_CADASIL_Control", "ERR15965983", "wild_type", "wild_type",
         "normal", "wild_type_control", "Patient_C",
         "iPSC_NC-VSMC", "Homo sapiens", "HiSeq4000_PE76", ""],
        ["WZ4_CADASIL_Control", "ERR15965985", "wild_type", "wild_type",
         "normal", "wild_type_control", "Patient_B",
         "iPSC_NC-VSMC", "Homo sapiens", "HiSeq4000_PE76",
         "Same individual as XZ7 (fingerprint Jaccard=0.46)"],
        ["WZ8_CADASIL_Control", "ERR15965986", "wild_type", "wild_type",
         "normal", "wild_type_control", "Patient_D",
         "iPSC_NC-VSMC", "Homo sapiens", "HiSeq4000_PE76", ""],
        ["WZ16_CADASIL_C224Y", "ERR15965984", "NOTCH3_p.C224Y", "NOTCH3_p.C224Y",
         "CADASIL", "cadasil_mutant", "Patient_A",
         "iPSC_NC-VSMC", "Homo sapiens", "HiSeq4000_PE76",
         "C224Y confirmed (VAF=0.45). Same individual as XZ5 and XZ6."],
        ["XZ5_CADASIL_R153C", "ERR15965987", "NOTCH3_p.R153C", "NOTCH3_p.C224Y",
         "CADASIL", "cadasil_mutant", "Patient_A",
         "iPSC_NC-VSMC", "Homo sapiens", "HiSeq4000_PE76",
         "MISLABELED: reported as R153C but carries C224Y (VAF=0.42). Same individual as WZ16."],
        ["XZ6_CADASIL_R153C_isoCtrl", "ERR15965988", "mutation_corrected_R153C",
         "mutation_corrected_C224Y",
         "normal", "isogenic_control", "Patient_A",
         "iPSC_NC-VSMC", "Homo sapiens", "HiSeq4000_PE76",
         "Isogenic control for Patient_A (C224Y corrected, VAF=0.00)."],
        ["XZ7_CADASIL_C224Y_isoCtrl", "ERR15965989", "mutation_corrected_C224Y",
         "wild_type",
         "normal", "wild_type_control", "Patient_B",
         "iPSC_NC-VSMC", "Homo sapiens", "HiSeq4000_PE76",
         "NOT isogenic for WZ16. Fingerprints match WZ4 (Jaccard=0.46), not WZ16 (0.06)."],
    ]

    with open(os.path.join(EXPORT, "metadata", "sample_metadata.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(meta_header)
        w.writerows(meta_rows)
    print("  Written: metadata/sample_metadata.csv")

    # =========================================================================
    # VARIANTS
    # =========================================================================
    print("\n=== Variant data ===")

    # NOTCH3 variant calls
    with open(os.path.join(WORK, "notch3_variant_summary.csv")) as f:
        reader = csv.reader(f)
        n3h = next(reader)
        n3r = list(reader)
    with open(os.path.join(EXPORT, "variants", "notch3_variant_calls.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(n3h)
        w.writerows(n3r)
    print("  Written: variants/notch3_variant_calls.csv")

    # Split transcript counts by gene
    with open(os.path.join(WORK, "transcript_variant_counts.csv")) as f:
        reader = csv.reader(f)
        tx_header_raw = next(reader)
        tx_rows = list(reader)

    tx_sample_header = ["transcript_id", "annotation"] + tx_header_raw[1:]

    gene_bins = {"FLT1": [], "FN1": [], "NOTCH3": []}
    for row in tx_rows:
        parts = row[0].split("|")
        tx_id, annotation = parts[0], parts[1] if len(parts) > 1 else ""
        clean_row = [tx_id, annotation] + row[1:]
        for gene in gene_bins:
            if gene in annotation:
                gene_bins[gene].append(clean_row)
                break

    descriptions = {
        "FLT1": "FLT1 transcript isoforms: mFlt1 (full-length 30 exons) vs sFlt1 (soluble 13-15 exons)",
        "FN1": "FN1 splice variants: EDA/EDB extra domains (47 exons=both, 44=neither)",
        "NOTCH3": "NOTCH3 transcript isoforms: canonical 33-exon plus truncated/NMD forms",
    }
    filenames = {
        "FLT1": "flt1_isoform_counts.csv",
        "FN1": "fn1_isoform_counts.csv",
        "NOTCH3": "notch3_isoform_counts.csv",
    }

    for gene, gene_rows in gene_bins.items():
        fpath = os.path.join(EXPORT, "variants", filenames[gene])
        with open(fpath, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow([f"# {descriptions[gene]}"])
            w.writerow(tx_sample_header)
            w.writerows(gene_rows)
        print(f"  Written: variants/{filenames[gene]} ({len(gene_rows)} transcripts)")

    # =========================================================================
    # README
    # =========================================================================
    n_final = len(annotated)
    readme = f"""E-MTAB-16303 Export Package — QC'd and Normalized
==================================================
CADASIL iPSC-derived NC-VSMC RNA-seq (Wang, Scholey et al., Science Advances 2025)
DOI: 10.1126/sciadv.aeb1134

Generated: April 14, 2026
Analyst: J. Kelle, Icahn School of Medicine at Mount Sinai

IMPORTANT: Sample provenance discrepancies identified — see metadata notes
and accompanying report (E-MTAB-16303_reanalysis_report.html).


Preprocessing Pipeline
----------------------
  1. Alignment: STAR 2.7.11b to GRCh38 + Gencode v44
  2. Gene counts: featureCounts (Subread 2.0.1), reverse-stranded, paired-end
  3. QC filtering (this script):
     a. Removed {n_mt_removed} mitochondrial genes (chrM)
     b. Removed {n_ribo_removed} ribosomal protein genes (RPS/RPL/MRPS/MRPL)
     c. Removed {n_biotype_removed} non-informative biotypes (pseudogenes, rRNA, sn/snoRNA, etc.)
     d. Removed {n_lowexpr_removed} lowly expressed genes (< 1 CPM in >= 2 samples)
     Final: {n_final} protein-coding and lncRNA genes retained (from {total_start} total)
  4. Normalization: CPM and log2(CPM + 1)


Directory Structure
-------------------

expression/
  raw_counts_unfiltered.csv   All {total_start} genes before filtering
                              Columns: ensembl_id, gene_name, gene_type, chromosome, [samples]

  raw_counts_filtered.csv     {n_final} genes after QC filtering
                              Columns: ensembl_id, gene_name, [samples]
                              Use for: DESeq2, edgeR (tools that expect raw counts)

  cpm_normalized.csv          Counts per million on filtered genes
                              Use for: visualization, cross-sample comparison

  log2cpm_normalized.csv      log2(CPM + 1) on filtered genes
                              Use for: heatmaps, signature scoring, GSEA, correlation

metadata/
  sample_metadata.csv         7 samples with corrected annotations
    Key columns:
      reported_genotype    — what the SDRF says
      observed_genotype    — what we found by variant calling
      disease_status       — normal or CADASIL
      sample_type          — wild_type_control, cadasil_mutant, or isogenic_control
      patient_id           — patient assignment from germline SNP fingerprinting
      notes                — discrepancy details

variants/
  notch3_variant_calls.csv    Ref/alt allele counts + VAF at:
                                p.R153C (chr19:15192182, genomic G>A)
                                p.C224Y (chr19:15191968, genomic C>T)

  flt1_isoform_counts.csv     FLT1 transcript-level counts (featureCounts -g transcript_id)
                                mFlt1 = ENST00000282397 (30 exons, membrane-bound)
                                sFlt1-i13 = ENST00000615840 (13 exons, soluble)
                                sFlt1-i14 = ENST00000639477 (14 exons, soluble)
                                sFlt1-e15a = ENST00000541932 (15 exons, soluble)

  fn1_isoform_counts.csv      FN1 splice variant counts
                                47 exons = EDA+/EDB+ (both extra domains)
                                44 exons = EDA-/EDB- (neither)

  notch3_isoform_counts.csv   NOTCH3 transcript counts
                                NOTCH3-201 = canonical (33 exons, MANE Select)

qc/
  filtering_summary.csv       Genes removed at each step
  library_sizes.csv           Pre- and post-filter library sizes per sample
  mitochondrial_fraction.csv  MT read fraction per sample (before removal)
  sample_correlation.csv      Pairwise Pearson correlation on log2(CPM)


Corrected Sample Groups
------------------------
CADASIL (C224Y):     WZ16_CADASIL_C224Y, XZ5_CADASIL_R153C  (both Patient_A)
Isogenic control:    XZ6_CADASIL_R153C_isoCtrl              (Patient_A, C224Y corrected)
Wild-type controls:  WZ12_CADASIL_Control (Patient_C)
                     WZ4_CADASIL_Control  (Patient_B)
                     WZ8_CADASIL_Control  (Patient_D)
                     XZ7_CADASIL_C224Y_isoCtrl (Patient_B — NOT isogenic for C224Y)
"""

    with open(os.path.join(EXPORT, "README.txt"), "w") as f:
        f.write(readme)
    print("\n  Written: README.txt")
    print(f"\n=== Export complete: {EXPORT} ===")


if __name__ == "__main__":
    main()
