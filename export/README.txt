E-MTAB-16303 Export Package — QC'd and Normalized
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
     a. Removed 37 mitochondrial genes (chrM)
     b. Removed 1881 ribosomal protein genes (RPS/RPL/MRPS/MRPL)
     c. Removed 18590 non-informative biotypes (pseudogenes, rRNA, sn/snoRNA, etc.)
     d. Removed 27639 lowly expressed genes (< 1 CPM in >= 2 samples)
     Final: 14607 protein-coding and lncRNA genes retained (from 62754 total)
  4. Normalization: CPM and log2(CPM + 1)


Directory Structure
-------------------

expression/
  raw_counts_unfiltered.csv   All 62754 genes before filtering
                              Columns: ensembl_id, gene_name, gene_type, chromosome, [samples]

  raw_counts_filtered.csv     14607 genes after QC filtering
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
