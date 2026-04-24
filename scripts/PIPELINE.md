# E-MTAB-16303 VSMC reanalysis: end-to-end pipeline

Takes the 7 public E-MTAB-16303 paired-end FASTQs (iPSC-derived neural-crest
vascular smooth muscle cells from Zhang et al., *Nature Communications*) and
returns: aligned BAMs, gene + transcript count matrices, per-sample NOTCH3
variant VAFs, NGSCheckMate common-SNP donor clusters, an empirically-corrected
sample metadata table, gene + isoform limma-voom DE results, and the
26-gene plasma-signature interrogation.

This is a validation cohort for the Keller-Lab iEC project; the same methods
live in `cadasil_elab_ipsc_ec/scripts/variants/`. Running on a new batch =
re-point the sample manifest and the SNP-panel symlink.

```
         ┌──▶ 00a salmon index (one-time, shared) ─┐
FASTQ ──┤                                          ├──▶ 00b salmon quant ──▶ 16 isoform DE
         └──▶ 01 align ─┬──▶ 02 count + pileup ──▶ 03 build matrices
                        └──▶ 07 genotype_snps ──▶ 08 cluster donors ──▶ 09 reassign
                                                                          │
             ┌────────────────────────────────────────────────────────────┘
             ▼
         16 isoform DE · 26 gene DE + signature interrogation · 24 annotate isoforms
```

| Sample | ERR | Submitted notch3_mutation | Inferred donor | Inferred genotype |
|---|---|---|---|---|
| WZ12_CADASIL_Control | ERR15965983 | none | (set after 09) | (set after 09) |
| WZ16_CADASIL_C224Y | ERR15965984 | p.C224Y | C224Y | CAD |
| WZ4_CADASIL_Control | ERR15965985 | none | (set after 09) | WT |
| WZ8_CADASIL_Control | ERR15965986 | none | (set after 09) | WT |
| XZ5_CADASIL_R153C | ERR15965987 | p.R153C | **C224Y** | **CAD** (errata: submitted R153C is wrong) |
| XZ6_CADASIL_R153C_isoCtrl | ERR15965988 | corrected_R153C | R153C | ISO |
| XZ7_CADASIL_C224Y_isoCtrl | ERR15965989 | corrected_C224Y | C224Y | ISO |

---

## 0. Prerequisites (one-time)

| Item | Location |
|---|---|
| Reference FASTA | `$SCRATCH/reference/GRCh38.primary_assembly.genome.fa` |
| GTF | `gencode.v44.primary_assembly.annotation.gtf` |
| STAR index | `$SCRATCH/reference/star_index` (sjdbOverhang=75) |
| Salmon quants | `data/salmon_quant/<sample>/quant.sf` (7 dirs) |
| NGSCheckMate SNP panel | `$SCRATCH/snp_panel/snp_targets.tsv` (~21k biallelic common SNPs, symlinked from iEC) |
| tx2gene map | `data/tx2gene_gencode_v44.csv` (symlinked from iEC) |
| Plasma signature | `data/signature_unified_betas.csv` (symlinked from iEC) |
| Submitted metadata | `sample_metadata.csv` (sample_name, err_id, notch3_mutation_reported, ...) |

`$SCRATCH = /sc/arion/scratch/kellej10/cadasil_public_ipsc`. Scratch is
purged after 14 days — final outputs land in `results/`, `data/`, `export/`.

---

## 1. Alignment and counting

```bash
bash   scripts/00_download_references.sh    # reference FASTA + GTF
bsub < scripts/01_build_star_index.lsf      # STAR index (sjdbOverhang=75)
bsub < scripts/02_align_and_count.lsf       # STAR align + featureCounts + pileup
bsub < scripts/02b_count_and_pileup.lsf     # transcript-level counts + NOTCH3 pileup
```

Outputs land in `$SCRATCH/aligned/` (BAMs) and `$SCRATCH/counts/`.

---

## 2. Build matrices

```bash
module purge && module load anaconda3/2024.06 && conda activate orion-minimal
python scripts/03_build_matrices.py
```

Writes `gene_expression_matrix.csv`, `transcript_variant_counts.csv`,
`notch3_variant_summary.csv`, `genes_of_interest_counts.csv` at the repo
root, then `04_build_export_package.py` bundles them into
`export/` for the manuscript supplement.

---

## 3. NGSCheckMate donor clustering (authoritative donor identity)

```bash
bsub < scripts/07_genotype_snps.lsf                 # bcftools mpileup + call across all 7 BAMs
python scripts/convert_notch3_long.py               # wide notch3_variant_summary.csv -> long
bsub < scripts/08_09_cluster_and_assign.lsf         # cluster + reassign (wraps 08 + 09)
```

### 07 – `07_genotype_snps.lsf`
`bcftools mpileup -T snp_panel/snp_targets.tsv.gz -a FORMAT/AD,FORMAT/DP
-q 20 -Q 20` → `bcftools call -m` at ~21k biallelic common SNPs. Writes
`$SCRATCH/genotype/vsmc_panel_genotypes.{vcf.gz,tsv}`.

### 08 – `08_cluster_donors.py`
For each pair of samples, computes IBS and Jaccard_het at panel sites
(DP≥10). Complete-linkage hierarchical cluster on `1 − Jaccard_het`, cut at
distance 0.15 (within-donor J ≈ 0.95; across-donor J ≈ 0.30). Writes
`results/variants/panel_{ibs,jaccard}_matrix.csv`,
`panel_donor_clusters.csv`, `panel_dendrogram_merges.csv`.

### 09 – `09_assign_genotype.py`
Each cluster is a donor. Picks the CAD driver variant (R153C or C224Y)
per cluster by requiring ≥1 sample with VAF ≥ 0.15 & DP ≥ 5; tie-break by
highest VAF. Within a CAD cluster, classifies each sample:
- VAF ≥ 0.15 & DP ≥ 5 → **CAD**
- VAF ≤ 0.05 & DP ≥ 5 → **ISO**
- 0.05 < VAF < 0.15 → **AMBIGUOUS**
- DP < 5 → **LOW_COV**

Pure-WT clusters are labelled `<donor>_noVariant` (multi-member) or
`SINGLETON_<donor>` (single-member). There is no KOLF2.1J reference set
in this cohort — WZ12/WZ4/WZ8 are independent WTs.

Writes `data/sample_metadata_reassigned.csv` — all downstream steps read
this, not `sample_metadata.csv`.

**Errata**: XZ5 clusters with the C224Y donor (not R153C), and its
NOTCH3 pileup shows VAF ≈ 0.42 for C224Y, near-zero for R153C. The
submitted `p.R153C` annotation is a sample-swap in the public metadata.

---

## 4. Differential expression

```bash
bsub < scripts/20_run_de.lsf                         # wraps 16 + 26 + 24
```

### 16 – `16_isoform_de.R`
tximport (txOut=TRUE, countsFromAbundance=scaledTPM) → edgeR `filterByExpr`
(min.count=5, min.total.count=10) → TMM → `voom` → `lmFit` → `eBayes`.
Two contrasts:
- `vsmc_pooled_CAD_vs_nonCAD` — all 7 samples, CAD vs everything else
- `vsmc_C224Y_cluster_CAD_vs_ISO` — 3 samples, tests within-donor isogenic pair

No donor blocking (donor-condition confounding is absolute in this cohort).
Writes per-contrast transcript DE + per-sample log2CPM to
`results/de_reassigned/isoform/`.

### 24 – `24_annotate_isoforms.py`
Joins the 16 panel output to `_isoform_annotations.ANNOT` (shared with
iEC repo: FLT1 mFlt1/sFlt1, FLT4 m/sFlt4, FN1 EIIIA/EIIIB/V alternative
exons) and writes `{FLT1,FLT4,FN1}_isoform_stats_annotated.csv`.

### 26 – `26_vsmc_signature_interrogation.R`
Gene-level limma-voom DE on the same two contrasts, then tests the 26-gene
plasma signature (`data/signature_unified_betas.csv`) for direction
concordance. Outputs `signature_interrogation_vsmc*.csv`.

---

## Outputs

```
data/
  sample_metadata_reassigned.csv   # authoritative metadata (step 09)
  salmon_quant/<sample>/quant.sf   # per-sample isoform quants
results/
  variants/
    panel_ibs_matrix.csv           # fingerprint: pairwise IBS
    panel_jaccard_matrix.csv       # fingerprint: pairwise Jaccard_het (primary)
    panel_donor_clusters.csv       # sample -> donor cluster
    panel_dendrogram_merges.csv    # hclust merge trace
    notch3_variant_summary.csv     # long-format VAF per sample x variant
    cohort_reassigned_summary.md   # human-readable cluster / genotype report
  de_reassigned/
    isoform/tx_de_<contrast>_{all,panel}.csv
    isoform/tx_logcpm_<contrast>.csv
    isoform/{FLT1,FLT4,FN1}_isoform_stats_annotated.csv
    de_<contrast>.csv              # gene-level DE (step 26)
    signature_interrogation_vsmc*.csv
```
