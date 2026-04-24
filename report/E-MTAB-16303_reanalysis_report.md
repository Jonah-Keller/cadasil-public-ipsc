# Independent Reanalysis of E-MTAB-16303: Sample Provenance Discrepancies in CADASIL iPSC RNA-seq Data

**Date:** April 14, 2026  
**Analyst:** J. Kelle, Icahn School of Medicine at Mount Sinai  
**Dataset:** E-MTAB-16303 (ArrayExpress/BioStudies)  
**Associated publication:** Wang, Scholey et al., "Selective vulnerability of cerebral vasculature to NOTCH3 variants in small vessel disease and rescue by phosphodiesterase-5 inhibitor," *Science Advances* (2025). DOI: 10.1126/sciadv.aeb1134  

---

## 1. Summary

During independent reanalysis of the RNA-seq dataset E-MTAB-16303, we identified discrepancies between the reported sample genotypes and the genotypes observed in the sequencing data. Specifically:

1. Sample XZ5, annotated as carrying NOTCH3 p.R153C (c.457C>T), instead carries the p.C224Y (c.671G>A) mutation — the same mutation attributed to sample WZ16.
2. Germline SNP fingerprinting reveals that XZ5 and WZ16 originate from the same patient, and that the reported isogenic control pairing (WZ16 ↔ XZ7) is not supported by the data.

These findings suggest sample mislabeling in the deposited metadata and may affect the interpretation of mutation-specific analyses in the associated publication.

---

## 2. Data Acquisition

### 2.1 Study overview

E-MTAB-16303 describes bulk RNA-seq of iPSC-derived neural crest vascular smooth muscle cells (NC-VSMCs) from CADASIL patients and controls. The study deposited 7 samples with paired-end 76 bp reads generated on an Illumina HiSeq 4000 at the Genomic Technologies Core Facility, University of Manchester.

### 2.2 Sample metadata as reported in SDRF

The following sample annotations were obtained directly from the study's SDRF file (E-MTAB-16303.sdrf.txt):

| Sample Name | ENA Run | Reported Genotype | Reported Disease | Factor Value |
|---|---|---|---|---|
| WZ12_CADASIL_Control | ERR15965983 | wild type genotype | normal | normal |
| WZ4_CADASIL_Control | ERR15965985 | wild type genotype | normal | normal |
| WZ8_CADASIL_Control | ERR15965986 | wild type genotype | normal | normal |
| WZ16_CADASIL_C224Y | ERR15965984 | NOTCH3 p.C224Y mutation | CADASIL | CADASIL |
| XZ5_CADASIL_R153C | ERR15965987 | NOTCH3 p.R153C mutaton [sic] | CADASIL | CADASIL |
| XZ6_CADASIL_R153C_isoCtrl | ERR15965988 | mutation corrected | normal | normal |
| XZ7_CADASIL_C224Y_isoCtrl | ERR15965989 | mutation corrected | normal | normal |

Note: The SDRF contains a typographical error in the XZ5 genotype field ("mutaton" instead of "mutation").

### 2.3 Raw data retrieval

Unaligned BAM files were retrieved from the European Nucleotide Archive (ENA) under project accession ERP185733:

```
ftp://ftp.sra.ebi.ac.uk/vol1/err/ERR159/083/ERR15965983/ERR15965983.bam
ftp://ftp.sra.ebi.ac.uk/vol1/err/ERR159/084/ERR15965984/ERR15965984.bam
ftp://ftp.sra.ebi.ac.uk/vol1/err/ERR159/085/ERR15965985/ERR15965985.bam
ftp://ftp.sra.ebi.ac.uk/vol1/err/ERR159/086/ERR15965986/ERR15965986.bam
ftp://ftp.sra.ebi.ac.uk/vol1/err/ERR159/087/ERR15965987/ERR15965987.bam
ftp://ftp.sra.ebi.ac.uk/vol1/err/ERR159/088/ERR15965988/ERR15965988.bam
ftp://ftp.sra.ebi.ac.uk/vol1/err/ERR159/089/ERR15965989/ERR15965989.bam
```

BAM header inspection confirmed these were unaligned BAMs (queryname-sorted, no @SQ reference sequences), requiring full alignment processing.

---

## 3. Alignment and Quantification Pipeline

### 3.1 Reference genome

- **Genome assembly:** GRCh38 primary assembly
- **Gene annotation:** Gencode v44 (gencode.v44.primary_assembly.annotation.gtf)
- **Source:** https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/

### 3.2 FASTQ extraction

Paired-end FASTQ files were extracted from unaligned BAMs using samtools v1.21:

```bash
samtools fastq -@ 4 -1 ${ERR}_1.fastq.gz -2 ${ERR}_2.fastq.gz \
    -0 /dev/null -s /dev/null ${ERR}.bam
```

Read counts per sample ranged from 59.7M to 77.4M total reads (29.8M–38.7M read pairs).

### 3.3 Alignment

Reads were aligned to the GRCh38 primary assembly using STAR v2.7.11b:

```bash
STAR --runMode alignReads \
     --runThreadN 8 \
     --genomeDir ${STAR_IDX} \
     --readFilesIn ${FQ1} ${FQ2} \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes NH HI AS NM MD \
     --quantMode GeneCounts TranscriptomeSAM \
     --sjdbOverhang 75
```

The STAR genome index was built with `--sjdbOverhang 75` (read length 76 - 1).

### 3.4 Gene-level quantification

Read counting was performed with featureCounts (Subread v2.0.1):

```bash
featureCounts -a ${GTF} -o gene_counts.txt \
    -g gene_id -t exon -p -B -s 2 -T 8 --primary ${BAMS}
```

Parameters: paired-end fragment counting (`-p -B`), reverse-stranded library (`-s 2`, consistent with TruSeq Stranded mRNA protocol), primary alignments only.

### 3.5 Alignment quality

| Sample | Total Alignments | Assigned Reads | Assignment Rate |
|---|---|---|---|
| WZ12_CADASIL_Control | 45,507,221 | 33,530,419 | 73.7% |
| WZ16_CADASIL_C224Y | 33,930,321 | 25,956,674 | 76.5% |
| WZ4_CADASIL_Control | 37,608,869 | 27,844,794 | 74.0% |
| WZ8_CADASIL_Control | 36,652,400 | 28,314,187 | 77.3% |
| XZ5_CADASIL_R153C | 37,790,891 | 31,702,955 | 83.9% |
| XZ6_CADASIL_R153C_isoCtrl | 35,497,383 | 30,039,444 | 84.6% |
| XZ7_CADASIL_C224Y_isoCtrl | 31,160,908 | 26,240,819 | 84.2% |

All samples showed high alignment rates and gene assignment rates (73.7–84.6%), indicating good library quality.

---

## 4. NOTCH3 Variant Detection

### 4.1 Approach

We performed targeted variant detection at known NOTCH3 mutation sites using samtools mpileup, following a targeted pileup strategy analogous to approaches used for somatic variant detection in RNA-seq data.

NOTCH3 (ENSG00000074181) is located on the minus strand of chromosome 19 (GRCh38: chr19:15,159,038–15,200,995). The canonical transcript is ENST00000263388.7 (NOTCH3-201, 33 exons, 8,680 bp).

### 4.2 Variant coordinates

Genomic coordinates for the two reported CADASIL mutations were obtained from the Ensembl Variant Effect Predictor (VEP) using transcript NM_000435.3:

| Mutation | cDNA Change | Genomic Position (GRCh38) | Exon | Genomic Ref/Alt (minus strand) |
|---|---|---|---|---|
| p.R153C | c.457C>T | chr19:15,192,182 | Exon 4 | G>A |
| p.C224Y | c.671G>A | chr19:15,191,968 | Exon 4 | C>T |

Both mutations fall within exon 4 of NOTCH3 (chr19:15,191,960–15,192,298, 339 bp). Because NOTCH3 is transcribed from the minus strand, cDNA nucleotide changes are reverse-complemented relative to the genomic reference.

### 4.3 Pileup methodology

Pileup was performed across the full NOTCH3 gene region for all 7 samples:

```bash
samtools mpileup -r chr19:15159038-15200995 \
    -f GRCh38.primary_assembly.genome.fa \
    -q 20 -Q 20 --no-BAQ ${BAM}
```

Variant allele frequency (VAF) was calculated as: alt_reads / total_depth.

### 4.4 Results at reported mutation sites

**Position chr19:15,192,182 (reported p.R153C site, c.457C>T, genomic G>A):**

| Sample | Reported Genotype | Depth | Ref (G) | Alt (A) | VAF |
|---|---|---|---|---|---|
| WZ12_CADASIL_Control | wild type | 174 | 174 | 0 | 0.000 |
| WZ4_CADASIL_Control | wild type | 27 | 27 | 0 | 0.000 |
| WZ8_CADASIL_Control | wild type | 162 | 162 | 0 | 0.000 |
| WZ16_CADASIL_C224Y | NOTCH3 p.C224Y | 162 | 162 | 0 | 0.000 |
| **XZ5_CADASIL_R153C** | **NOTCH3 p.R153C** | **485** | **485** | **0** | **0.000** |
| XZ6_CADASIL_R153C_isoCtrl | mutation corrected | 90 | 90 | 0 | 0.000 |
| XZ7_CADASIL_C224Y_isoCtrl | mutation corrected | 320 | 320 | 0 | 0.000 |

**Finding:** Despite excellent coverage (485x), sample XZ5 shows **zero variant alleles** at the reported R153C position. All reads match the reference (G), with no evidence of the expected A allele.

**Position chr19:15,191,968 (reported p.C224Y site, c.671G>A, genomic C>T):**

| Sample | Reported Genotype | Depth | Ref (C) | Alt (T) | VAF |
|---|---|---|---|---|---|
| WZ12_CADASIL_Control | wild type | 283 | 283 | 0 | 0.000 |
| WZ4_CADASIL_Control | wild type | 46 | 46 | 0 | 0.000 |
| WZ8_CADASIL_Control | wild type | 291 | 291 | 0 | 0.000 |
| **WZ16_CADASIL_C224Y** | **NOTCH3 p.C224Y** | **251** | **139** | **112** | **0.454** |
| **XZ5_CADASIL_R153C** | **NOTCH3 p.R153C** | **464** | **271** | **193** | **0.420** |
| XZ6_CADASIL_R153C_isoCtrl | mutation corrected | 72 | 72 | 0 | 0.000 |
| XZ7_CADASIL_C224Y_isoCtrl | mutation corrected | 270 | 270 | 0 | 0.000 |

**Finding:** Both CADASIL samples (WZ16 and XZ5) carry the **C224Y mutation** (VAF ~0.45, consistent with heterozygosity). XZ5, annotated as carrying p.R153C, instead carries p.C224Y. The isogenic controls (XZ6 and XZ7) correctly show no variant alleles at either position, confirming successful CRISPR correction.

### 4.5 Genome-wide scan for R153C

To rule out coordinate errors, we scanned the entire NOTCH3 gene region in XZ5 for any position with a heterozygous variant allele fraction (VAF 0.3–0.7, depth ≥ 20). No position within the NOTCH3 coding sequence showed a heterozygous variant unique to XZ5 that was absent from WZ16 and controls, confirming that XZ5 does not carry a private NOTCH3 coding mutation distinct from C224Y.

---

## 5. Sample Provenance Verification by Germline SNP Fingerprinting

### 5.1 Rationale

If XZ5 (labeled R153C) and WZ16 (labeled C224Y) carry the same NOTCH3 mutation, they may originate from the same patient. Conversely, the reported isogenic control pairings (XZ5 ↔ XZ6 and WZ16 ↔ XZ7) should be verifiable by shared germline heterozygous SNPs, since CRISPR correction of a single locus preserves the genome-wide SNP background.

### 5.2 Method

We performed samtools mpileup across 12 gene regions on different chromosomes (chr1, chr2, chr3, chr5, chr7, chr9, chr11, chr12, chr15, chr17, chr19, chr22) for all 7 samples. Heterozygous germline SNP positions were identified as sites with depth ≥ 15 and VAF between 0.2–0.8 with at least 5 non-reference reads.

Pairwise sample similarity was quantified using the Jaccard index:

```
J(A, B) = |shared het positions| / |union of het positions|
```

### 5.3 Results

| Sample Pair | Shared | Union | Jaccard | Relationship |
|---|---|---|---|---|
| **WZ16_C224Y ↔ XZ5_R153C** | **13** | **26** | **0.500** | **Same individual** |
| **XZ5_R153C ↔ XZ6_R153C_isoCtrl** | **12** | **25** | **0.480** | **Same cell line (isogenic pair)** |
| **WZ4_Control ↔ XZ7_C224Y_isoCtrl** | **18** | **39** | **0.462** | **Same individual** (unexpected) |
| WZ16_C224Y ↔ XZ6_R153C_isoCtrl | 7 | 21 | 0.333 | Moderate (consistent with WZ16 ≈ XZ5 ≈ XZ6) |
| WZ16_C224Y ↔ XZ7_C224Y_isoCtrl | 3 | 47 | **0.064** | **Different individuals** |
| WZ12 ↔ WZ4 | 8 | 44 | 0.182 | Different individuals |
| WZ12 ↔ WZ8 | 1 | 41 | 0.024 | Different individuals |
| All other control pairs | — | — | <0.19 | Different individuals |

### 5.4 Interpretation

The fingerprinting data reveals the following sample relationships:

**Cluster 1 (Patient A):** WZ16_C224Y, XZ5_R153C, and XZ6_R153C_isoCtrl are derived from the **same individual**. WZ16 and XZ5 share the highest Jaccard similarity (0.500), and XZ6 clusters with both. All three carry or carried the C224Y mutation. The "R153C" label on XZ5 and XZ6 is not supported by the data.

**Cluster 2 (Patient B):** WZ4_Control and XZ7_C224Y_isoCtrl are derived from the **same individual** (Jaccard 0.462). This is unexpected because XZ7 is annotated as an isogenic control for C224Y (implying derivation from WZ16), but the SNP fingerprint matches WZ4 instead.

**Independent samples:** WZ12_Control and WZ8_Control are each from different individuals, with low Jaccard scores against all other samples.

**The reported pairing WZ16 ↔ XZ7 is NOT supported** (Jaccard 0.064), indicating these samples are from different individuals.

---

## 6. Summary of Discrepancies

| Issue | Reported | Observed | Evidence |
|---|---|---|---|
| XZ5 genotype | NOTCH3 p.R153C | **NOTCH3 p.C224Y** | VAF = 0.42 at C224Y site; VAF = 0.00 at R153C site |
| XZ5 provenance | Independent patient from WZ16 | **Same patient as WZ16** | Germline SNP Jaccard = 0.500 |
| XZ6 label | Isogenic control for R153C | Isogenic control for **C224Y** (same patient as WZ16/XZ5) | Germline SNP Jaccard with XZ5 = 0.480 |
| XZ7 provenance | Isogenic control from WZ16 (C224Y) | Derived from **WZ4** (wild type) | Jaccard with WZ4 = 0.462; Jaccard with WZ16 = 0.064 |

### 6.1 Revised sample relationships

```
Patient A:  WZ16_C224Y ──── XZ5 (mislabeled "R153C", actually C224Y) ──── XZ6 (isoCtrl)
            All three from the same individual. Both WZ16 and XZ5 carry p.C224Y.
            XZ6 is the CRISPR-corrected isogenic control (confirmed: no variant).

Patient B:  WZ4_Control ──── XZ7 (labeled "C224Y_isoCtrl", actually WT from WZ4 background)
            Both from the same individual. Neither carries a NOTCH3 mutation.

Patient C:  WZ12_Control (independent wild-type)

Patient D:  WZ8_Control (independent wild-type)
```

### 6.2 Impact on study conclusions

The dataset provides:
- **2 biological replicates** of C224Y mutant NC-VSMCs (WZ16 and XZ5), not one C224Y + one R153C as reported
- **1 true isogenic control** for C224Y (XZ6), not two as reported
- **4 wild-type controls** (WZ12, WZ4, WZ8, and XZ7), rather than 3 WT + 2 isogenic
- **No R153C data** in this dataset

Differential expression analyses comparing "CADASIL vs. control" remain valid (both mutation samples are true C224Y heterozygotes vs. true wild-type/corrected controls). However, any mutation-specific comparisons between R153C and C224Y, or analyses leveraging two independent NOTCH3 mutations as biological replicates of the CADASIL phenotype, would need to be revisited.

---

## 7. Methods Reproducibility

All analyses were performed on the Minerva HPC cluster (Icahn School of Medicine at Mount Sinai).

### Software versions
- samtools 1.21
- STAR 2.7.11b
- Subread (featureCounts) 2.0.1
- Python 3.x (standard library only)

### Reference data
- GRCh38 primary assembly genome (Gencode)
- Gencode v44 primary assembly annotation GTF

### Code availability
All analysis scripts are available at:
`/sc/arion/work/kellej10/cadasil_public_ipsc/scripts/`

Key scripts:
- `00_download_references.sh` — Reference genome and annotation download
- `01_build_star_index.lsf` — STAR genome index generation
- `02_align_and_count.lsf` — FASTQ extraction, alignment, and quantification
- `02b_count_and_pileup.lsf` — featureCounts and NOTCH3 variant pileup
- `03_build_matrices.py` — Expression matrix and metadata assembly
- `check_notch3_variants.sh` — NOTCH3 variant position scanning
- `compare_het_positions.sh` — Mutation-specific heterozygous position comparison
- `sample_fingerprint.lsf` — Germline SNP fingerprinting and pairwise concordance

### Data outputs
- `gene_expression_matrix.csv` — 62,754 genes x 7 samples (raw counts)
- `sample_metadata.csv` — Sample annotations with reported and observed genotypes
- `transcript_variant_counts.csv` — Transcript-level counts for FLT1, FN1, NOTCH3 isoforms
- `notch3_variant_summary.csv` — Per-sample VAF at R153C and C224Y positions
- `genes_of_interest_counts.csv` — Counts for key VSMC and CADASIL-related genes

---

## 8. Contact

For questions regarding this reanalysis, please contact:

J. Kelle  
Icahn School of Medicine at Mount Sinai  
Vascular Brain Health Laboratory  

---

*This report was generated on April 14, 2026.*
