# Task document — port `cadasil_elab_ipsc_ec` strategies into this repo

**Status:** drafted 2026-04-24, for a follow-on Claude agent to execute.
**Driver:** Jonah will run a fresh Claude against this file; the agent should
read every section, then work top-down and open a PR at the end.

---

## 0. Background: the two repos

| Repo | Role | Cohort |
|---|---|---|
| `cadasil_public_ipsc/` (this repo) | **Progenitor** — first pass on the Zhang et al. public VSMC dataset. Established the STAR + featureCounts + mpileup variant approach. | 7 public samples, E-MTAB-16303 |
| `cadasil_elab_ipsc_ec/` | **Reference implementation** — built by porting this repo forward, then extended with SNP-panel donor inference, Salmon closed-loop, limma-voom DE, cross-cohort validation. Has since re-run **this** 7-sample cohort under `vsmc_validation/` as a validation branch. | 33 iEC samples (main) + the same 7 VSMC samples (validation) |

**Goal of this port:** bring this repo up to parity with the reference
implementation on three fronts — (a) **donor identity via NGSCheckMate SNP
panel**, (b) **Salmon-quant-based isoform DE**, (c) **limma-voom DE with a
shared annotation module** — so findings from the public data (the XZ5 R153C
→ C224Y errata call; sFLT1-i13 behavior; FN1 EDA splice) can be reproduced
from FASTQ with the same rigor as the lab cohort.

**Scope:** this repo stays a single-cohort repo. Do NOT import the `cohorts/`
split — overkill for 7 samples. But DO adopt the same **scripts/variants/**
numbering so both repos look alike.

---

## 1. Current state of this repo (read before editing)

Scripts present in `scripts/`:

```
00_download_references.sh           ok — keep
01_build_star_index.lsf             ok — keep (but see §3 on Salmon index)
02_align_and_count.lsf              keep; will depend on a new Salmon quant step
02b_count_and_pileup.lsf            keep (feeds NOTCH3 focal-site reports)
03_build_matrices.py                keep structure; see §4 on annotation import
04_build_export_package.py          keep — this is the publish-ready bundle
05_isoform_comparison.py            REPLACE (ad-hoc CPM ratios; switch to
                                    proper limma-voom on Salmon quants)
sample_fingerprint.lsf              REPLACE (housekeeping-het Jaccard;
                                    switch to NGSCheckMate SNP panel flow)
build_fingerprint_csv.sh            DELETE after replacement lands
check_notch3_variants.sh            keep (useful ad-hoc sanity check)
check_r153c_position.sh             keep
compare_het_positions.sh            DELETE (obsolete after SNP-panel port)
extract_transcripts.sh              keep
RUN_PIPELINE.txt                    empty → write a real pipeline.md (§7)
```

Reference implementation to port FROM:
`/sc/arion/work/kellej10/cadasil_elab_ipsc_ec/scripts/variants/`
**Read `PIPELINE.md` in that directory first. It is the source of truth for
what the target architecture looks like.**

---

## 2. Task list (do in order; each task opens a self-contained PR)

### Task 1 — Initialize git + bring settings forward

1. `git init` at the repo root; `.gitignore` should exclude: `*.bam*`, `*.bai`, `*.fastq*`, `aligned/`, `counts/`, `logs/`, `snp_panel/raw/`, `salmon_quant/`, `*.zip` (the export bundles are large binaries and should not go in git; keep the regeneration scripts instead).
2. Copy `CLAUDE.md` conventions from the iEC repo if they're helpful.
3. First commit: `Initial: port planning; prior loose files tracked as baseline`.

### Task 2 — Port NGSCheckMate SNP-panel genotyping

**Replaces:** `sample_fingerprint.lsf` + `build_fingerprint_csv.sh` + `compare_het_positions.sh`.

**Why:** the current housekeeping-het Jaccard conflates expression noise with genotype. NGSCheckMate uses a curated panel of ~21k biallelic common SNPs, giving a clean IBS / Jaccard_het matrix that resolves same-donor pairs at J > 0.9 and unrelated pairs at J < 0.4. The iEC repo uses this to independently reproduce the Zhang errata call (XZ5 clusters with WZ16 on the C224Y allele, not with XZ6 as R153C would require).

**Source files to port (from iEC repo):**
- `scripts/variants/07_genotype_snps.lsf` — `bcftools mpileup -T snp_targets.tsv.gz -b bam_list.txt -a FORMAT/AD,FORMAT/DP` then `bcftools call -m -T snp_targets.tsv.gz`. Output: multi-sample VCF + wide TSV `panel_genotypes.tsv`.
- `scripts/variants/08_cluster_donors.py` — parses the wide TSV; keeps sites with DP ≥ 10 per sample AND called in ≥ 3 samples; encodes `0 / 1 / 2 / -1`; computes pairwise IBS and Jaccard_het; complete-linkage hierarchical cluster on `1 - Jaccard_het` cut at `distance = 0.15`; writes `panel_ibs_matrix.csv`, `panel_jaccard_matrix.csv`, `panel_call_counts.csv`, `panel_dendrogram_merges.csv`, `panel_donor_clusters.csv`.
- `scripts/variants/09_assign_genotype.py` — joins cluster assignments with the focal-site VAF table (NOTCH3 R153C + C224Y here); writes `sample_metadata_reassigned.csv` with columns `donor_inferred`, `genotype_inferred`, `driver_variant`, `driver_vaf`, `driver_depth`, `donor_mismatch`, `genotype_mismatch`. Rules: VAF ≥ 0.15 at DP ≥ 5 → CAD; VAF ≤ 0.05 at DP ≥ 5 → ISO; in-between → AMBIGUOUS; DP < 5 → LOW_COV; no driver in cluster → WT.
- `scripts/variants/08_09_cluster_and_assign.lsf` — the LSF wrapper that chains 08 + 09.

**NGSCheckMate SNP panel file:** the iEC repo pulls it from `$SCRATCH/snp_panel/snp_targets.tsv.gz` (~21k biallelic common SNPs). It is NOT in git. Symlink from the iEC scratch tree:
```
ln -s /sc/arion/scratch/kellej10/cadasil_elab_ipsc_ec/snp_panel \
      /sc/arion/scratch/kellej10/cadasil_public_ipsc/snp_panel
```
(both cohorts use the exact same panel on GRCh38; no regeneration needed.)

**Cohort-specific tweaks when porting:**
- `COHORT_DONOR_VARIANTS="R153C C224Y"` (this cohort has no R544C/C87R/G420C donors).
- 7 samples, not 33 — the cluster-cut threshold 0.15 is still correct, but verify the dendrogram visually by running `08` on the first pass.
- `09_assign_genotype.py` has a KOLF2.1J reference hard-coded for the iEC cohort (`A17 A18 A26 A31`). Strip that path for this cohort, or guard it with a list that defaults empty. WZ12 / WZ4 / WZ8 are independent WTs, not isogenic controls — they should cluster each in their own singleton and get `WT` label with no `_noVariant` suffix.

**Acceptance test:** after running 07 → 08_09, `sample_metadata_reassigned.csv` should show:
- XZ5 → `donor_inferred=C224Y`, `genotype_inferred=CAD`, `driver_variant=p.C224Y`, `driver_vaf≈0.42`.
- WZ16 → same cluster as XZ5 and XZ6 at Jaccard_het > 0.9 (all three are the same donor).
- XZ7 + any associated cluster members → `donor_inferred=C224Y`, `genotype_inferred=ISO`.
- WZ12, WZ4, WZ8 → three independent clusters, each `WT`.

This reproduces the Zhang errata **independently from the paper text**.

### Task 3 — Port Salmon closed-loop (index + per-sample quant)

**Why:** featureCounts at the transcript level distributes ambiguous reads equally, so any pair of transcripts that share all their exons gets near-identical counts. `05_isoform_comparison.py` currently consumes these counts, which means the sFLT1/mFLT1 and EDA+/EDA- ratios it reports are **partially an artifact of count distribution, not true isoform-resolved abundance**. Salmon's EM algorithm fixes this.

**Source files to port:**
- `scripts/variants/00a_build_salmon_index.lsf` — decoy-aware gencode v44, k=31, `module load salmon/1.4.0`. **This is a one-time build; the index built for the iEC repo is at `$PUB_REF_DIR/salmon_index_v44` and is fine to reuse.** If the index directory already exists at `/sc/arion/scratch/kellej10/cadasil_public_ipsc/reference/salmon_index_v44`, you can skip the rebuild and just document the path.
- `scripts/variants/00b_salmon_quant.lsf` — array job, one quant per sample, `--libType A`, outputs `data/salmon_quant/<sample>/quant.sf`. Drop the cohort-config layer — just hard-code the 7 samples.

**Library type caveat:** this cohort is stranded reverse (TruSeq stranded mRNA). Salmon `--libType A` auto-detects, but verify `lib_format_counts.json` reports `ISR` for the first sample before running the whole array.

**Acceptance test:** `data/salmon_quant/<sample>/quant.sf` exists for all 7 samples and has ~250k transcripts per file.

### Task 4 — Port limma-voom DE on Salmon quants (replaces 05)

**Replaces:** `05_isoform_comparison.py`.

**Why:** proper isoform DE needs per-sample variance modeling. limma-voom on
tximport-imported Salmon quants (with `txOut=TRUE`,
`countsFromAbundance="scaledTPM"`, `ignoreAfterBar=TRUE`) is what the iEC repo
uses. It also gives you `logFC`, `P.Value`, `adj.P.Val` — the fields the
`_isoform_annotations.py` annotator expects.

**Source files to port:**
- `scripts/variants/16_isoform_de.R` — **take the VSMC branch only** (the `else if (cohort == "vsmc")` block at lines 131-151). The two contrasts it runs are exactly what this cohort should compute:
  - `vsmc_pooled_CAD_vs_nonCAD` — 2 CAD vs 5 non-CAD
  - `vsmc_C224Y_cluster_CAD_vs_ISO` — 2 CAD vs 1 ISO (C224Y cluster only)
  
  Strip the iEC branch entirely. Pull `qc_drop` from a simple hard-coded empty default (no QC drops needed here).
- `scripts/variants/26_vsmc_signature_interrogation.R` — gene-level DE + signature interrogation. **Optional for the first port — only worth doing if you also port the 26-gene unified signature CSV.** The plasma signature lives at `/sc/arion/work/kellej10/cadasil_elab_ipsc_ec/data/signature_unified_betas.csv`; if the manuscript needs the public repo to emit its own interrogation, port it. Otherwise skip.

**Key tximport call (lift verbatim):**
```r
txi <- tximport(quant_files, type = "salmon", txOut = TRUE,
                countsFromAbundance = "scaledTPM",
                ignoreAfterBar = TRUE)
```
The `ignoreAfterBar = TRUE` is load-bearing — the closed-loop Salmon index uses pipe-delimited gencode headers (`ENST00000456328.2|ENSG...|DDX11L2-202|...`), and without this flag the tx IDs in `rownames(txi$counts)` don't match `tx2gene.transcript_id`.

**tx2gene file:** reuse the one the iEC repo already built at `/sc/arion/work/kellej10/cadasil_elab_ipsc_ec/data/tx2gene_gencode_v44.csv` (cohort-independent). Just point to that absolute path, same as the iEC 16 script did for its `else if` branch.

**Acceptance test:** `results/de_reassigned/isoform/tx_de_vsmc_pooled_CAD_vs_nonCAD_panel.csv` exists with `logFC`, `P.Value`, `adj.P.Val` columns and rows for the FLT1/FLT4/FN1 panel transcripts. Sanity: sFLT1-i13 (`ENST00000541932`) should have `logFC ≈ +0.07, P ≈ 0.96` (non-replicating in VSMC — this is the headline finding from the cross-cohort work).

### Task 5 — Single-source isoform annotations

**Why:** the same FLT1/FLT4/FN1 labels should come out of both manuscripts. Today `03_build_matrices.py:65-77` has a `FN1_TRANSCRIPTS` dict and `04_build_export_package.py` / `05_isoform_comparison.py` use ad-hoc isoform names like `mFlt1_full_length` and `sFlt1` that don't match the iEC repo's `FLT1-201` / `FLT1-204` etc. In the iEC repo this was fixed by a shared module at `/sc/arion/work/kellej10/cadasil_elab_ipsc_ec/scripts/variants/_isoform_annotations.py`.

**Action:** copy that file here verbatim:
```
cp /sc/arion/work/kellej10/cadasil_elab_ipsc_ec/scripts/variants/_isoform_annotations.py \
   scripts/_isoform_annotations.py
```
Then update `03_build_matrices.py`, `04_build_export_package.py`, and the new `05_isoform_de.R`-derivatives to import `ANNOT` and `strip_version` from it (drop the local duplicates). A Python import snippet you can reuse:
```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _isoform_annotations import ANNOT, strip_version
```

**Acceptance test:** running the pipeline end-to-end emits CSVs where the `ensembl_name` column for a given `transcript_id` matches the iEC repo's outputs **exactly**. In particular `ENST00000541932 → FLT1-204`, `ENST00000432072 → FN1-209`, `ENST00000336916 → FN1-202`. These three were the genuine bugs in the iEC repo's first-pass table and are now fixed by the shared module — don't reintroduce them.

### Task 6 — Rename isoform output columns to match

The iEC export schema is:
```
transcript_id, ensembl_name, isoform_label, aa_length, n_exons,
biotype, notes, contrast, AveExpr, logFC, t, P.Value, adj.P.Val
```
(see `24_annotate_isoforms.py:FIELDS` in the iEC repo). Port `24_annotate_isoforms.py` here, wired against the new `05` output, and have it write `{FLT1,FLT4,FN1}_isoform_stats_annotated.csv` to `results/de_reassigned/supp_tables/`. The existing `04_build_export_package.py` bundler should then pick these up by filename.

### Task 7 — Write `PIPELINE.md`

The iEC repo's `PIPELINE.md` at `scripts/variants/PIPELINE.md` is the template. Mirror its structure:
- §0 Prerequisites
- §1 Sample-level alignment + counting
- §2 Joint SNP genotyping + donor clustering (the "correction" core)
- §3 Cohort QC
- §4 DE
- §5 Figures / supp tables
- §6 How to run (for this cohort it's literally "run the LSFs in numeric order")
- §7 FLT1 isoform nomenclature (**copy verbatim** from the iEC PIPELINE.md §7 — the sFLT1-i13 / sFLT1-i15 / FLT1-204 / FLT1-207 table is the disambiguation Jonah wants preserved)
- §8 File inventory

Overwrite the empty `RUN_PIPELINE.txt` with a one-line pointer to `PIPELINE.md` (or delete it).

### Task 8 — Rewrite `sample_metadata.csv`

The current `sample_metadata.csv` has a manually-authored `notch3_mutation_observed` column reflecting the Zhang errata. After Task 2 runs, the **output** of step 09 (`sample_metadata_reassigned.csv`) supersedes this file. Keep both:
- `sample_metadata.csv` = submitted labels (input only, for audit trail)
- `sample_metadata_reassigned.csv` = canonical downstream input (output of 09)

All scripts post-Task 2 should read `sample_metadata_reassigned.csv`, NOT `sample_metadata.csv`.

---

## 3. Do-not-port list

Things from the iEC repo that are deliberately NOT coming over:
- **`cohorts/` directory + `COHORT` env switch** — this repo is single-cohort; overhead not justified.
- **`QC_DROP` env var** — no samples need dropping here (already verified: all 7 have clean coverage).
- **`tier1` / `tier2_block` / `duplicateCorrelation` logic** — this cohort doesn't have multiple donors per condition (only XZ5+WZ16 share a donor, and that's the entire CAD arm). Pooled DE + C224Y-cluster DE is the whole story.
- **`13_de_tier1.R` / `14_de_tier2.R` / `15_signature_interrogation.R`** — these are iEC-tier-specific. Only 16 and 26 have VSMC-relevant branches.

---

## 4. Order of operations

Dependencies:
- Task 1 (git init) before any commit.
- Task 2 (SNP panel) before Task 4 (DE uses reassigned metadata).
- Task 3 (Salmon) before Task 4 (DE uses Salmon quants).
- Task 5 (shared annotations) can happen in parallel with Task 4, but Task 6 needs Task 5 done.
- Task 7 (docs) last.

Rough branch plan:
```
feat/git-init            -> Task 1
feat/snp-panel           -> Task 2
feat/salmon-closed-loop  -> Task 3
feat/limma-voom-de       -> Task 4 + 5 + 6 (these tie together)
docs/pipeline-md         -> Task 7
```

---

## 5. Verification: cross-repo reproducibility

Once all tasks land, running this repo end-to-end should produce results **identical** to the iEC repo's `vsmc_validation/` subtree (which re-ran the same 7 samples under the closed-loop pipeline). Spot-check:

| File in this repo | Same file in iEC | Should match? |
|---|---|---|
| `results/variants/panel_jaccard_matrix.csv` | `vsmc_validation/results/variants/panel_jaccard_matrix.csv` | yes, exactly |
| `data/sample_metadata_reassigned.csv` | `vsmc_validation/data/sample_metadata_reassigned.csv` | yes, exactly |
| `results/de_reassigned/isoform/tx_de_vsmc_pooled_CAD_vs_nonCAD_panel.csv` | `vsmc_validation/results/de_reassigned/isoform/tx_de_vsmc_pooled_CAD_vs_nonCAD_panel.csv` | yes, to numerical tolerance (same Salmon seed, same filterByExpr thresholds) |

If any of these diverge beyond numerical noise, something in the port is off — stop and debug before proceeding.

---

## 6. Things to ask Jonah before starting

If the follow-on agent hits any of these, pause and ask:

1. Should the fingerprint/Jaccard outputs keep their original filenames for backward compatibility with the existing `export/` bundle, or can file names be aligned with the iEC repo (`panel_jaccard_matrix.csv` vs the current `fingerprint_concordance.csv`)? Default: rename + keep a symlink for 1 release cycle.
2. The existing `E-MTAB-16303_export.zip` bundle is what gets shipped with the public-repo manuscript. After the port, regenerate and diff — any existing downstream consumer (figure script in the manuscript? supp table reviewer?) that depends on the old column names?
3. Does the public-repo manuscript need its own 26-gene signature interrogation (Task 4 Option B), or is the iEC repo the only place that table needs to live?

---

## 7. Commit hygiene

Use the same commit style as the iEC repo (imperative, explain the "why"):
```
Port NGSCheckMate SNP-panel donor genotyping

Replaces the housekeeping-het Jaccard fingerprint with the bcftools
mpileup/call flow used in cadasil_elab_ipsc_ec. Cluster-cut at
1-Jaccard_het=0.15 reproduces the XZ5 → C224Y errata independently.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
```

Open one PR per task where possible. Final PR body should link back to this
task document.
