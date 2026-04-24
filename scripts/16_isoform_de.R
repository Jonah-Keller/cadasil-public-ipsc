#!/usr/bin/env Rscript
# 16_isoform_de.R -- transcript-level limma-voom DE for FN1, FLT4, FLT1
# using Salmon quants (closed-loop from FASTQ).
#
# Replaces 05_isoform_comparison.py (featureCounts transcript-level CPM
# ratios). featureCounts distributes ambiguous reads equally across
# overlapping transcripts, so isoform counts are not independent; Salmon's
# EM assigns reads probabilistically, giving statistically tractable
# isoform-level DE via tximport(txOut=TRUE) + limma-voom.
#
# VSMC cohort: 7 samples, donor-condition confounding is absolute (each
# donor cluster is all-CAD or all-ISO/WT), so no donor blocking. The
# C224Y-cluster contrast is gated on minimum sample counts.
#
# Inputs:
#   data/sample_metadata_reassigned.csv  (from 09_assign_genotype.py)
#   data/salmon_quant/<sample>/quant.sf  (from Salmon step; symlinked or fresh)
#   data/tx2gene_gencode_v44.csv
# Output: results/de_reassigned/isoform/tx_de_<contrast>_{all,panel}.csv
#                                       tx_logcpm_<contrast>.csv

Sys.setenv(OMP_NUM_THREADS = "4", OPENBLAS_NUM_THREADS = "4")
set.seed(42)
suppressPackageStartupMessages({
  library(data.table); library(tximport); library(edgeR); library(limma)
})

root   <- Sys.getenv("PROJ_ROOT", "/sc/arion/work/kellej10/cadasil_public_ipsc")
outdir <- file.path(root, "results/de_reassigned/isoform")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- fread(file.path(root, "data/sample_metadata_reassigned.csv"))
qc_drop <- strsplit(Sys.getenv("QC_DROP", ""), "\\s+")[[1]]
qc_drop <- qc_drop[nzchar(qc_drop)]
meta <- meta[!sample_id %in% qc_drop]

tx2g <- fread(file.path(root, "data/tx2gene_gencode_v44.csv"))

quant_files <- file.path(root, "data/salmon_quant", meta$sample_id, "quant.sf")
names(quant_files) <- meta$sample_id
stopifnot(all(file.exists(quant_files)))
txi <- tximport(quant_files, type = "salmon", txOut = TRUE,
                countsFromAbundance = "scaledTPM",
                ignoreAfterBar = TRUE)
cnts <- round(txi$counts)
storage.mode(cnts) <- "integer"

tx_info <- tx2g[transcript_id %in% rownames(cnts)]
setkey(tx_info, transcript_id)
panel_tx <- tx_info[gene_name %in% c("FN1","FLT4","FLT1"), transcript_id]
cat(sprintf("Panel transcripts in tx2gene: %d\n", length(panel_tx)))

run_tx_de <- function(m, name, design_rhs = "~ condition") {
  cm  <- cnts[, m$sample_id, drop = FALSE]
  dge <- DGEList(counts = cm)
  keep <- filterByExpr(dge, group = m$condition,
                       min.count = 5, min.total.count = 10)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  design <- stats::model.matrix(as.formula(design_rhs), data = m)
  v   <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  stopifnot("conditionCAD" %in% colnames(design))
  tt  <- topTable(fit, coef = "conditionCAD", number = Inf, sort.by = "none")
  tt$transcript_id <- rownames(tt)
  ti  <- tx_info[match(tt$transcript_id, transcript_id)]
  tt$gene_name       <- ti$gene_name
  tt$transcript_type <- ti$transcript_type
  tt$contrast        <- name
  fwrite(tt, file.path(outdir, sprintf("tx_de_%s_all.csv", name)))
  panel <- tt[tt$gene_name %in% c("FN1","FLT4","FLT1"), ]
  panel <- panel[order(panel$gene_name, panel$P.Value), ]
  fwrite(panel, file.path(outdir, sprintf("tx_de_%s_panel.csv", name)))

  ids <- intersect(rownames(v$E), panel_tx)
  e  <- v$E[ids, , drop = FALSE]
  long <- data.table(
    sample_id     = rep(colnames(e), each = nrow(e)),
    transcript_id = rep(rownames(e), times = ncol(e)),
    log2cpm       = as.vector(e))
  long[, gene_name       := tx_info[long$transcript_id, gene_name,       on = "transcript_id"]]
  long[, transcript_type := tx_info[long$transcript_id, transcript_type, on = "transcript_id"]]
  long <- merge(long, m[, .(sample_id, condition, donor_inferred)],
                by = "sample_id", all.x = TRUE)
  long[, contrast := name]
  fwrite(long, file.path(outdir, sprintf("tx_logcpm_%s.csv", name)))

  cat(sprintf("  [%s] panel hits nominal: %d / %d;  FDR<0.05: %d\n",
              name,
              sum(panel$P.Value < 0.05, na.rm=TRUE),
              nrow(panel),
              sum(panel$adj.P.Val < 0.05, na.rm=TRUE)))
  print(head(panel[, c("gene_name","transcript_id","transcript_type",
                       "logFC","P.Value","adj.P.Val")], 12))
  invisible(panel)
}

# Contrasts: pooled CAD vs non-CAD, and C224Y-cluster CAD vs ISO.
contrasts_list <- list(

  vsmc_pooled_CAD_vs_nonCAD = list(
    select = function(meta) {
      m <- copy(meta)
      m[, condition := factor(
          ifelse(genotype_inferred == "CAD", "CAD", "WT"),
          levels = c("WT", "CAD"))]
      m
    },
    design_rhs = "~ condition"
  ),

  vsmc_C224Y_cluster_CAD_vs_ISO = list(
    select = function(meta) {
      m <- meta[donor_inferred == "C224Y"]
      if (nrow(m) < 3 ||
          sum(m$genotype_inferred == "CAD") < 1 ||
          sum(m$genotype_inferred == "ISO") < 1) {
        return(NULL)
      }
      m[, condition := factor(
          ifelse(genotype_inferred == "CAD", "CAD", "WT"),
          levels = c("WT", "CAD"))]
      m
    },
    design_rhs = "~ condition"
  )
)

for (nm in names(contrasts_list)) {
  spec <- contrasts_list[[nm]]
  m    <- spec$select(copy(meta))
  if (is.null(m) || nrow(m) == 0L) {
    cat(sprintf("\n[%s] skipped (selector returned empty / gate failed)\n", nm))
    next
  }
  cat(sprintf("\n=== Isoform DE: %s (n=%d) ===\n", nm, nrow(m)))
  run_tx_de(m, nm, design_rhs = spec$design_rhs)
}

cat(sprintf("\nIsoform DE outputs -> %s\n", outdir))
