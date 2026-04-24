#!/usr/bin/env Rscript
# 26_vsmc_signature_interrogation.R -- VSMC gene DE + unified plasma-signature
# interrogation. Writes to results/de_reassigned/.
#
# Source of plasma signature (26 genes, effect priors from iEC cohort) lives
# in data/signature_unified_betas.csv (symlink to shared iEC file).
#
# Counts: this repo's existing gene matrix (export/expression/raw_counts.csv)
# already has gene_id + gene_name columns, so we read directly from there.

Sys.setenv(OMP_NUM_THREADS = "4", OPENBLAS_NUM_THREADS = "4")
set.seed(42)
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma)
})

root   <- Sys.getenv("PROJ_ROOT", "/sc/arion/work/kellej10/cadasil_public_ipsc")
sig_f  <- file.path(root, "data/signature_unified_betas.csv")
meta_f <- file.path(root, "data/sample_metadata_reassigned.csv")
cnt_f  <- file.path(root, "export/expression/raw_counts.csv")
outdir <- file.path(root, "results/de_reassigned")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- fread(meta_f)
cnt  <- fread(cnt_f)
setnames(cnt, "ensembl_id", "gene_id")  # align with iEC schema
sig  <- fread(sig_f)
sig_genes <- sig[, .(gene, direction_plasma,
                     beta_prior_iec = beta_iec_primary,
                     p_prior_iec    = p_iec_primary)]

gene_ids <- cnt$gene_id
gene_nm  <- cnt$gene_name
M <- as.matrix(cnt[, !c("gene_id","gene_name")])
rownames(M) <- gene_ids
M <- M[, meta$sample_id, drop = FALSE]

fit_de <- function(samples, name) {
  m <- meta[sample_id %in% samples]
  m[, condition := factor(ifelse(genotype_inferred == "CAD","CAD","WT"),
                          levels = c("WT","CAD"))]
  cm  <- M[, m$sample_id, drop = FALSE]
  dge <- DGEList(counts = cm, genes = data.frame(gene_id = rownames(cm),
                                                 gene_name = gene_nm,
                                                 stringsAsFactors = FALSE))
  keep <- filterByExpr(dge, group = m$condition,
                       min.count = 5, min.total.count = 10)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  design <- stats::model.matrix(~ condition, data = m)
  v   <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  tt  <- topTable(fit, coef = "conditionCAD", number = Inf, sort.by = "none")
  tt$gene_id   <- dge$genes$gene_id
  tt$gene_name <- dge$genes$gene_name
  tt$contrast  <- name
  fwrite(tt, file.path(outdir, sprintf("de_%s.csv", name)))
  cat(sprintf("  [%s] %d genes tested; nominal p<0.05: %d; FDR<0.05: %d\n",
              name, nrow(tt),
              sum(tt$P.Value < 0.05, na.rm=TRUE),
              sum(tt$adj.P.Val < 0.05, na.rm=TRUE)))
  invisible(tt)
}

cat("\n=== Gene DE: VSMC pooled CAD vs non-CAD ===\n")
tt_pool <- fit_de(meta$sample_id, "vsmc_pooled_CAD_vs_nonCAD")

cat("\n=== Gene DE: VSMC C224Y cluster CAD vs ISO ===\n")
c224_ids <- meta[donor_inferred == "C224Y", sample_id]
tt_iso   <- fit_de(c224_ids, "vsmc_C224Y_cluster_CAD_vs_ISO")

de_tables <- list(VSMC_pool  = as.data.table(tt_pool),
                  VSMC_C224Y = as.data.table(tt_iso))
short     <- names(de_tables)

out <- sig_genes
for (s in short) {
  sub <- de_tables[[s]][gene_name %in% sig_genes$gene,
                        .(gene = gene_name, logFC, P.Value, adj.P.Val)]
  setnames(sub, c("logFC","P.Value","adj.P.Val"),
                c(sprintf("lfc_%s", s), sprintf("p_%s", s), sprintf("fdr_%s", s)))
  out <- merge(out, sub, by = "gene", all.x = TRUE, sort = FALSE)
}
out[, expected_sign := ifelse(direction_plasma == "Up",   1,
                      ifelse(direction_plasma == "Down", -1, NA_integer_))]
for (s in short) {
  lfc <- sprintf("lfc_%s", s)
  p   <- sprintf("p_%s",   s)
  fdr <- sprintf("fdr_%s", s)
  out[, paste0("dir_match_", s) := sign(get(lfc)) == expected_sign]
  out[, paste0("nom_",       s) := get(p)   < 0.05]
  out[, paste0("sig_",       s) := get(fdr) < 0.05]
}

setcolorder(out, c("gene","direction_plasma","expected_sign",
                   "beta_prior_iec","p_prior_iec"))
fwrite(out, file.path(outdir, "signature_interrogation_vsmc.csv"))
cat(sprintf("\nwrote %s\n", file.path(outdir, "signature_interrogation_vsmc.csv")))

cat("\n=== VSMC per-tier signature performance ===\n")
summ <- data.table(
  tier = short,
  genes_tested = sapply(short, function(s)
                   sum(!is.na(out[[sprintf("lfc_%s", s)]]))),
  direction_match = sapply(short, function(s)
                     sum(out[[sprintf("dir_match_%s", s)]], na.rm = TRUE)),
  nominal_p05 = sapply(short, function(s)
                 sum(out[[sprintf("nom_%s", s)]], na.rm = TRUE)),
  nominal_p05_concordant = sapply(short, function(s)
                 sum(out[[sprintf("nom_%s", s)]] &
                     out[[sprintf("dir_match_%s", s)]], na.rm = TRUE)),
  fdr05 = sapply(short, function(s)
           sum(out[[sprintf("sig_%s", s)]], na.rm = TRUE)),
  fdr05_concordant = sapply(short, function(s)
           sum(out[[sprintf("sig_%s", s)]] &
               out[[sprintf("dir_match_%s", s)]], na.rm = TRUE))
)
fwrite(summ, file.path(outdir, "signature_interrogation_vsmc_summary.csv"))
print(summ)

cat("\n=== Binomial sign-concordance test vs plasma direction ===\n")
for (s in short) {
  n_tested <- sum(!is.na(out[[sprintf("lfc_%s", s)]]) &
                  !is.na(out$expected_sign))
  n_match  <- sum(out[[sprintf("dir_match_%s", s)]], na.rm = TRUE)
  if (n_tested == 0) next
  bt <- binom.test(n_match, n_tested, p = 0.5, alternative = "greater")
  cat(sprintf("  %-12s  %d/%d concordant (%.0f%%)  p=%.4f\n",
              s, n_match, n_tested, 100*n_match/n_tested, bt$p.value))
}

pretty <- out[, .(
  gene, direction = direction_plasma,
  prior_iec = signif(beta_prior_iec, 3),
  VSMC_pool   = sprintf("%+.2f %s", lfc_VSMC_pool,
                        ifelse(sig_VSMC_pool,"***",
                        ifelse(nom_VSMC_pool,"*",""))),
  VSMC_C224Y  = sprintf("%+.2f %s", lfc_VSMC_C224Y,
                        ifelse(sig_VSMC_C224Y,"***",
                        ifelse(nom_VSMC_C224Y,"*",""))),
  dir_hits = rowSums(cbind(dir_match_VSMC_pool, dir_match_VSMC_C224Y),
                     na.rm = TRUE),
  nom_hits = rowSums(cbind(nom_VSMC_pool   & dir_match_VSMC_pool,
                           nom_VSMC_C224Y  & dir_match_VSMC_C224Y),
                     na.rm = TRUE)
)]
setorder(pretty, -nom_hits, -dir_hits)
fwrite(pretty, file.path(outdir, "signature_interrogation_vsmc_pretty.csv"))
cat("\n=== VSMC per-gene interrogation (* = p<0.05, *** = FDR<0.05) ===\n")
print(pretty, nrows = Inf)

cat("\ndone\n")
