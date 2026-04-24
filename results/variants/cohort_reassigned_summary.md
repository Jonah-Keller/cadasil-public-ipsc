# VSMC cohort (E-MTAB-16303) -- empirical donor / genotype reassignment

Donor identity inferred from NGSCheckMate common-SNP genotypes (pairwise Jaccard_het, complete linkage, cut at 1-J = 0.15). Genotype (CAD/ISO/WT) inferred from NOTCH3 variant pileup at R153C / C224Y sites.

## Per-donor sample counts

| donor_inferred | CAD | ISO | WT | AMBIGUOUS | LOW_COV | total |
|---|---|---|---|---|---|---|
| C224Y | 2 | 1 | 0 | 0 | 0 | 3 |
| SINGLETON_WT | 0 | 0 | 2 | 0 | 0 | 2 |
| WT_noVariant | 0 | 0 | 2 | 0 | 0 | 2 |

## Per-cluster detail

### C224Y (cluster: [2])

| sample | inferred genotype | driver VAF (dp) | meta donor | meta geno |
|---|---|---|---|---|
| WZ16_CADASIL_C224Y | CAD | 0.454 (dp=251) | C224Y | WT |
| XZ5_CADASIL_R153C | CAD | 0.416 (dp=464) | R153C | WT |
| XZ6_CADASIL_R153C_isoCtrl | ISO | 0.000 (dp=72) | R153C | ISO |

### SINGLETON_WT (clusters: [1, 4])

| sample | inferred genotype | driver VAF (dp) | meta donor | meta geno |
|---|---|---|---|---|
| WZ12_CADASIL_Control | WT | --- | WT | WT |
| WZ8_CADASIL_Control | WT | --- | WT | WT |

### WT_noVariant (cluster: [3])

| sample | inferred genotype | driver VAF (dp) | meta donor | meta geno |
|---|---|---|---|---|
| WZ4_CADASIL_Control | WT | --- | WT | WT |
| XZ7_CADASIL_C224Y_isoCtrl | WT | --- | C224Y | ISO |

## Metadata mismatches

Total: 7 / 7 samples differ from submitted metadata.

| sample | meta donor | meta geno | inferred donor | inferred geno | driver | VAF | dp |
|---|---|---|---|---|---|---|---|
| WZ12_CADASIL_Control | WT | WT | SINGLETON_WT | WT |  |  |  |
| WZ16_CADASIL_C224Y | C224Y | WT | C224Y | CAD | p.C224Y | 0.454 | 251 |
| WZ4_CADASIL_Control | WT | WT | WT_noVariant | WT |  |  |  |
| WZ8_CADASIL_Control | WT | WT | SINGLETON_WT | WT |  |  |  |
| XZ5_CADASIL_R153C | R153C | WT | C224Y | CAD | p.C224Y | 0.416 | 464 |
| XZ6_CADASIL_R153C_isoCtrl | R153C | ISO | C224Y | ISO | p.C224Y | 0.000 | 72 |
| XZ7_CADASIL_C224Y_isoCtrl | C224Y | ISO | WT_noVariant | WT |  |  |  |
