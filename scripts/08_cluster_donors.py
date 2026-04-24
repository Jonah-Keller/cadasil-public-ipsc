#!/usr/bin/env python
# ============================================================================
# 08_cluster_donors.py -- infer donor identity from panel genotypes
# ----------------------------------------------------------------------------
# Ported from cadasil_elab_ipsc_ec. Consumes the multi-sample genotype TSV
# written by 07_genotype_snps.lsf. For each pair of samples:
#   - IBS = fraction of biallelic sites with identical genotype (both called)
#   - Jaccard_het = |het_A AND het_B| / |het_A OR het_B|
# Hierarchical cluster (complete linkage, 1-Jaccard_het) into donor groups.
#
# Writes:
#   results/variants/panel_ibs_matrix.csv
#   results/variants/panel_jaccard_matrix.csv
#   results/variants/panel_call_counts.csv
#   results/variants/panel_donor_clusters.csv
#   results/variants/panel_dendrogram_merges.csv
# ============================================================================
import os, sys, csv
from collections import defaultdict
from pathlib import Path

PROJ_ROOT = Path("/sc/arion/work/kellej10/cadasil_public_ipsc")
RESULTS_DIR = PROJ_ROOT / "results/variants"
GT_TSV = Path(os.environ.get(
    "PANEL_GT_TSV",
    "/sc/arion/scratch/kellej10/cadasil_public_ipsc/genotype/vsmc_panel_genotypes.tsv"))

MIN_DP = 10
# Within-donor Jaccard_het is observed 0.95-0.98; across-donor 0.28-0.38 in
# the iEC cohort. Cut at distance = 1 - 0.85 = 0.15 on Jaccard_het gives
# clean donor groups. (IBS inflates across all samples due to monomorphic
# sites, regardless of donor.)
JACCARD_DIST_CUT = 0.15


def parse_gt_field(field: str):
    """Parse 'GT:DP:AD' -> (gt_code, dp).

    gt_code: 0 = hom_ref, 1 = het, 2 = hom_alt, -1 = missing/low_dp.
    """
    if field in (".", "./.", ".|.", ".:0:0,0", ""):
        return -1, 0
    parts = field.split(":")
    gt = parts[0]
    dp = 0
    for p in parts[1:]:
        if p.isdigit():
            dp = int(p); break
    if dp < MIN_DP:
        return -1, dp
    g = gt.replace("|", "/")
    if g in ("0/0",): return 0, dp
    if g in ("1/1",): return 2, dp
    if g in ("0/1", "1/0"): return 1, dp
    return -1, dp


def main():
    if not GT_TSV.exists():
        sys.exit(f"[08] ERROR: genotype TSV not found: {GT_TSV}")

    print(f"[08] reading {GT_TSV}")
    with open(GT_TSV) as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header[4:]
        n_samp = len(samples)
        gt_mat = [[] for _ in range(n_samp)]
        n_sites = 0
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 4 + n_samp:
                continue
            fields = parts[4:]
            site_codes = [parse_gt_field(fl)[0] for fl in fields]
            if sum(1 for c in site_codes if c >= 0) < 3:
                continue
            for i, c in enumerate(site_codes):
                gt_mat[i].append(c)
            n_sites += 1
    print(f"[08] kept {n_sites} panel sites across {n_samp} samples")

    call_totals = []
    for i, s in enumerate(samples):
        gts = gt_mat[i]
        ncall = sum(1 for c in gts if c >= 0)
        nhet = sum(1 for c in gts if c == 1)
        call_totals.append((s, ncall, nhet, n_sites))

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_DIR / "panel_call_counts.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample", "n_called", "n_het", "n_sites_total",
                    "call_rate", "het_rate_of_called"])
        for s, nc, nh, nt in call_totals:
            cr = f"{nc/nt:.4f}" if nt else "0"
            hr = f"{nh/nc:.4f}" if nc else "0"
            w.writerow([s, nc, nh, nt, cr, hr])

    # Pairwise IBS + Jaccard_het
    ibs = [[1.0]*n_samp for _ in range(n_samp)]
    jac = [[1.0]*n_samp for _ in range(n_samp)]
    for i in range(n_samp):
        gi = gt_mat[i]
        for j in range(i+1, n_samp):
            gj = gt_mat[j]
            same = 0; both = 0
            hi_or_j = 0; hi_and_j = 0
            for a, b in zip(gi, gj):
                if a < 0 or b < 0: continue
                both += 1
                if a == b: same += 1
                het_a = (a == 1); het_b = (b == 1)
                if het_a or het_b: hi_or_j += 1
                if het_a and het_b: hi_and_j += 1
            ibs_v = (same / both) if both else 0.0
            jac_v = (hi_and_j / hi_or_j) if hi_or_j else 0.0
            ibs[i][j] = ibs[j][i] = ibs_v
            jac[i][j] = jac[j][i] = jac_v

    def write_matrix(path, mat):
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["sample"] + samples)
            for i, s in enumerate(samples):
                w.writerow([s] + [f"{v:.4f}" for v in mat[i]])

    write_matrix(RESULTS_DIR / "panel_ibs_matrix.csv", ibs)
    write_matrix(RESULTS_DIR / "panel_jaccard_matrix.csv", jac)

    # Hierarchical cluster on 1 - Jaccard_het, complete linkage
    clusters = {i: [i] for i in range(n_samp)}

    def cluster_dist_complete(a_idx, b_idx):
        worst = 0.0
        for x in clusters[a_idx]:
            for y in clusters[b_idx]:
                d = 1.0 - jac[x][y]
                if d > worst: worst = d
        return worst

    active = set(clusters.keys())
    merges = []
    while len(active) > 1:
        best_d = 9e9; best_pair = None
        alist = sorted(active)
        for ia in range(len(alist)):
            for ib in range(ia+1, len(alist)):
                d = cluster_dist_complete(alist[ia], alist[ib])
                if d < best_d:
                    best_d = d; best_pair = (alist[ia], alist[ib])
        a, b = best_pair
        merges.append((a, b, best_d, len(clusters[a])+len(clusters[b])))
        if best_d > JACCARD_DIST_CUT:
            break
        clusters[a] = clusters[a] + clusters[b]
        active.remove(b)
        del clusters[b]

    cluster_id = {}
    for cid, (root, members) in enumerate(sorted(clusters.items()), start=1):
        for m in members:
            cluster_id[samples[m]] = cid

    # Annotate clusters with NOTCH3 variant evidence
    variant_csv = RESULTS_DIR / "notch3_variant_summary.csv"
    sample_variant = defaultdict(dict)
    if variant_csv.exists():
        with open(variant_csv) as f:
            for row in csv.DictReader(f):
                try:
                    vaf = float(row["vaf"]) if row["vaf"] else 0.0
                except ValueError:
                    vaf = 0.0
                sample_variant[row["sample"]][row["variant"]] = (vaf, int(row["depth"] or 0))

    cluster_summary = defaultdict(lambda: {"R153C": [], "C224Y": [], "members": []})
    for s in samples:
        cid = cluster_id[s]
        cluster_summary[cid]["members"].append(s)
        for v in ("p.R153C", "p.C224Y"):
            vaf, dp = sample_variant.get(s, {}).get(v, (0.0, 0))
            if vaf >= 0.15 and dp >= 10:
                cluster_summary[cid][v.replace("p.", "")].append(f"{s}(VAF={vaf:.2f},dp={dp})")

    # Sample -> cluster assignment (joined with submitted metadata)
    meta_csv = PROJ_ROOT / "sample_metadata.csv"
    meta = {}
    if meta_csv.exists():
        with open(meta_csv) as f:
            for row in csv.DictReader(f):
                # This repo's metadata keys on sample_name; alias via err_id too
                meta[row.get("sample_name", row.get("sample_id", ""))] = row

    with open(RESULTS_DIR / "panel_donor_clusters.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample", "cluster_id", "meta_donor", "meta_genotype",
                    "meta_protocol", "n_called", "n_het"])
        call_lookup = {s: (nc, nh) for s, nc, nh, _ in call_totals}
        for s in samples:
            nc, nh = call_lookup.get(s, (0, 0))
            m = meta.get(s, {})
            # Derive donor from submitted NOTCH3 mutation (fallback to WT)
            submitted_mut = m.get("notch3_mutation_reported", "none")
            meta_donor = submitted_mut.replace("p.", "") if submitted_mut not in ("none", "", "corrected_R153C", "corrected_C224Y") \
                         else (submitted_mut.replace("corrected_", "").replace("p.", "") if submitted_mut.startswith("corrected_") \
                         else "WT")
            meta_geno = m.get("genotype", "")
            w.writerow([s, cluster_id[s], meta_donor, meta_geno,
                        m.get("sample_type", ""), nc, nh])

    print("\n[08] cluster summary")
    for cid in sorted(cluster_summary.keys()):
        d = cluster_summary[cid]
        print(f"  cluster {cid}: {len(d['members'])} samples")
        print(f"    members: {', '.join(d['members'])}")
        if d["R153C"]: print(f"    R153C hits: {d['R153C']}")
        if d["C224Y"]: print(f"    C224Y hits: {d['C224Y']}")

    with open(RESULTS_DIR / "panel_dendrogram_merges.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["step", "clusterA_root", "clusterB_root",
                    "distance_1minusJaccardHet", "merged_size"])
        for step, (a, b, d, sz) in enumerate(merges, 1):
            w.writerow([step, samples[a], samples[b], f"{d:.4f}", sz])

    print(f"[08] done. outputs in {RESULTS_DIR}")


if __name__ == "__main__":
    main()
