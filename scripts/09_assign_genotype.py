#!/usr/bin/env python
# ============================================================================
# 09_assign_genotype.py -- reassign donor / genotype labels empirically
# ----------------------------------------------------------------------------
# Ported from cadasil_elab_ipsc_ec. Cohort-specific adjustments:
#   - COHORT_DONOR_VARIANTS = R153C C224Y (no R544C/C87R/G420C in VSMC cohort)
#   - No KOLF2.1J reference set (WZ12/WZ4/WZ8 are independent WTs)
#
# Inputs:
#   results/variants/panel_donor_clusters.csv   (from 08_cluster_donors.py)
#   results/variants/notch3_variant_summary.csv (from 03_build_matrices.py)
#
# Rules:
#   1. Each cluster from 08 is treated as a single donor line.
#   2. For each cluster, pick the CAD driver variant with >=1 sample at
#      VAF >= 0.15 and DP >= 5. Tie-break by highest observed VAF.
#   3. Within each CADASIL donor cluster:
#        VAF >= 0.15 (DP>=5)  -> CAD
#        VAF <= 0.05 (DP>=5)  -> ISO (isogenic control, same donor)
#        0.05 < VAF < 0.15    -> AMBIGUOUS
#        DP < 5               -> LOW_COV
#   4. Pure-WT cluster (no driver detected): donor = <submitted>_noVariant
#      (or SINGLETON_<submitted> for single-member clusters), genotype = WT.
#
# Output: data/sample_metadata_reassigned.csv
# ============================================================================
import os, csv
from collections import defaultdict
from pathlib import Path

PROJ_ROOT = Path("/sc/arion/work/kellej10/cadasil_public_ipsc")
RESULTS = PROJ_ROOT / "results/variants"
DATA = PROJ_ROOT / "data"
DATA.mkdir(parents=True, exist_ok=True)

COHORT_VARIANTS = os.environ.get("COHORT_DONOR_VARIANTS",
                                 "R153C C224Y").split()

CAD_VAF_MIN = 0.15
ISO_VAF_MAX = 0.05
MIN_DEPTH = 5


def load_variant_vafs():
    out = defaultdict(dict)
    with open(RESULTS / "notch3_variant_summary.csv") as f:
        for r in csv.DictReader(f):
            try:
                vaf = float(r["vaf"]) if r["vaf"] else 0.0
            except ValueError:
                vaf = 0.0
            try:
                dp = int(r["depth"]) if r["depth"] else 0
            except ValueError:
                dp = 0
            out[r["sample"]][r["variant"]] = (vaf, dp)
    return out


def load_clusters():
    with open(RESULTS / "panel_donor_clusters.csv") as f:
        return list(csv.DictReader(f))


def load_meta():
    """Load submitted metadata keyed by sample_name (= BAM-filename prefix)."""
    out = {}
    with open(PROJ_ROOT / "sample_metadata.csv") as f:
        for r in csv.DictReader(f):
            key = r.get("sample_name") or r.get("sample_id")
            if key:
                out[key] = r
    return out


def cluster_donor_call(samples, vafs):
    best = {}
    for v_name in COHORT_VARIANTS:
        v_key = f"p.{v_name}"
        maxv = 0.0; n_hits = 0
        for s in samples:
            vaf, dp = vafs.get(s, {}).get(v_key, (0.0, 0))
            if dp < MIN_DEPTH: continue
            if vaf >= CAD_VAF_MIN: n_hits += 1
            if vaf > maxv: maxv = vaf
        best[v_name] = (maxv, n_hits)

    hits = [(n, mv, v) for v, (mv, n) in best.items() if n >= 1]
    if hits:
        hits.sort(key=lambda x: (-x[0], -x[1]))
        n, mv, v = hits[0]
        return v, f"p.{v}"
    return None, None


def classify_sample(sample, driver_variant, vafs):
    if driver_variant is None:
        return "WT", "", ""
    vaf, dp = vafs.get(sample, {}).get(driver_variant, (0.0, 0))
    if dp < MIN_DEPTH:
        return "LOW_COV", f"{vaf:.3f}", dp
    if vaf >= CAD_VAF_MIN:
        return "CAD", f"{vaf:.3f}", dp
    if vaf <= ISO_VAF_MAX:
        return "ISO", f"{vaf:.3f}", dp
    return "AMBIGUOUS", f"{vaf:.3f}", dp


def meta_donor_for(m):
    """Derive submitted donor label from sample_metadata.csv row."""
    submitted = m.get("notch3_mutation_reported", "none")
    if submitted in ("", "none", None):
        return "WT"
    if submitted.startswith("corrected_"):
        return submitted.replace("corrected_", "")
    return submitted.replace("p.", "")


def main():
    vafs = load_variant_vafs()
    cluster_rows = load_clusters()
    meta = load_meta()

    by_cluster = defaultdict(list)
    for r in cluster_rows:
        by_cluster[int(r["cluster_id"])].append(r["sample"])

    cluster_donor = {}
    for cid, members in by_cluster.items():
        donor_label, driver = cluster_donor_call(members, vafs)
        cluster_donor[cid] = (donor_label, driver)

    # Pure-WT clusters: no KOLF2.1J reference set in this cohort (WZ12/WZ4/WZ8
    # are independent WTs). Label by majority metadata donor + _noVariant, or
    # SINGLETON_<donor> for single-member clusters.
    for cid, (donor_label, driver) in list(cluster_donor.items()):
        if donor_label is not None:
            continue
        members = by_cluster[cid]
        meta_donors = [meta_donor_for(meta.get(s, {})) for s in members]
        meta_donors = [d for d in meta_donors if d]
        if meta_donors:
            counts = defaultdict(int)
            for d in meta_donors:
                counts[d] += 1
            best = max(counts.items(), key=lambda x: x[1])[0]
            if len(members) == 1:
                cluster_donor[cid] = (f"SINGLETON_{best}", None)
            else:
                cluster_donor[cid] = (f"{best}_noVariant", None)
        else:
            cluster_donor[cid] = ("UNKNOWN_WT", None)

    out_rows = []
    for r in cluster_rows:
        s = r["sample"]
        cid = int(r["cluster_id"])
        donor_inferred, driver = cluster_donor[cid]
        genotype, vaf_str, dp = classify_sample(s, driver, vafs)
        m = meta.get(s, {})
        meta_donor = meta_donor_for(m)
        # meta genotype: normalize CAD/ISO/WT
        raw_mg = m.get("genotype", "") or m.get("sample_type", "")
        if "CAD" in raw_mg.upper(): meta_geno = "CAD"
        elif "iso" in raw_mg.lower() or "mutation_corrected" in raw_mg: meta_geno = "ISO"
        else: meta_geno = "WT"

        mismatch_donor = donor_inferred != meta_donor
        mismatch_geno = genotype != meta_geno
        out_rows.append({
            "sample_id": s,
            "err_id": m.get("err_id", ""),
            "cluster_id": cid,
            "donor_inferred": donor_inferred,
            "genotype_inferred": genotype,
            "driver_variant": driver or "",
            "driver_vaf": vaf_str,
            "driver_depth": dp,
            "meta_donor": meta_donor,
            "meta_genotype": meta_geno,
            "donor_mismatch": int(mismatch_donor),
            "genotype_mismatch": int(mismatch_geno),
            "sample_type": m.get("sample_type", ""),
            "cell_type": m.get("cell_type", ""),
            "n_called": r["n_called"],
            "n_het": r["n_het"],
        })

    out_path = DATA / "sample_metadata_reassigned.csv"
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)

    print(f"[09] wrote {out_path}")
    print("\n[09] reassignment summary")
    donor_counts = defaultdict(lambda: defaultdict(int))
    for r in out_rows:
        donor_counts[r["donor_inferred"]][r["genotype_inferred"]] += 1
    print(f"  {'donor':<20} {'CAD':>5} {'ISO':>5} {'WT':>5} {'AMB':>5} {'LOW':>5}")
    for d, gcounts in sorted(donor_counts.items()):
        print(f"  {d:<20} {gcounts['CAD']:>5} {gcounts['ISO']:>5} "
              f"{gcounts['WT']:>5} {gcounts['AMBIGUOUS']:>5} {gcounts['LOW_COV']:>5}")

    mm = [r for r in out_rows if r["donor_mismatch"] or r["genotype_mismatch"]]
    print(f"\n[09] metadata mismatches: {len(mm)}")
    for r in mm:
        print(f"  {r['sample_id']}: meta=({r['meta_donor']}, {r['meta_genotype']}) "
              f"-> inferred=({r['donor_inferred']}, {r['genotype_inferred']}) "
              f"driver={r['driver_variant']} VAF={r['driver_vaf']} dp={r['driver_depth']}")

    # Human-readable summary
    md = RESULTS / "cohort_reassigned_summary.md"
    by_donor = defaultdict(list)
    for r in out_rows:
        by_donor[r["donor_inferred"]].append(r)
    with open(md, "w") as f:
        f.write("# VSMC cohort (E-MTAB-16303) -- empirical donor / genotype reassignment\n\n")
        f.write("Donor identity inferred from NGSCheckMate common-SNP genotypes "
                "(pairwise Jaccard_het, complete linkage, cut at 1-J = 0.15). "
                "Genotype (CAD/ISO/WT) inferred from NOTCH3 variant pileup at "
                "R153C / C224Y sites.\n\n")
        f.write("## Per-donor sample counts\n\n")
        f.write("| donor_inferred | CAD | ISO | WT | AMBIGUOUS | LOW_COV | total |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for d, rows in sorted(by_donor.items()):
            gc = defaultdict(int)
            for r in rows: gc[r["genotype_inferred"]] += 1
            f.write(f"| {d} | {gc['CAD']} | {gc['ISO']} | {gc['WT']} "
                    f"| {gc['AMBIGUOUS']} | {gc['LOW_COV']} | {len(rows)} |\n")
        f.write("\n## Per-cluster detail\n\n")
        for d, rows in sorted(by_donor.items()):
            cids = sorted(set(r['cluster_id'] for r in rows))
            f.write(f"### {d} (cluster{'s' if len(cids)>1 else ''}: {cids})\n\n")
            f.write("| sample | inferred genotype | driver VAF (dp) | meta donor | meta geno |\n")
            f.write("|---|---|---|---|---|\n")
            for r in rows:
                vaf_cell = f"{r['driver_vaf']} (dp={r['driver_depth']})" if r['driver_variant'] else "---"
                f.write(f"| {r['sample_id']} | {r['genotype_inferred']} | {vaf_cell} "
                        f"| {r['meta_donor']} | {r['meta_genotype']} |\n")
            f.write("\n")
        f.write("## Metadata mismatches\n\n")
        f.write(f"Total: {len(mm)} / {len(out_rows)} samples differ from submitted metadata.\n\n")
        f.write("| sample | meta donor | meta geno | inferred donor | inferred geno | driver | VAF | dp |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        for r in mm:
            f.write(f"| {r['sample_id']} | {r['meta_donor']} | {r['meta_genotype']} "
                    f"| {r['donor_inferred']} | {r['genotype_inferred']} "
                    f"| {r['driver_variant']} | {r['driver_vaf']} | {r['driver_depth']} |\n")
    print(f"[09] wrote {md}")


if __name__ == "__main__":
    main()
