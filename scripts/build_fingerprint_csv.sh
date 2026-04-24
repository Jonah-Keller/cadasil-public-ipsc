#!/bin/bash
# Build fingerprint_concordance.csv from het_snps files
# Run after sample_fingerprint.lsf completes

FP_DIR="/sc/arion/scratch/kellej10/cadasil_public_ipsc/counts/fingerprint"
OUT="/sc/arion/work/kellej10/cadasil_public_ipsc/fingerprint_concordance.csv"

echo "sample_1,sample_2,shared,union,jaccard" > "$OUT"

for F1 in ${FP_DIR}/*_het_snps.txt; do
  S1=$(basename "$F1" | sed 's/_het_snps.txt//')
  for F2 in ${FP_DIR}/*_het_snps.txt; do
    S2=$(basename "$F2" | sed 's/_het_snps.txt//')
    if [[ "$S1" < "$S2" ]]; then
      SHARED=$(awk '{print $1"_"$2}' "$F1" "$F2" | sort | uniq -d | wc -l)
      UNION=$(awk '{print $1"_"$2}' "$F1" "$F2" | sort -u | wc -l)
      if [ "$UNION" -gt 0 ]; then
        JACCARD=$(echo "scale=4; $SHARED / $UNION" | bc)
      else
        JACCARD="0"
      fi
      echo "${S1},${S2},${SHARED},${UNION},${JACCARD}" >> "$OUT"
    fi
  done
done

echo "Written: $OUT"
cat "$OUT"
