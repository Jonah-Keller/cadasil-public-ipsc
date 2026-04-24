#!/bin/bash
# Scan NOTCH3 pileup for positions with significant non-reference alleles
COUNTS_DIR="/sc/arion/scratch/kellej10/cadasil_public_ipsc/counts"

for SAMPLE in XZ5_CADASIL_R153C WZ16_CADASIL_C224Y WZ12_CADASIL_Control; do
    echo "=== $SAMPLE — top variant positions in NOTCH3 ==="
    awk '$4 >= 10 {
        n=split($5,a,"");
        ref=0;
        for(i=1;i<=n;i++){if(a[i]=="." || a[i]==",") ref++}
        alt=n-ref;
        if(alt>5) printf "%s\t%s\t%s\t%s\t%s\t%s\t%.4f\n",$1,$2,$3,$4,ref,alt,alt/$4
    }' "${COUNTS_DIR}/${SAMPLE}_NOTCH3_pileup.txt" | sort -t$'\t' -k7 -rn | head -10
    echo ""
done

echo "=== Direct check at expected variant positions ==="
echo "Position 15192182 (expected R153C c.457C>T, genomic G>A):"
for f in ${COUNTS_DIR}/*_NOTCH3_pileup.txt; do
    SAMPLE=$(basename "$f" | sed 's/_NOTCH3_pileup.txt//')
    LINE=$(awk '$2==15192182' "$f")
    if [ -n "$LINE" ]; then
        echo "  $SAMPLE: $LINE"
    else
        echo "  $SAMPLE: no coverage"
    fi
done

echo ""
echo "Position 15191968 (expected C224Y c.671G>A, genomic C>T):"
for f in ${COUNTS_DIR}/*_NOTCH3_pileup.txt; do
    SAMPLE=$(basename "$f" | sed 's/_NOTCH3_pileup.txt//')
    LINE=$(awk '$2==15191968' "$f")
    if [ -n "$LINE" ]; then
        echo "  $SAMPLE: $LINE"
    else
        echo "  $SAMPLE: no coverage"
    fi
done
