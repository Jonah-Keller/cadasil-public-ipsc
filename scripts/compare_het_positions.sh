#!/bin/bash
# Compare heterozygous positions between CADASIL samples and controls
# to distinguish real mutations from common SNPs
COUNTS_DIR="/sc/arion/scratch/kellej10/cadasil_public_ipsc/counts"

echo "=== Heterozygous positions (VAF 0.3-0.7, depth>=20) unique to XZ5_R153C ==="
echo "(Positions NOT heterozygous in controls or WZ16_C224Y)"
echo ""

# Get het positions for XZ5
awk '$4 >= 20 {
    n=split($5,a,""); ref=0;
    for(i=1;i<=n;i++){if(a[i]=="." || a[i]==",") ref++}
    alt=n-ref; vaf=alt/$4;
    if(vaf > 0.3 && vaf < 0.7) print $2
}' ${COUNTS_DIR}/XZ5_CADASIL_R153C_NOTCH3_pileup.txt > /tmp/xz5_het.txt

# Get het positions for WZ16
awk '$4 >= 20 {
    n=split($5,a,""); ref=0;
    for(i=1;i<=n;i++){if(a[i]=="." || a[i]==",") ref++}
    alt=n-ref; vaf=alt/$4;
    if(vaf > 0.3 && vaf < 0.7) print $2
}' ${COUNTS_DIR}/WZ16_CADASIL_C224Y_NOTCH3_pileup.txt > /tmp/wz16_het.txt

# Get het positions for control WZ12
awk '$4 >= 20 {
    n=split($5,a,""); ref=0;
    for(i=1;i<=n;i++){if(a[i]=="." || a[i]==",") ref++}
    alt=n-ref; vaf=alt/$4;
    if(vaf > 0.3 && vaf < 0.7) print $2
}' ${COUNTS_DIR}/WZ12_CADASIL_Control_NOTCH3_pileup.txt > /tmp/wz12_het.txt

echo "Het positions: XZ5=$(wc -l < /tmp/xz5_het.txt) WZ16=$(wc -l < /tmp/wz16_het.txt) WZ12_ctrl=$(wc -l < /tmp/wz12_het.txt)"

echo ""
echo "=== Positions heterozygous in BOTH CADASIL samples but NOT controls ==="
comm -12 <(sort /tmp/xz5_het.txt) <(sort /tmp/wz16_het.txt) | while read pos; do
    if ! grep -q "^${pos}$" /tmp/wz12_het.txt; then
        echo -n "pos=$pos  "
        echo -n "XZ5: "
        awk -v p="$pos" '$2==p {printf "depth=%s ref=%s ", $4, $3}' ${COUNTS_DIR}/XZ5_CADASIL_R153C_NOTCH3_pileup.txt
        echo -n "WZ16: "
        awk -v p="$pos" '$2==p {printf "depth=%s ref=%s ", $4, $3}' ${COUNTS_DIR}/WZ16_CADASIL_C224Y_NOTCH3_pileup.txt
        echo ""
    fi
done

echo ""
echo "=== Positions heterozygous ONLY in XZ5 (not in WZ16 or controls) ==="
comm -23 <(sort /tmp/xz5_het.txt) <(sort -u /tmp/wz16_het.txt /tmp/wz12_het.txt) | while read pos; do
    awk -v p="$pos" '$2==p {
        n=split($5,a,""); ref=0;
        for(i=1;i<=n;i++){if(a[i]=="." || a[i]==",") ref++}
        alt=n-ref; vaf=alt/$4;
        printf "pos=%s ref_base=%s depth=%s ref_reads=%s alt_reads=%s VAF=%.4f\n", $2, $3, $4, ref, alt, vaf
    }' ${COUNTS_DIR}/XZ5_CADASIL_R153C_NOTCH3_pileup.txt
done
