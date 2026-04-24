#!/bin/bash
# Check which exon position 15192182 falls in and look for R153C
GTF=/sc/arion/scratch/kellej10/cadasil_public_ipsc/reference/gencode.v44.primary_assembly.annotation.gtf

echo "=== NOTCH3 exons on chr19 (ENST00000263388, canonical) ==="
grep 'ENST00000263388' "$GTF" | awk -F'\t' '$3=="exon"' | awk -F'\t' '{print NR, $1, $4, $5, $5-$4+1, $7}' | sort -k3 -n

echo ""
echo "=== Which exon contains position 15192182? ==="
grep 'ENST00000263388' "$GTF" | awk -F'\t' '$3=="exon" && 15192182 >= $4 && 15192182 <= $5 {print "Position 15192182 is in:", $0}'

echo ""
echo "=== Which exon contains position 15191968? ==="
grep 'ENST00000263388' "$GTF" | awk -F'\t' '$3=="exon" && 15191968 >= $4 && 15191968 <= $5 {print "Position 15191968 is in:", $0}'

echo ""
echo "=== Check reference base at both positions ==="
module load samtools/1.21 2>/dev/null
REF=/sc/arion/scratch/kellej10/cadasil_public_ipsc/reference/GRCh38.primary_assembly.genome.fa
samtools faidx "$REF" chr19:15192182-15192182
samtools faidx "$REF" chr19:15191968-15191968

echo ""
echo "=== Scan XZ5 pileup for any position with VAF between 0.3-0.7 (heterozygous) ==="
echo "These are candidate R153C positions:"
awk '$4 >= 20 {
    n=split($5,a,"");
    ref=0;
    for(i=1;i<=n;i++){if(a[i]=="." || a[i]==",") ref++}
    alt=n-ref;
    vaf=alt/$4;
    if(vaf > 0.3 && vaf < 0.7) printf "%s\t%s\t%s\tdepth=%s\tref=%s\talt=%s\tVAF=%.4f\n",$1,$2,$3,$4,ref,alt,vaf
}' /sc/arion/scratch/kellej10/cadasil_public_ipsc/counts/XZ5_CADASIL_R153C_NOTCH3_pileup.txt
