#!/bin/bash
GTF=/sc/arion/scratch/kellej10/cadasil_public_ipsc/reference/gencode.v44.primary_assembly.annotation.gtf

echo "=== Transcript annotations for FLT1, FN1, NOTCH3 ==="
echo ""
grep -P '\ttranscript\t' "$GTF" | grep -E 'gene_name "(FLT1|FN1|NOTCH3)"' | perl -ne '/transcript_id "([^"]+)"/; $tid=$1; /transcript_name "([^"]+)"/; $tname=$1; /gene_name "([^"]+)"/; $gene=$1; /transcript_type "([^"]+)"/; $ttype=$1; @tags=($_ =~ /tag "([^"]+)"/g); printf "%s\t%s\t%s\t%s\t%s\n",$gene,$tid,$tname,$ttype,join(";",@tags);'

echo ""
echo "=== Exon counts and total exon length per transcript ==="
echo ""
awk -F'\t' '$3=="exon" && ($0 ~ /gene_name "FLT1"/ || $0 ~ /gene_name "FN1"/ || $0 ~ /gene_name "NOTCH3"/) {match($0,/transcript_id "([^"]+)"/,tid);match($0,/transcript_name "([^"]+)"/,tname);match($0,/gene_name "([^"]+)"/,gname);id=tid[1];names[id]=tname[1];genes[id]=gname[1];exon_count[id]++;exon_len[id]+=($5-$4+1)} END{for(id in exon_count) printf "%s\t%s\t%s\t%d\t%d\n",genes[id],id,names[id],exon_count[id],exon_len[id]}' "$GTF" | sort -k1,1 -k5,5rn
