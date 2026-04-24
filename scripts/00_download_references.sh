#!/bin/bash
# Download Gencode GRCh38 reference genome and annotation
# Run on login node (needs internet) or interactive job

set -euo pipefail

REF_DIR="/sc/arion/scratch/kellej10/cadasil_public_ipsc/reference"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

echo "=== Downloading Gencode v44 GRCh38 primary assembly ==="
# Genome FASTA
wget -nc https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
# GTF annotation
wget -nc https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz
# Transcriptome FASTA (for Salmon, if needed later)
wget -nc https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz

echo "=== Decompressing ==="
gunzip -k GRCh38.primary_assembly.genome.fa.gz 2>/dev/null || true
gunzip -k gencode.v44.primary_assembly.annotation.gtf.gz 2>/dev/null || true
gunzip -k gencode.v44.transcripts.fa.gz 2>/dev/null || true

echo "=== Done. Files in $REF_DIR ==="
ls -lh "$REF_DIR"
