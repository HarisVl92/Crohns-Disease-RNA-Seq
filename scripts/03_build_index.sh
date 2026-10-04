#!/usr/bin/env bash
# =============================================================================
# 03_build_index.sh
# Εξαγωγή splice sites από το GTF + χτίσιμο HISAT2 index (plain HFM).
#
# ΑΠΟΦΑΣΗ ΜΝΗΜΗΣ (64GB RAM):
#   Το πλήρες HGFM index (hisat2-build --ss --exon) απαιτεί ~160GB RAM για
#   ανθρώπινο γονιδίωμα -> ΑΔΥΝΑΤΟ εδώ. Αντ' αυτού χτίζουμε το ΑΠΛΟ index
#   (~6GB RAM) και δίνουμε τα known splice sites στο ALIGNMENT βήμα μέσω
#   --known-splicesite-infile. Ίδιο αποτέλεσμα spliced alignment, ρεαλιστική μνήμη.
# =============================================================================
set -euo pipefail

PROJ="/home/haris/Desktop/Bioinformatics_Project"
REF_FASTA="$PROJ/data/reference/fasta/GRCh38_r110/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
REF_GTF="$PROJ/data/reference/gtf/GRCh38_r110/Homo_sapiens.GRCh38.110.chr.gtf"
IDX_DIR="$PROJ/data/reference/hisat2_index/GRCh38_r110"
LOG="$PROJ/logs"
THREADS=16

mkdir -p "$IDX_DIR" "$LOG"
source ~/anaconda3/etc/profile.d/conda.sh
conda activate rnaseq

# --- 1) Splice sites από το GTF (χρησιμοποιούνται στο alignment, όχι στο index) ---
echo "[$(date +%T)] Εξαγωγή splice sites από το GTF..."
hisat2_extract_splice_sites.py "$REF_GTF" > "$IDX_DIR/splicesites.txt"
echo "[$(date +%T)] Splice sites: $(wc -l < "$IDX_DIR/splicesites.txt") junctions"

# --- 2) Plain HISAT2 index ---
echo "[$(date +%T)] Χτίσιμο HISAT2 index (plain HFM, -p $THREADS)..."
hisat2-build -p "$THREADS" "$REF_FASTA" "$IDX_DIR/genome" 2> "$LOG/03_hisat2_build.log"

echo "[$(date +%T)] ΟΛΟΚΛΗΡΩΘΗΚΕ. Αρχεία index:"
ls -lh "$IDX_DIR/"
