#!/usr/bin/env bash
# =============================================================================
# 04_align.sh
# HISAT2 spliced alignment (single-end) -> sorted+indexed BAM ανά δείγμα.
#
#  * --known-splicesite-infile: δίνει τα known junctions στο alignment (η
#    memory-friendly εναλλακτική του --ss/--exon στο index).
#  * single-end -> -U (όχι -1/-2).
#  * pipe απευθείας στο samtools sort (χωρίς ενδιάμεσο SAM στον δίσκο).
#  * Read groups (--rg) για σωστά annotated BAMs.
# =============================================================================
set -euo pipefail

PROJ="/home/haris/Desktop/Bioinformatics_Project"
IDX="$PROJ/data/reference/hisat2_index/GRCh38_r110/genome"
SS="$PROJ/data/reference/hisat2_index/GRCh38_r110/splicesites.txt"
TRIM="$PROJ/data/trimmed"
OUT="$PROJ/results/alignment"
LOG="$PROJ/logs/align"
THREADS=16
SORT_THREADS=6

mkdir -p "$OUT" "$LOG"
source ~/anaconda3/etc/profile.d/conda.sh
conda activate rnaseq

SAMPLES=(SRR1782721 SRR1782830 SRR1782732 SRR1782738 SRR1782739 SRR1782761)

for S in "${SAMPLES[@]}"; do
  echo "[$(date +%T)] Alignment: $S"
  hisat2 -p "$THREADS" \
         -x "$IDX" \
         --known-splicesite-infile "$SS" \
         -U "$TRIM/${S}_trimmed.fq.gz" \
         --rg-id "$S" --rg "SM:${S}" --rg "PL:ILLUMINA" \
         2> "$LOG/${S}.hisat2.log" \
    | samtools sort -@ "$SORT_THREADS" -o "$OUT/${S}.sorted.bam" -
  samtools index -@ "$SORT_THREADS" "$OUT/${S}.sorted.bam"
  rate=$(grep "overall alignment rate" "$LOG/${S}.hisat2.log" || true)
  echo "[$(date +%T)]   -> ${S}.sorted.bam | $rate"
done

echo "[$(date +%T)] ALIGNMENT DONE"
