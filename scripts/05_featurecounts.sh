#!/usr/bin/env bash
# =============================================================================
# 05_featurecounts.sh - gene-level read counting with featureCounts (Subread).
#
# Step 1: detect library strandedness empirically: count one sample with
#         -s 0/1/2 and keep the setting with the most assigned reads.
# Step 2: count all samples with that setting.
# Step 3: MultiQC report of the alignment and counting statistics.
#
#   * single-end data -> no -p (that flag counts fragments of paired-end reads)
#   * -t exon -g gene_id -> gene-level counts, as DESeq2 expects
#   * a gzipped copy of the count matrix is kept under version control
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
use_env featureCounts multiqc

BAM_DIR="$PROJ/results/alignment"
OUT="$PROJ/results/counts"
LOG="$PROJ/logs"
mkdir -p "$OUT" "$LOG"

BAMS=()
for S in "${SAMPLES[@]}"; do
  [ -s "$BAM_DIR/${S}.sorted.bam" ] || die "missing BAM for $S (run 04_align.sh first)"
  BAMS+=("$BAM_DIR/${S}.sorted.bam")
done

# --- Step 1: strandedness test on the first sample ---
log "Strandedness test on $(basename "${BAMS[0]}")"
ASSIGNED=()
for s in 0 1 2; do
  featureCounts -T "$THREADS" -s "$s" -t exon -g gene_id -a "$REF_GTF" \
    -o "$OUT/_strandtest_s${s}.txt" "${BAMS[0]}" > "$LOG/05_strandtest_s${s}.log" 2>&1
  ASSIGNED[s]=$(awk '$1 == "Assigned" {print $2}' "$OUT/_strandtest_s${s}.txt.summary")
  log "  -s $s: assigned = ${ASSIGNED[s]}"
done
BEST_S=0
for s in 1 2; do
  if (( ASSIGNED[s] > ASSIGNED[BEST_S] )); then BEST_S=$s; fi
done
log "Selected strandedness: -s $BEST_S"
rm -f "$OUT"/_strandtest_s*.txt "$OUT"/_strandtest_s*.txt.summary

# --- Step 2: count all samples ---
log "featureCounts, all samples (-s $BEST_S)"
featureCounts -T "$THREADS" -s "$BEST_S" -t exon -g gene_id -a "$REF_GTF" \
  -o "$OUT/counts.txt" "${BAMS[@]}" 2> "$LOG/05_featurecounts.log"
echo "$BEST_S" > "$OUT/strandedness.used.txt"
gzip -9 -n -k -f "$OUT/counts.txt"

# --- Step 3: MultiQC report (HISAT2 logs + featureCounts summary) ---
multiqc --force --quiet --cl-config "show_analysis_paths: false" \
  --filename alignment_multiqc.html --outdir "$PROJ/results/multiqc" \
  "$PROJ/logs/align" "$OUT/counts.txt.summary"

log "DONE. Assignment summary:"
awk 'NR == 1 || $2 + $3 + $4 + $5 + $6 + $7 > 0' "$OUT/counts.txt.summary"
