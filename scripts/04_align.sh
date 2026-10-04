#!/usr/bin/env bash
# =============================================================================
# 04_align.sh - HISAT2 spliced alignment (single-end) -> sorted, indexed BAM.
#
#   * --known-splicesite-infile supplies the known junctions (the memory-friendly
#     alternative to an index built with --ss/--exon, see 03_build_index.sh)
#   * single-end reads -> -U (not -1/-2)
#   * piped straight into samtools sort: no intermediate SAM file on disk
#   * read groups (--rg-id/--rg) so that the BAM files are properly annotated
#   * idempotent: samples whose BAM is already indexed are skipped
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
use_env hisat2 samtools

IDX="$HISAT2_INDEX_DIR/genome"
SPLICE_SITES="$HISAT2_INDEX_DIR/splicesites.txt"
TRIM="$PROJ/data/trimmed"
OUT="$PROJ/results/alignment"
LOG="$PROJ/logs/align"
SORT_THREADS="${SORT_THREADS:-6}"
mkdir -p "$OUT" "$LOG"

for S in "${SAMPLES[@]}"; do
  BAM="$OUT/${S}.sorted.bam"
  if [ -s "$BAM.bai" ]; then
    log "$S: $(basename "$BAM") already indexed -> skip"
    continue
  fi
  [ -s "$TRIM/${S}_trimmed.fq.gz" ] || die "missing trimmed reads for $S (run 01b_qc_trim.sh first)"
  log "Alignment: $S"
  hisat2 -p "$THREADS" -x "$IDX" \
         --known-splicesite-infile "$SPLICE_SITES" \
         -U "$TRIM/${S}_trimmed.fq.gz" \
         --rg-id "$S" --rg "SM:${S}" --rg "PL:ILLUMINA" \
         2> "$LOG/${S}.hisat2.log" \
    | samtools sort -@ "$SORT_THREADS" -o "$BAM" -
  samtools index -@ "$SORT_THREADS" "$BAM"
  log "  -> $(basename "$BAM") | $(grep 'overall alignment rate' "$LOG/${S}.hisat2.log" || true)"
done
log "ALIGNMENT DONE"
