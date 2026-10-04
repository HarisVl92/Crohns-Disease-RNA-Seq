#!/usr/bin/env bash
# =============================================================================
# 01b_qc_trim.sh - read QC and adapter/quality trimming.
#
#   FastQC (raw) -> Trim Galore -> FastQC (trimmed) -> MultiQC reports
#
# The Trim Galore parameters reproduce the original run, as recorded in
# data/trimmed/*_trimming_report.json: Trim Galore v2.2.0 with
# --quality 20 --stringency 1 --length 20 and the auto-detected Illumina adapter.
#
# Note: Trim Galore v2 auto-enabled poly-G trimming for 5 of 6 samples. These
# data come from a HiSeq 2000 (4-colour chemistry), so --no_poly_g is the cleaner
# choice for a new run; the effect on the present results is negligible.
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
use_env fastqc trim_galore multiqc

RAW="$PROJ/data/raw"
TRIM="$PROJ/data/trimmed"
QC_RAW="$PROJ/results/fastqc_raw"
QC_TRIM="$PROJ/results/fastqc_trimmed"
mkdir -p "$TRIM" "$QC_RAW" "$QC_TRIM" "$PROJ/logs"

RAW_FASTQ=()
TRIMMED=()
TO_TRIM=()
for S in "${SAMPLES[@]}"; do
  [ -s "$RAW/${S}.fastq.gz" ] || die "missing $RAW/${S}.fastq.gz (run 01_download_data.sh first)"
  RAW_FASTQ+=("$RAW/${S}.fastq.gz")
  TRIMMED+=("$TRIM/${S}_trimmed.fq.gz")
  [ -s "$TRIM/${S}_trimmed.fq.gz" ] || TO_TRIM+=("$RAW/${S}.fastq.gz")
done

log "FastQC on raw reads (${#RAW_FASTQ[@]} files)"
fastqc --quiet --threads "$THREADS" --outdir "$QC_RAW" "${RAW_FASTQ[@]}" > /dev/null
multiqc --force --quiet --cl-config "show_analysis_paths: false" \
  --filename multiqc_raw.html --outdir "$QC_RAW" "$QC_RAW"

if [ "${#TO_TRIM[@]}" -eq 0 ]; then
  log "Trimmed reads already present -> skip Trim Galore"
else
  log "Trim Galore (${#TO_TRIM[@]} files)"
  trim_galore --quality 20 --stringency 1 --length 20 --cores 4 \
    --output_dir "$TRIM" "${TO_TRIM[@]}" \
    > "$PROJ/logs/trim_galore_$(date +%Y%m%d_%H%M%S).log" 2>&1
fi

log "FastQC on trimmed reads"
fastqc --quiet --threads "$THREADS" --outdir "$QC_TRIM" "${TRIMMED[@]}" > /dev/null
multiqc --force --quiet --cl-config "show_analysis_paths: false" \
  --filename multiqc_trimmed.html --outdir "$QC_TRIM" "$QC_TRIM"
log "DONE: trimmed reads in $TRIM; QC reports in $QC_RAW and $QC_TRIM"
