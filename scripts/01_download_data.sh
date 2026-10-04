#!/usr/bin/env bash
# =============================================================================
# 01_download_data.sh - download the raw single-end reads from ENA.
#
# Dataset: GEO GSE57945 / SRA SRP042228 (Haberman et al. 2014, RISK cohort),
# treatment-naive paediatric ileal biopsies: 3 Crohn's disease, 3 non-IBD controls.
#
#   * FASTQ URLs and MD5 checksums are requested fresh from the ENA API, so
#     the script keeps working if ENA reorganises its FTP paths.
#   * Idempotent: files already present with the right MD5 are skipped.
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

RAW="$PROJ/data/raw"
ENA="https://www.ebi.ac.uk/ena/portal/api/filereport"
LOG="$PROJ/logs/download_data_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$RAW" "$PROJ/logs"
exec > >(tee -a "$LOG") 2>&1

# ENA adds run_accession as the first column, so fields are located by header name.
column_of() {  # $1 = TSV header line, $2 = field name -> 1-based column number
  printf '%s\n' "$1" | tr '\t' '\n' | grep -nx "$2" | cut -d: -f1
}

log "Destination: $RAW (${#SAMPLES[@]} samples)"
for run in "${SAMPLES[@]}"; do
  resp=$(curl -s --max-time 60 -G "$ENA" \
      --data-urlencode "accession=$run" \
      --data-urlencode "result=read_run" \
      --data-urlencode "fields=fastq_ftp,fastq_md5" \
      --data-urlencode "format=tsv")
  hdr=$(printf '%s\n' "$resp" | head -1)
  row=$(printf '%s\n' "$resp" | sed -n '2p')
  iftp=$(column_of "$hdr" fastq_ftp) && imd5=$(column_of "$hdr" fastq_md5) \
    || die "$run: unexpected ENA response: $resp"
  url=$(printf '%s\n' "$row" | cut -f"$iftp")
  md5=$(printf '%s\n' "$row" | cut -f"$imd5")
  [ -n "$url" ] || die "$run: empty fastq_ftp field from ENA"
  out="$RAW/${run}.fastq.gz"

  if [ -f "$out" ] && [ "$(md5sum "$out" | cut -d' ' -f1)" = "$md5" ]; then
    log "$run: already present, MD5 OK -> skip"
    continue
  fi
  log "$run: download -> $out"
  wget -q -c -O "$out" "https://$url" || die "$run: wget failed"
  got=$(md5sum "$out" | cut -d' ' -f1)
  [ "$got" = "$md5" ] || die "$run: MD5 mismatch (expected $md5, got $got)"
  log "$run: MD5 OK, size $(du -h "$out" | cut -f1)"
done
log "DONE: $(ls -1 "$RAW"/*.fastq.gz 2>/dev/null | wc -l) FASTQ files in $RAW"
