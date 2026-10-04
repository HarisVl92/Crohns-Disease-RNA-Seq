#!/usr/bin/env bash
# 01_download_data.sh
# Κατέβασμα host ileal RNA-Seq reads — GEO GSE57945 / SRA SRP042228 (Haberman 2014, RISK cohort)
# 3x Crohn (CD, macroscopic inflammation) + 3x Not-IBD control · single-end · treatment-naive
# Πρότυπο: ίδια φιλοσοφία με το reference-download script (strict mode, logging, idempotency, checksum).
set -euo pipefail

PROJ="/home/haris/Desktop/Bioinformatics_Project"
RAW="$PROJ/data/raw"
ENA="https://www.ebi.ac.uk/ena/portal/api/filereport"
LOG="$PROJ/logs/download_data_$(date +%Y%m%d_%H%M%S).log"

RUNS=(SRR1782721 SRR1782830 SRR1782732 SRR1782738 SRR1782739 SRR1782761)

mkdir -p "$RAW" "$PROJ/logs"
exec > >(tee -a "$LOG") 2>&1
log(){ echo "[$(date +%H:%M:%S)] $*"; }
die(){ echo "ERROR: $*" >&2; exit 1; }

log "Προορισμός: $RAW  (συνολικά ${#RUNS[@]} samples)"
for run in "${RUNS[@]}"; do
  # Φρέσκο, authoritative URL + md5 από το ENA (όχι hardcoded — επιβίωση σε αλλαγές path)
  # ΠΡΟΣΟΧΗ: το ENA προσθέτει αυτόματα στήλη run_accession ΠΡΩΤΗ -> μην βασίζεσαι σε θέση.
  # Parse με βάση το header name (αλάνθαστο ανεξαρτήτως σειράς στηλών).
  resp=$(curl -s --max-time 60 -G "$ENA" \
      --data-urlencode "accession=$run" \
      --data-urlencode "result=read_run" \
      --data-urlencode "fields=fastq_ftp,fastq_md5" \
      --data-urlencode "format=tsv")
  hdr=$(printf '%s\n' "$resp" | head -1)
  row=$(printf '%s\n' "$resp" | sed -n '2p')
  iftp=$(printf '%s\n' "$hdr" | tr '\t' '\n' | grep -nx 'fastq_ftp' | cut -d: -f1)
  imd5=$(printf '%s\n' "$hdr" | tr '\t' '\n' | grep -nx 'fastq_md5' | cut -d: -f1)
  url=$(printf '%s\n' "$row" | cut -f"$iftp")
  md5=$(printf '%s\n' "$row" | cut -f"$imd5")
  [ -n "$url" ] || die "$run: κενό fastq_ftp από ENA"
  out="$RAW/${run}.fastq.gz"

  if [ -f "$out" ] && [ "$(md5sum "$out" | cut -d' ' -f1)" = "$md5" ]; then
    log "$run: υπάρχει ήδη + md5 OK → skip"; continue
  fi
  log "$run: download → $out"
  wget -q -c -O "$out" "https://$url" || die "$run: wget απέτυχε"
  got=$(md5sum "$out" | cut -d' ' -f1)
  [ "$got" = "$md5" ] || die "$run: md5 mismatch (exp=$md5 got=$got)"
  log "$run: md5 OK · μέγεθος $(du -h "$out" | cut -f1)"
done
log "ΟΛΟΚΛΗΡΩΘΗΚΕ — $(ls -1 "$RAW"/*.fastq.gz 2>/dev/null | wc -l) αρχεία στο $RAW"
