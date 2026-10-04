#!/usr/bin/env bash
# 02_download_reference.sh — GRCh38 reference genome (Ensembl release 110)
#   FASTA: dna_sm.primary_assembly  (soft-masked + primary -> αποφυγή alt-loci multimapping)
#   GTF:   .chr                      (μόνο 1-22,X,Y,MT -> ταιριάζει με το primary_assembly)
# Επαλήθευση: τα Ensembl CHECKSUMS είναι BSD `sum` (ΟΧΙ md5).
# Αναμενόμενες τιμές (από Session 10):  FASTA -> 19039 920395  |  GTF -> 04474 53052
set -euo pipefail

PROJ="/home/haris/Desktop/Bioinformatics_Project"
REL=110
FA_DIR="$PROJ/data/reference/fasta/GRCh38_r${REL}"
GTF_DIR="$PROJ/data/reference/gtf/GRCh38_r${REL}"
FA_BASE="http://ftp.ensembl.org/pub/release-${REL}/fasta/homo_sapiens/dna"
GTF_BASE="http://ftp.ensembl.org/pub/release-${REL}/gtf/homo_sapiens"
FA="Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
GTF="Homo_sapiens.GRCh38.${REL}.chr.gtf.gz"
LOG="$PROJ/logs/reference_download_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$FA_DIR" "$GTF_DIR" "$PROJ/logs"
exec > >(tee -a "$LOG") 2>&1
log(){ echo "[$(date +%H:%M:%S)] $*"; }
die(){ echo "ERROR: $*" >&2; exit 1; }

verify(){ # $1=local file  $2=CHECKSUMS url  — φρέσκο CHECKSUMS κάθε φορά (όχι stale local)
  local f="$1" url="$2" base want got
  base=$(basename "$f")
  want=$(curl -s --max-time 60 "$url" | grep " ${base}\$" | awk '{print $1, $2}')
  [ -n "$want" ] || die "δεν βρέθηκε checksum για $base στο $url"
  got=$(sum "$f" | awk '{print $1, $2}')
  [ "$got" = "$want" ] || die "checksum mismatch $base (want='$want' got='$got')"
  log "$base: BSD sum OK -> $got"
}
fetch(){ # $1=dir  $2=base_url  $3=filename  — idempotent (skip download αν υπάρχει, αλλά πάντα verify)
  local dir="$1" burl="$2" name="$3" f="$1/$3"
  if [ -f "$f" ]; then log "$name: υπάρχει ήδη -> verify"
  else log "download $name"; wget -q -c -O "$f" "$burl/$name" || die "wget απέτυχε: $name"; fi
  verify "$f" "$burl/CHECKSUMS"
}

log "Reference GRCh38 Ensembl r${REL} -> $PROJ/data/reference"
fetch "$FA_DIR"  "$FA_BASE"  "$FA"
fetch "$GTF_DIR" "$GTF_BASE" "$GTF"
log "ΟΛΟΚΛΗΡΩΘΗΚΕ: FASTA + GTF (release ${REL})"
ls -lh "$FA_DIR/$FA" "$GTF_DIR/$GTF"
