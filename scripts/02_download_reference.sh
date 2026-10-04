#!/usr/bin/env bash
# =============================================================================
# 02_download_reference.sh - GRCh38 genome and annotation (Ensembl release 110).
#
#   FASTA  dna_sm.primary_assembly: soft-masked, primary assembly only
#          (no alt loci, which would otherwise inflate multi-mapping)
#   GTF    .chr: chromosomes 1-22, X, Y and MT only, matching the primary assembly
#
# Ensembl CHECKSUMS files use BSD `sum`, not MD5. The script also writes the
# uncompressed copies read by hisat2-build and featureCounts, and the gene
# annotation table (gene_id -> gene_name, biotype) used by the DESeq2 report.
# Idempotent: existing files are verified instead of downloaded again.
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

FA_BASE="http://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna"
GTF_BASE="http://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/homo_sapiens"
ANNOTATION="$REF_DIR/gene_annotation.tsv"
LOG="$PROJ/logs/reference_download_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$(dirname "$REF_FASTA")" "$(dirname "$REF_GTF")" "$PROJ/logs"
exec > >(tee -a "$LOG") 2>&1

verify() {  # $1 = local file, $2 = CHECKSUMS URL (fetched fresh every time)
  local f="$1" url="$2" base want got
  base=$(basename "$f")
  want=$(curl -s --max-time 60 "$url" | grep " ${base}\$" | awk '{print $1, $2}')
  [ -n "$want" ] || die "no checksum for $base in $url"
  got=$(sum "$f" | awk '{print $1, $2}')
  [ "$got" = "$want" ] || die "checksum mismatch for $base (expected '$want', got '$got')"
  log "$base: BSD sum OK ($got)"
}

fetch() {  # $1 = local .gz path, $2 = base URL; download unless present, always verify
  local f="$1" url="$2" name
  name=$(basename "$f")
  if [ -f "$f" ]; then
    log "$name: already present -> verify"
  else
    log "download $name"
    wget -q -c -O "$f" "$url/$name" || die "wget failed: $name"
  fi
  verify "$f" "$url/CHECKSUMS"
}

log "Reference GRCh38, Ensembl release ${ENSEMBL_RELEASE} -> $REF_DIR"
fetch "${REF_FASTA}.gz" "$FA_BASE"
fetch "${REF_GTF}.gz" "$GTF_BASE"

# hisat2-build and featureCounts read the uncompressed files.
for f in "$REF_FASTA" "$REF_GTF"; do
  if [ -s "$f" ]; then
    log "$(basename "$f"): uncompressed copy present"
  else
    log "decompress $(basename "$f").gz"
    gunzip -k -f "$f.gz"
  fi
done

# Gene annotation table for the DESeq2 report; genes without a name keep their Ensembl ID.
if [ -s "$ANNOTATION" ]; then
  log "$(basename "$ANNOTATION"): present -> skip"
else
  log "build $(basename "$ANNOTATION") from the GTF"
  awk -F'\t' '
    function attr(s, key) {
      if (match(s, key " \"[^\"]*\""))
        return substr(s, RSTART + length(key) + 2, RLENGTH - length(key) - 3)
      return ""
    }
    BEGIN { OFS = "\t"; print "gene_id", "gene_name", "biotype" }
    $3 == "gene" {
      id = attr($9, "gene_id"); name = attr($9, "gene_name")
      print id, (name == "" ? id : name), attr($9, "gene_biotype")
    }' "$REF_GTF" > "$ANNOTATION"
fi
log "DONE: FASTA, GTF and gene annotation (Ensembl ${ENSEMBL_RELEASE})"
