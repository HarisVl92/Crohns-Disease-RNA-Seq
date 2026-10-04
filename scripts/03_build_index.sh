#!/usr/bin/env bash
# =============================================================================
# 03_build_index.sh - HISAT2 index for GRCh38 plus a known splice-site file.
#
# Memory decision (64 GB RAM): a graph index that includes splice sites and
# exons (hisat2-build --ss --exon) needs ~160 GB RAM for the human genome.
# Instead, a plain index is built (~6 GB RAM, ~8 min with 16 threads) and the
# known splice sites are supplied at alignment time with
# --known-splicesite-infile: the same spliced alignment at a realistic cost.
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
use_env hisat2-build hisat2_extract_splice_sites.py

[ -s "$REF_FASTA" ] && [ -s "$REF_GTF" ] || die "reference missing (run 02_download_reference.sh first)"
mkdir -p "$HISAT2_INDEX_DIR" "$PROJ/logs"

log "Extract known splice sites from the GTF"
hisat2_extract_splice_sites.py "$REF_GTF" > "$HISAT2_INDEX_DIR/splicesites.txt"
log "Splice sites: $(wc -l < "$HISAT2_INDEX_DIR/splicesites.txt") junctions"

log "Build the plain HISAT2 index (-p $THREADS)"
hisat2-build -p "$THREADS" "$REF_FASTA" "$HISAT2_INDEX_DIR/genome" 2> "$PROJ/logs/03_hisat2_build.log"
log "DONE: index files in $HISAT2_INDEX_DIR"
