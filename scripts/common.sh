#!/usr/bin/env bash
# =============================================================================
# common.sh - shared settings and helpers, sourced by every shell step.
#
#   PROJ      repository root (default: the parent folder of scripts/)
#   THREADS   CPU threads for multi-threaded tools (default: 16)
#   MANIFEST  sample sheet; the run accessions are read from its first column
#
# Any of them can be overridden from the environment, e.g.
#   THREADS=8 bash scripts/04_align.sh
# =============================================================================
PROJ="${PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# Work from the repository root with relative paths from here on, so that logs
# and outputs (count-matrix headers, MultiQC reports) contain no machine-specific paths.
cd "$PROJ"
PROJ="."
THREADS="${THREADS:-16}"
MANIFEST="${MANIFEST:-$PROJ/data/samples_manifest.tsv}"

# Reference genome: Ensembl GRCh38, release 110
ENSEMBL_RELEASE=110
REF_DIR="$PROJ/data/reference"
REF_FASTA="$REF_DIR/fasta/GRCh38_r${ENSEMBL_RELEASE}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
REF_GTF="$REF_DIR/gtf/GRCh38_r${ENSEMBL_RELEASE}/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.chr.gtf"
HISAT2_INDEX_DIR="$REF_DIR/hisat2_index/GRCh38_r${ENSEMBL_RELEASE}"

log() { echo "[$(date +%H:%M:%S)] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

# The sample list lives in exactly one place: the manifest (comment lines,
# the header and empty lines are skipped).
[ -f "$MANIFEST" ] || die "sample manifest not found: $MANIFEST"
mapfile -t SAMPLES < <(awk -F'\t' '!/^#/ && NF && $1 != "run_accession" {print $1}' "$MANIFEST")
[ "${#SAMPLES[@]}" -gt 0 ] || die "no samples listed in $MANIFEST"

# Activate the 'rnaseq' conda environment (when conda is available) and check
# that the required tools are on PATH.   Usage: use_env tool [tool ...]
use_env() {
  if [ "${CONDA_DEFAULT_ENV:-}" != "rnaseq" ] && command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    set +u  # conda activation scripts may read unset variables
    conda activate rnaseq
    set -u
  fi
  local tool
  for tool in "$@"; do
    command -v "$tool" >/dev/null 2>&1 \
      || die "'$tool' not found; create the environment with: conda env create -f environment.yml"
  done
}
