#!/usr/bin/env bash
# =============================================================================
# 05_featurecounts.sh
# Gene-level counting με featureCounts (subread).
#
# Βήμα 1: ΑΝΙΧΝΕΥΣΗ STRANDEDNESS εμπειρικά -> τρέχουμε -s 0/1/2 σε ΕΝΑ δείγμα
#         και διαλέγουμε αυτό με το υψηλότερο assignment rate.
# Βήμα 2: Τελικό counting ΟΛΩΝ των δειγμάτων με το σωστό -s.
#
#  * single-end -> ΧΩΡΙΣ -p (το -p είναι fragment counting για paired-end).
#  * -t exon -g gene_id -> gene-level counts (default, σωστό για DESeq2).
# =============================================================================
set -euo pipefail

PROJ="/home/haris/Desktop/Bioinformatics_Project"
GTF="$PROJ/data/reference/gtf/GRCh38_r110/Homo_sapiens.GRCh38.110.chr.gtf"
BAM="$PROJ/results/alignment"
OUT="$PROJ/results/counts"
LOG="$PROJ/logs"
THREADS=16

mkdir -p "$OUT" "$LOG"
source ~/anaconda3/etc/profile.d/conda.sh
conda activate rnaseq

SAMPLES=(SRR1782721 SRR1782830 SRR1782732 SRR1782738 SRR1782739 SRR1782761)
BAMS=(); for S in "${SAMPLES[@]}"; do BAMS+=("$BAM/${S}.sorted.bam"); done

# --- Βήμα 1: ανίχνευση strandedness σε 1 δείγμα ---
TEST_BAM="$BAM/${SAMPLES[0]}.sorted.bam"
echo "[$(date +%T)] Ανίχνευση strandedness σε $(basename "$TEST_BAM")..."
declare -A ASSIGNED
for s in 0 1 2; do
  featureCounts -T "$THREADS" -s "$s" -t exon -g gene_id \
    -a "$GTF" -o "$OUT/_strandtest_s${s}.txt" "$TEST_BAM" \
    > "$LOG/05_strandtest_s${s}.log" 2>&1
  a=$(grep -P "^Assigned" "$OUT/_strandtest_s${s}.txt.summary" | awk '{print $2}')
  ASSIGNED[$s]=$a
  echo "   -s ${s}: assigned = ${a}"
done
# Διάλεξε το s με το max assigned
BEST_S=0; BEST_V=${ASSIGNED[0]}
for s in 1 2; do if (( ${ASSIGNED[$s]} > BEST_V )); then BEST_V=${ASSIGNED[$s]}; BEST_S=$s; fi; done
echo "[$(date +%T)] Επιλεγμένο strandedness: -s ${BEST_S}"
rm -f "$OUT"/_strandtest_s*.txt "$OUT"/_strandtest_s*.txt.summary

# --- Βήμα 2: τελικό counting όλων των δειγμάτων ---
echo "[$(date +%T)] featureCounts (όλα τα δείγματα, -s ${BEST_S})..."
featureCounts -T "$THREADS" -s "$BEST_S" -t exon -g gene_id \
  -a "$GTF" -o "$OUT/counts.txt" "${BAMS[@]}" \
  2> "$LOG/05_featurecounts.log"

echo "${BEST_S}" > "$OUT/strandedness.used.txt"
echo "[$(date +%T)] DONE."
echo "=== Assignment summary ==="
cat "$OUT/counts.txt.summary" | awk 'NR==1 || $2+$3+$4+$5+$6+$7>0'
