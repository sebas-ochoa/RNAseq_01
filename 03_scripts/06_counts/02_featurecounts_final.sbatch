#!/usr/bin/env bash
#SBATCH --job-name=fc_final
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --output=04_logs/slurm/%x_%j.out
#SBATCH --error=04_logs/slurm/%x_%j.err

set -euo pipefail
umask 007

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"
source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

module purge
module load "$STACK_MODULE"
module load "$GCC_MODULE"
module load "$SUBREAD_MODULE"

mkdir -p "$SLURM_LOG_DIR" "$COUNTS_DIR"
chmod -R u+rwX,g+rwX "$SLURM_LOG_DIR" "$COUNTS_DIR" 2>/dev/null || true

mapfile -t BAMS < <(find "$ALIGN_DIR" -maxdepth 1 -type f -name '*.sorted.bam' | sort)

OUT_COUNTS="${COUNTS_DIR}/gene_counts.tsv"
STRAND_MODE=0

echo "[INFO] Running final featureCounts with -s ${STRAND_MODE}"
featureCounts \
  -T "$SLURM_CPUS_PER_TASK" \
  -p --countReadPairs -B -C \
  -t exon -g gene_id \
  -s "$STRAND_MODE" \
  -a "$GTF" \
  -o "$OUT_COUNTS" \
  "${BAMS[@]}"

MATRIX_OUT="${COUNTS_DIR}/gene_counts.matrix.tsv"
awk 'BEGIN{OFS="\t"}
  NR==2 {
    printf "Geneid"
    for (i=7; i<=NF; i++) {
      sample=$i
      sub(/^.*\//, "", sample)
      sub(/\.sorted\.bam$/, "", sample)
      printf OFS sample
    }
    printf "\n"
    next
  }
  NR>2 {
    printf "%s", $1
    for (i=7; i<=NF; i++) {
      printf OFS $i
    }
    printf "\n"
  }
' "$OUT_COUNTS" > "$MATRIX_OUT"

echo "[INFO] Final counts done: $OUT_COUNTS"
echo "[INFO] Matrix done: $MATRIX_OUT"
