#!/usr/bin/env bash
#SBATCH --job-name=qualimap_bamqc
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3G
#SBATCH --time=01:00:00
#SBATCH --output=04_logs/slurm/%x_%A_%a.out
#SBATCH --error=04_logs/slurm/%x_%A_%a.err
#SBATCH --array=1-12

set -euo pipefail
umask 007

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"
source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

QUALIMAP_DIR="${QUALIMAP_DIR:-${PROJECT_ROOT}/02_work/05_align/qualimap}"
QUALIMAP_MEM="${QUALIMAP_MEM:-10G}"

module purge
module load "$STACK_MODULE"
module load "$QUALIMAP_MODULE"

mkdir -p "$SLURM_LOG_DIR" "$QUALIMAP_DIR"
chmod -R u+rwX,g+rwX "$SLURM_LOG_DIR" "$QUALIMAP_DIR" 2>/dev/null || true

LINE="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES_TSV")"
[[ -n "$LINE" ]] || exit 0

SAMPLE="$(echo "$LINE" | cut -f1)"
BAM="${ALIGN_DIR}/${SAMPLE}.sorted.bam"

[[ -s "$BAM" ]] || { echo "[ERROR] BAM not found: $BAM" >&2; exit 1; }

OUT_DIR="${QUALIMAP_DIR}/${SAMPLE}"
[[ -s "${OUT_DIR}/qualimapReport.html" ]] && exit 0

mkdir -p "$OUT_DIR"

qualimap bamqc \
  -bam "$BAM" \
  -nt "$SLURM_CPUS_PER_TASK" \
  -outdir "$OUT_DIR" \
  -outformat HTML \
  --java-mem-size="$QUALIMAP_MEM" \
  -c \
  --outside-stats

echo "[INFO] Done: ${OUT_DIR}/qualimapReport.html"
