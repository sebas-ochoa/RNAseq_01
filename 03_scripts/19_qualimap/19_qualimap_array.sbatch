#!/usr/bin/env bash
#SBATCH --job-name=19_qualimap
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

module purge
module load "$STACK_MODULE"
export PATH="${JAVA_HOME}/bin:${PATH}"

mkdir -p "$QUALIMAP19_DIR" "04_logs/slurm"
chmod -R u+rwX,g+rwX "$QUALIMAP19_DIR" 2>/dev/null || true

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES_TSV")
[[ -n "$LINE" ]] || exit 0
SAMPLE=$(echo "$LINE" | cut -f1)

BAM="${ALIGN17_DIR}/${SAMPLE}.sorted.bam"
OUT_DIR="${QUALIMAP19_DIR}/${SAMPLE}"
mkdir -p "$OUT_DIR"

echo "[INFO] Sample: $SAMPLE"

"${QUALIMAP_HOME}/qualimap" bamqc \
  -bam "$BAM" \
  -nt "$SLURM_CPUS_PER_TASK" \
  -outdir "$OUT_DIR" \
  -outformat HTML \
  --java-mem-size="$QUALIMAP_MEM" \
  -c \
  --outside-stats

echo "[INFO] Done: ${OUT_DIR}/qualimapReport.html"
