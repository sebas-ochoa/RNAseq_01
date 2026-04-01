#!/usr/bin/env bash
#SBATCH --job-name=15_subsample
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
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
module load "$SEQTK_MODULE"

mkdir -p "$SUB15_DIR" "04_logs/slurm"
chmod -R u+rwX,g+rwX "$SUB15_DIR" 2>/dev/null || true

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES_TSV")
[[ -n "$LINE" ]] || exit 0
SAMPLE=$(echo "$LINE" | cut -f1)

R1="${TRIM14_DIR}/${SAMPLE}_1.trim.fq.gz"
R2="${TRIM14_DIR}/${SAMPLE}_2.trim.fq.gz"
OUT_R1="${SUB15_DIR}/${SAMPLE}_1.sub.fq.gz"
OUT_R2="${SUB15_DIR}/${SAMPLE}_2.sub.fq.gz"

echo "[INFO] Sample: $SAMPLE — subsampling to ${SUBSAMPLE_READS} reads"

seqtk sample -s 100 "$R1" "$SUBSAMPLE_READS" | gzip -c > "$OUT_R1"
seqtk sample -s 100 "$R2" "$SUBSAMPLE_READS" | gzip -c > "$OUT_R2"

echo "[INFO] Done: $SAMPLE"
