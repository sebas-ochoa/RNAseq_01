#!/usr/bin/env bash
#SBATCH --job-name=16_fastqc
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00
#SBATCH --output=04_logs/slurm/%x_%A_%a.out
#SBATCH --error=04_logs/slurm/%x_%A_%a.err
#SBATCH --array=1-12

set -euo pipefail
umask 007

module purge
module load stack/2024-06
module load fastqc

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"
source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

OUT="${QC16_DIR}/fastqc"
mkdir -p "$OUT" "04_logs/slurm"
chmod -R u+rwX,g+rwX "$OUT" 2>/dev/null || true

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES_TSV")
[[ -n "$LINE" ]] || exit 0
SAMPLE=$(echo "$LINE" | cut -f1)

R1="${SUB15_DIR}/${SAMPLE}_1.sub.fq.gz"
R2="${SUB15_DIR}/${SAMPLE}_2.sub.fq.gz"

fastqc -t "$SLURM_CPUS_PER_TASK" -o "$OUT" "$R1" "$R2"

echo "[INFO] Done: $SAMPLE"
