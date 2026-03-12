#!/usr/bin/env bash
#SBATCH --job-name=11_fastqc
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --time=02:00:00
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

OUT="${PROJECT_ROOT}/02_work/11_qc_raw/fastqc"
mkdir -p "$OUT" "04_logs/slurm"
chmod -R u+rwX,g+rwX "$OUT" 2>/dev/null || true

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${PROJECT_ROOT}/00_admin/metadata/samples.tsv")
[[ -n "$LINE" ]] || exit 0
SAMPLE=$(echo "$LINE" | cut -f1)
R1=$(echo "$LINE" | cut -f2)
R2=$(echo "$LINE" | cut -f3)

fastqc -t "$SLURM_CPUS_PER_TASK" -o "$OUT" "$R1" "$R2"

echo "[INFO] Done: $SAMPLE"
