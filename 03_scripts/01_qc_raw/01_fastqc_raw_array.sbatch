#!/usr/bin/env bash
#SBATCH --job-name=fastqc_raw
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --time=02:00:00
#SBATCH --output=04_logs/slurm/%x_%A_%a.out
#SBATCH --error=04_logs/slurm/%x_%A_%a.err
#SBATCH --array=1-12

set -euo pipefail

module purge
module load stack/2024-06
module load fastqc

# Absolute project root (where you ran sbatch)
PROJECT_ROOT="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is empty}"
cd "$PROJECT_ROOT"

OUT_FASTQC="${PROJECT_ROOT}/02_work/01_qc_raw/fastqc"
mkdir -p "$OUT_FASTQC"

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${PROJECT_ROOT}/00_admin/metadata/samples.tsv")
SAMPLE=$(echo "$LINE" | cut -f1)
R1=$(echo "$LINE" | cut -f2)
R2=$(echo "$LINE" | cut -f3)

fastqc -t "$SLURM_CPUS_PER_TASK" -o "$OUT_FASTQC" "$R1" "$R2"

echo "[INFO] Done: $SAMPLE"
