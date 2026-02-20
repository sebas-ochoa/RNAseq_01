#!/usr/bin/env bash
#SBATCH --job-name=fastqc_trim
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

PROJECT_ROOT="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is empty}"
cd "$PROJECT_ROOT"

TRIM_DIR="${PROJECT_ROOT}/02_work/02_trim/fastp"
OUT_FASTQC="${PROJECT_ROOT}/02_work/03_qc_trim/fastqc"
mkdir -p "$OUT_FASTQC"

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${PROJECT_ROOT}/00_admin/metadata/samples.tsv")
if [[ -z "$LINE" ]]; then
  echo "[WARN] No sample found for task ${SLURM_ARRAY_TASK_ID}; exiting."
  exit 0
fi

SAMPLE=$(echo "$LINE" | cut -f1)
R1="${TRIM_DIR}/${SAMPLE}_1.trim.fq.gz"
R2="${TRIM_DIR}/${SAMPLE}_2.trim.fq.gz"

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
  echo "[ERROR] Trimmed FASTQ not found for sample: $SAMPLE"
  echo "[ERROR] Expected: $R1"
  echo "[ERROR] Expected: $R2"
  exit 1
fi

fastqc -t "$SLURM_CPUS_PER_TASK" -o "$OUT_FASTQC" "$R1" "$R2"

echo "[INFO] Done: $SAMPLE"
