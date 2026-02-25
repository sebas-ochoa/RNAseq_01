#!/usr/bin/env bash
#SBATCH --job-name=fastqc_trim
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

TRIM_DIR="${PROJECT_ROOT}/02_work/02_trim/fastp"
OUT_BASE="${PROJECT_ROOT}/02_work/03_qc_trim"
OUT_FASTQC="${OUT_BASE}/fastqc"
mkdir -p "$OUT_BASE" "$OUT_FASTQC" "${PROJECT_ROOT}/04_logs/slurm"
chmod -R u+rwX,g+rwX "$OUT_BASE" "$OUT_FASTQC" "${PROJECT_ROOT}/04_logs/slurm" 2>/dev/null || true

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${PROJECT_ROOT}/00_admin/metadata/samples.tsv")
[[ -n "$LINE" ]] || exit 0

SAMPLE=$(echo "$LINE" | cut -f1)
R1="${TRIM_DIR}/${SAMPLE}_1.trim.fq.gz"
R2="${TRIM_DIR}/${SAMPLE}_2.trim.fq.gz"

echo "[INFO] Sample: $SAMPLE"
fastqc -t "$SLURM_CPUS_PER_TASK" -o "$OUT_FASTQC" "$R1" "$R2"

echo "[INFO] Done: $SAMPLE"
