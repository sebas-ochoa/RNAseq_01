#!/usr/bin/env bash
#SBATCH --job-name=14_fastp
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=02:00:00
#SBATCH --output=04_logs/slurm/%x_%A_%a.out
#SBATCH --error=04_logs/slurm/%x_%A_%a.err
#SBATCH --array=1-12

set -euo pipefail
umask 007

module purge
module load stack/2024-06
module load fastp/0.23.4

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"
source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

mkdir -p "$TRIM14_DIR" "04_logs/slurm"
chmod -R u+rwX,g+rwX "$TRIM14_DIR" 2>/dev/null || true

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES_TSV")
[[ -n "$LINE" ]] || exit 0
SAMPLE=$(echo "$LINE" | cut -f1)

R1="${SORTMERNA_DIR}/${SAMPLE}_fwd.fq.gz"
R2="${SORTMERNA_DIR}/${SAMPLE}_rev.fq.gz"
OUT_R1="${TRIM14_DIR}/${SAMPLE}_1.trim.fq.gz"
OUT_R2="${TRIM14_DIR}/${SAMPLE}_2.trim.fq.gz"
HTML="${TRIM14_DIR}/${SAMPLE}.fastp.html"

echo "[INFO] Sample: $SAMPLE"

fastp \
  -i "$R1" -I "$R2" \
  -o "$OUT_R1" -O "$OUT_R2" \
  --detect_adapter_for_pe \
  --cut_tail \
  --cut_window_size 4 --cut_mean_quality 20 \
  --length_required 25 \
  --thread "$SLURM_CPUS_PER_TASK" \
  --html "$HTML"

echo "[INFO] Done: $SAMPLE"
