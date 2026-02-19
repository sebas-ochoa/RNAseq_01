#!/usr/bin/env bash
#SBATCH --job-name=trim_fastp
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=02:00:00
#SBATCH --output=04_logs/slurm/%x_%A_%a.out
#SBATCH --error=04_logs/slurm/%x_%A_%a.err
#SBATCH --array=1-12

set -euo pipefail

module purge
module load stack/2024-06
module load fastp/0.23.4

cd "$SLURM_SUBMIT_DIR"

OUT_DIR="02_work/02_trim/fastp"
mkdir -p "$OUT_DIR"

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" 00_admin/metadata/samples.tsv)
SAMPLE=$(echo "$LINE" | cut -f1)
R1=$(echo "$LINE" | cut -f2)
R2=$(echo "$LINE" | cut -f3)

OUT_R1="${OUT_DIR}/${SAMPLE}_1.trim.fq.gz"
OUT_R2="${OUT_DIR}/${SAMPLE}_2.trim.fq.gz"
HTML="${OUT_DIR}/${SAMPLE}.fastp.html"

echo "[INFO] Sample: $SAMPLE"
echo "[INFO] R1: $R1"
echo "[INFO] R2: $R2"

fastp \
  -i "$R1" -I "$R2" \
  -o "$OUT_R1" -O "$OUT_R2" \
  --detect_adapter_for_pe \
  --trim_front1 8 --trim_front2 8 \
  --cut_tail \
  --cut_window_size 4 --cut_mean_quality 20 \
  --length_required 25 \
  --thread "$SLURM_CPUS_PER_TASK" \
  --html "$HTML"

echo "[INFO] Done: $SAMPLE"
