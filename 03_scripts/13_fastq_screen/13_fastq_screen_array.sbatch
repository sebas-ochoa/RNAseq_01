#!/usr/bin/env bash
#SBATCH --job-name=13_fastq_screen
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
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
module load "$BOWTIE2_MODULE"

mkdir -p "$FASTQ_SCREEN_DIR" "04_logs/slurm"
chmod -R u+rwX,g+rwX "$FASTQ_SCREEN_DIR" 2>/dev/null || true

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES_TSV")
[[ -n "$LINE" ]] || exit 0
SAMPLE=$(echo "$LINE" | cut -f1)

R1="${SORTMERNA_DIR}/${SAMPLE}_fwd.fq.gz"

echo "[INFO] Sample: $SAMPLE"

"${FASTQ_SCREEN_HOME}/fastq_screen" \
  --conf "$FASTQ_SCREEN_CONF" \
  --aligner bowtie2 \
  --threads "$SLURM_CPUS_PER_TASK" \
  --outdir "$FASTQ_SCREEN_DIR" \
  "$R1"

echo "[INFO] Done: $SAMPLE"
