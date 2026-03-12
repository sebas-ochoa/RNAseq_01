#!/usr/bin/env bash
#SBATCH --job-name=12_sortmerna
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=08:00:00
#SBATCH --output=04_logs/slurm/%x_%A_%a.out
#SBATCH --error=04_logs/slurm/%x_%A_%a.err
#SBATCH --array=1-12

set -euo pipefail
umask 007

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"
source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

mkdir -p "$SORTMERNA_DIR" "$TMP_ROOT" "04_logs/slurm"
chmod -R u+rwX,g+rwX "$SORTMERNA_DIR" 2>/dev/null || true

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES_TSV")
[[ -n "$LINE" ]] || exit 0
SAMPLE=$(echo "$LINE" | cut -f1)
R1=$(echo "$LINE" | cut -f2)
R2=$(echo "$LINE" | cut -f3)

WORKDIR="${TMP_ROOT}/sortmerna_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$WORKDIR"

echo "[INFO] Sample: $SAMPLE"

"${SORTMERNA_HOME}/bin/sortmerna" \
  --ref "${SORTMERNA_DB_DIR}/silva-arc-16s-id95.fasta" \
  --ref "${SORTMERNA_DB_DIR}/silva-arc-23s-id98.fasta" \
  --ref "${SORTMERNA_DB_DIR}/silva-bac-16s-id90.fasta" \
  --ref "${SORTMERNA_DB_DIR}/silva-bac-23s-id98.fasta" \
  --ref "${SORTMERNA_DB_DIR}/silva-euk-18s-id95.fasta" \
  --ref "${SORTMERNA_DB_DIR}/silva-euk-28s-id98.fasta" \
  --ref "${SORTMERNA_DB_DIR}/rfam-5s-database-id98.fasta" \
  --ref "${SORTMERNA_DB_DIR}/rfam-5.8s-database-id98.fasta" \
  --reads "$R1" \
  --reads "$R2" \
  --fastx \
  --paired-in \
  --out2 \
  --aligned "${WORKDIR}/${SAMPLE}.rRNA" \
  --other   "${SORTMERNA_DIR}/${SAMPLE}" \
  --threads "$SLURM_CPUS_PER_TASK" \
  --workdir "$WORKDIR"

rm -rf "$WORKDIR"

echo "[INFO] Done: ${SORTMERNA_DIR}/${SAMPLE}_fwd.fq.gz"
