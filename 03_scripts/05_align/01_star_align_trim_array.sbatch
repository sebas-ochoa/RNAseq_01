#!/usr/bin/env bash
#SBATCH --job-name=star_align_trim
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

module purge
module load "$STACK_MODULE"
module load "$GCC_MODULE"
module load "$STAR_MODULE"
module load "$SAMTOOLS_MODULE"

mkdir -p "$SLURM_LOG_DIR" "$ALIGN_DIR" "$TMP_ROOT"
chmod -R u+rwX,g+rwX "$SLURM_LOG_DIR" "$ALIGN_DIR" "$TMP_ROOT" 2>/dev/null || true

[[ -s "${STAR_INDEX_DIR}/SA" ]] || { echo "[ERROR] STAR index missing in $STAR_INDEX_DIR"; exit 1; }

LINE="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES_TSV")"
[[ -n "$LINE" ]] || exit 0

SAMPLE="$(echo "$LINE" | cut -f1)"
R1="${TRIM_DIR}/${SAMPLE}_1.trim.fq.gz"
R2="${TRIM_DIR}/${SAMPLE}_2.trim.fq.gz"

FINAL_BAM="${ALIGN_DIR}/${SAMPLE}.sorted.bam"
FINAL_BAI="${FINAL_BAM}.bai"
[[ -s "$FINAL_BAM" && -s "$FINAL_BAI" ]] && exit 0

JOB_TMP="${TMP_ROOT}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$JOB_TMP"
chmod -R u+rwX,g+rwX "$JOB_TMP" 2>/dev/null || true

STAR \
  --runThreadN "$SLURM_CPUS_PER_TASK" \
  --genomeDir "$STAR_INDEX_DIR" \
  --readFilesIn "$R1" "$R2" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${ALIGN_DIR}/${SAMPLE}." \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM XS \
  --outTmpDir "${JOB_TMP}/star_${SAMPLE}" \
  --limitBAMsortRAM "$STAR_ALIGN_SORT_RAM"

mv -f "${ALIGN_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam" "$FINAL_BAM"
samtools index -@ "$SLURM_CPUS_PER_TASK" "$FINAL_BAM"
rm -rf "$JOB_TMP" 2>/dev/null || true

echo "[INFO] Done: $SAMPLE"
