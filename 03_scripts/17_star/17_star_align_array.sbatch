#!/usr/bin/env bash
#SBATCH --job-name=17_star
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

mkdir -p "$ALIGN17_DIR" "$TMP_ROOT" "04_logs/slurm"
chmod -R u+rwX,g+rwX "$ALIGN17_DIR" "$TMP_ROOT" 2>/dev/null || true

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES_TSV")
[[ -n "$LINE" ]] || exit 0
SAMPLE=$(echo "$LINE" | cut -f1)

R1="${SUB15_DIR}/${SAMPLE}_1.sub.fq.gz"
R2="${SUB15_DIR}/${SAMPLE}_2.sub.fq.gz"
FINAL_BAM="${ALIGN17_DIR}/${SAMPLE}.sorted.bam"
JOB_TMP="${TMP_ROOT}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$JOB_TMP"

echo "[INFO] Sample: $SAMPLE"

STAR \
  --runThreadN "$SLURM_CPUS_PER_TASK" \
  --genomeDir "$STAR_INDEX_DIR" \
  --readFilesIn "$R1" "$R2" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${ALIGN17_DIR}/${SAMPLE}." \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM XS \
  --outTmpDir "${JOB_TMP}/star_${SAMPLE}" \
  --limitBAMsortRAM "$STAR_ALIGN_SORT_RAM"

mv -f "${ALIGN17_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam" "$FINAL_BAM"
samtools index -@ "$SLURM_CPUS_PER_TASK" "$FINAL_BAM"
rm -rf "$JOB_TMP"

echo "[INFO] Done: $FINAL_BAM"
