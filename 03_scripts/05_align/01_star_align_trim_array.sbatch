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

on_error() {
  local line="$1"
  local code="$2"
  echo "[ERROR] Task ${SLURM_ARRAY_TASK_ID:-NA} failed at line ${line}: ${BASH_COMMAND}" >&2
  exit "$code"
}
trap 'on_error $LINENO $?' ERR

ensure_writable_dir() {
  local dir="$1"
  mkdir -p "$dir"
  chmod 2770 "$dir" 2>/dev/null || true
  if [[ ! -d "$dir" || ! -w "$dir" ]]; then
    echo "[ERROR] Directory missing or not writable: $dir" >&2
    ls -ld "$dir" >&2 || true
    exit 1
  fi
  local probe="${dir}/.write_test_${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-NA}"
  if ! touch "$probe" 2>/dev/null; then
    echo "[ERROR] Could not create file in: $dir" >&2
    ls -ld "$dir" >&2 || true
    exit 1
  fi
  rm -f "$probe"
}

PROJECT_ROOT="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is empty}"
cd "$PROJECT_ROOT"
source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

module purge
module load "$STACK_MODULE"
module load "$GCC_MODULE"
module load "$STAR_MODULE"
module load "$SAMTOOLS_MODULE"

ensure_writable_dir "$SLURM_LOG_DIR"
ensure_writable_dir "$ALIGN_DIR"
ensure_writable_dir "$TMP_ROOT"

if [[ ! -d "$STAR_INDEX_DIR" || ! -s "${STAR_INDEX_DIR}/SA" ]]; then
  echo "[ERROR] STAR index not found/complete in: $STAR_INDEX_DIR" >&2
  echo "[ERROR] Run 03_scripts/04_ref/02_star_index.sbatch first." >&2
  exit 1
fi

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES_TSV")
if [[ -z "$LINE" ]]; then
  echo "[WARN] No sample found for task ${SLURM_ARRAY_TASK_ID}; exiting."
  exit 0
fi

SAMPLE=$(echo "$LINE" | cut -f1)
R1="${TRIM_DIR}/${SAMPLE}_1.trim.fq.gz"
R2="${TRIM_DIR}/${SAMPLE}_2.trim.fq.gz"

if [[ ! -s "$R1" || ! -s "$R2" ]]; then
  echo "[ERROR] Missing trimmed FASTQ for sample: $SAMPLE" >&2
  echo "[ERROR] Expected: $R1" >&2
  echo "[ERROR] Expected: $R2" >&2
  exit 1
fi

FINAL_BAM="${ALIGN_DIR}/${SAMPLE}.sorted.bam"
FINAL_BAI="${FINAL_BAM}.bai"
if [[ -s "$FINAL_BAM" && -s "$FINAL_BAI" ]]; then
  echo "[INFO] BAM + BAI already exist for $SAMPLE; skipping alignment."
  exit 0
fi

JOB_TMP="${TMP_ROOT}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
ensure_writable_dir "$JOB_TMP"

echo "[INFO] Sample: $SAMPLE"
echo "[INFO] R1: $R1"
echo "[INFO] R2: $R2"

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

STAR_BAM="${ALIGN_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam"
if [[ ! -s "$STAR_BAM" ]]; then
  echo "[ERROR] STAR output BAM missing for sample: $SAMPLE" >&2
  exit 1
fi

mv -f "$STAR_BAM" "$FINAL_BAM"
samtools index -@ "$SLURM_CPUS_PER_TASK" "$FINAL_BAM"

echo "[INFO] Done: $SAMPLE"
echo "[INFO] BAM: $FINAL_BAM"
