#!/usr/bin/env bash
#SBATCH --job-name=15_subsample
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=4G
#SBATCH --time=10:00:00
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
module load "$SEQTK_MODULE"

mkdir -p "$SUB15_DIR" "04_logs/slurm"
chmod -R u+rwX,g+rwX "$SUB15_DIR" 2>/dev/null || true

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES_TSV")
[[ -n "$LINE" ]] || exit 0
SAMPLE=$(echo "$LINE" | cut -f1)

R1="${TRIM14_DIR}/${SAMPLE}_1.trim.fq.gz"
R2="${TRIM14_DIR}/${SAMPLE}_2.trim.fq.gz"
OUT_R1="${SUB15_DIR}/${SAMPLE}_1.sub.fq.gz"
OUT_R2="${SUB15_DIR}/${SAMPLE}_2.sub.fq.gz"

# Count reads in R1 (lines / 4)
N_READS=$(zcat "$R1" | awk 'END{print int(NR/4)}')
echo "[INFO] Sample: $SAMPLE — $N_READS reads, target ${SUBSAMPLE_READS}"

if [[ "$N_READS" -le "$SUBSAMPLE_READS" ]]; then
  echo "[INFO] Already at or below target — copying as-is"
  cp "$R1" "$OUT_R1"
  cp "$R2" "$OUT_R2"
else
  # Use fraction instead of absolute count (constant memory)
  FRAC=$(awk "BEGIN{printf \"%.6f\", ${SUBSAMPLE_READS}/${N_READS}}")
  echo "[INFO] Subsampling fraction: $FRAC"
  seqtk sample -s 100 "$R1" "$FRAC" | gzip -c > "$OUT_R1"
  seqtk sample -s 100 "$R2" "$FRAC" | gzip -c > "$OUT_R2"
fi

N_OUT=$(zcat "$OUT_R1" | awk 'END{print int(NR/4)}')
echo "[INFO] Done: $SAMPLE — output $N_OUT reads"
