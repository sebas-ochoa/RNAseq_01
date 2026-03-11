#!/usr/bin/env bash
#SBATCH --job-name=multiqc_full
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:45:00
#SBATCH --output=04_logs/slurm/%x_%j.out
#SBATCH --error=04_logs/slurm/%x_%j.err

set -euo pipefail
umask 007

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"
source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

MULTIQC_ALL_DIR="${MULTIQC_ALL_DIR:-${PROJECT_ROOT}/02_work/multiqc_all}"

module purge
module load "$STACK_MODULE"
module load "$PYTHON_MODULE"

mkdir -p "$SLURM_LOG_DIR" "$MULTIQC_ALL_DIR" "$TMP_ROOT"
chmod -R u+rwX,g+rwX "$SLURM_LOG_DIR" "$MULTIQC_ALL_DIR" "$TMP_ROOT" 2>/dev/null || true

source "$MULTIQC_VENV/bin/activate"

TMPDIR="${TMP_ROOT}/${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
chmod -R u+rwX,g+rwX "$TMPDIR" 2>/dev/null || true
export TMPDIR TEMP="$TMPDIR" TMP="$TMPDIR"

multiqc 02_work/ \
  -o "$MULTIQC_ALL_DIR" \
  --filename multiqc_report.html \
  --ignore "*/star_index/*" \
  -f

deactivate
echo "[INFO] Done: ${MULTIQC_ALL_DIR}/multiqc_report.html"
