#!/usr/bin/env bash
#SBATCH --job-name=multiqc_align
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:45:00
#SBATCH --output=04_logs/slurm/%x_%j.out
#SBATCH --error=04_logs/slurm/%x_%j.err

set -euo pipefail
umask 007

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"
source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

module purge
module load "$STACK_MODULE"
module load "$PYTHON_MODULE"

mkdir -p "$SLURM_LOG_DIR" "$ALIGN_QC_DIR" "$TMP_ROOT"
chmod -R u+rwX,g+rwX "$SLURM_LOG_DIR" "$ALIGN_QC_DIR" "$TMP_ROOT" 2>/dev/null || true

source "$MULTIQC_VENV/bin/activate"

TMPDIR="${TMP_ROOT}/${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
chmod -R u+rwX,g+rwX "$TMPDIR" 2>/dev/null || true
export TMPDIR
export TEMP="$TMPDIR"
export TMP="$TMPDIR"

multiqc "$ALIGN_DIR" -o "$ALIGN_QC_DIR" --filename multiqc_report.html

deactivate
echo "[INFO] Done: $ALIGN_QC_DIR/multiqc_report.html"
