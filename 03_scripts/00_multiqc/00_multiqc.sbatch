#!/usr/bin/env bash
#SBATCH --job-name=00_multiqc
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

module purge
module load "$STACK_MODULE"
module load "$PYTHON_MODULE"

mkdir -p "$MULTIQC00_DIR" "$TMP_ROOT" "04_logs/slurm"
chmod -R u+rwX,g+rwX "$MULTIQC00_DIR" 2>/dev/null || true

source "${MULTIQC_VENV}/bin/activate"

TMPDIR="${TMP_ROOT}/${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
export TMPDIR TEMP="$TMPDIR" TMP="$TMPDIR"

multiqc 02_work/ \
  -o "$MULTIQC00_DIR" \
  --filename multiqc_report.html \
  --ignore "*/star_index/*" \
  -f

deactivate
echo "[INFO] Done: ${MULTIQC00_DIR}/multiqc_report.html"
