#!/usr/bin/env bash
#SBATCH --job-name=multiqc_align
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:45:00
#SBATCH --output=04_logs/slurm/%x_%j.out
#SBATCH --error=04_logs/slurm/%x_%j.err

set -euo pipefail
umask 007

on_error() {
  local line="$1"
  local code="$2"
  echo "[ERROR] Job ${SLURM_JOB_ID:-NA} failed at line ${line}: ${BASH_COMMAND}" >&2
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
  local probe="${dir}/.write_test_${SLURM_JOB_ID:-manual}"
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
module load "$PYTHON_MODULE"

ensure_writable_dir "$SLURM_LOG_DIR"
ensure_writable_dir "$ALIGN_QC_DIR"
ensure_writable_dir "$TMP_ROOT"

if [[ ! -f "$MULTIQC_VENV/bin/activate" ]]; then
  echo "[ERROR] MultiQC venv not found: $MULTIQC_VENV" >&2
  exit 1
fi

if ! compgen -G "${ALIGN_DIR}/*.Log.final.out" > /dev/null; then
  echo "[ERROR] No STAR log files found in ${ALIGN_DIR}." >&2
  echo "[ERROR] Run 03_scripts/05_align/01_star_align_trim_array.sbatch first." >&2
  exit 1
fi

source "$MULTIQC_VENV/bin/activate"

TMPDIR="${TMP_ROOT}/${SLURM_JOB_ID}"
ensure_writable_dir "$TMPDIR"
export TMPDIR
export TEMP="$TMPDIR"
export TMP="$TMPDIR"

echo "[INFO] Running MultiQC on STAR outputs..."
multiqc "$ALIGN_DIR" -o "$ALIGN_QC_DIR" --filename multiqc_report.html

echo "[INFO] Done: $ALIGN_QC_DIR/multiqc_report.html"
deactivate
