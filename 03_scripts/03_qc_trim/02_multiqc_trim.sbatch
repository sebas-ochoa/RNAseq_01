#!/usr/bin/env bash
#SBATCH --job-name=multiqc_trim
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:30:00
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

module purge
module load stack/2024-06
module load python/3.11.6

cd "$SLURM_SUBMIT_DIR"

IN_DIR="02_work/03_qc_trim/fastqc"
OUT_BASE="02_work/03_qc_trim"
OUT_DIR="${OUT_BASE}/multiqc"
TMPDIR="$SLURM_SUBMIT_DIR/90_scratch/tmp/${SLURM_JOB_ID}"

mkdir -p "$OUT_BASE" "$OUT_DIR" "$TMPDIR"
chmod 2770 "$OUT_BASE" "$OUT_DIR" "$TMPDIR" 2>/dev/null || true

for d in "$OUT_DIR" "$TMPDIR"; do
  if [[ ! -d "$d" || ! -w "$d" ]]; then
    echo "[ERROR] Directory missing or not writable: $d" >&2
    ls -ld "$OUT_BASE" "$OUT_DIR" "$TMPDIR" >&2 || true
    exit 1
  fi
  probe="${d}/.write_test_${SLURM_JOB_ID}"
  if ! touch "$probe" 2>/dev/null; then
    echo "[ERROR] Could not create file in: $d" >&2
    ls -ld "$OUT_BASE" "$OUT_DIR" "$TMPDIR" >&2 || true
    exit 1
  fi
  rm -f "$probe"
done

# Activate your user-installed MultiQC
source 00_admin/env/multiqc_venv/bin/activate

export TMPDIR
export TEMP="$TMPDIR"
export TMP="$TMPDIR"

if ! compgen -G "${IN_DIR}/*_fastqc.zip" > /dev/null; then
  echo "[ERROR] No FastQC zip files found in ${IN_DIR}." >&2
  echo "[ERROR] Run 03_scripts/03_qc_trim/01_fastqc_trim_array.sbatch first." >&2
  deactivate
  exit 1
fi

echo "[INFO] TMPDIR: $TMPDIR"
echo "[INFO] Running MultiQC..."
multiqc "$IN_DIR" -o "$OUT_DIR" --filename multiqc_report.html

echo "[INFO] Done: $OUT_DIR/multiqc_report.html"
deactivate
