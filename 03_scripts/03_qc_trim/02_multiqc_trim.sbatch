#!/usr/bin/env bash
#SBATCH --job-name=multiqc_trim
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:30:00
#SBATCH --output=04_logs/slurm/%x_%j.out
#SBATCH --error=04_logs/slurm/%x_%j.err

set -euo pipefail

module purge
module load stack/2024-06
module load python/3.11.6

cd "$SLURM_SUBMIT_DIR"

# Activate your user-installed MultiQC
source 00_admin/env/multiqc_venv/bin/activate

# Force a writable temp directory (avoids /scratch permission issues)
TMPDIR="$SLURM_SUBMIT_DIR/90_scratch/tmp/${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
export TMPDIR
export TEMP="$TMPDIR"
export TMP="$TMPDIR"

IN_DIR="02_work/03_qc_trim/fastqc"
OUT_DIR="02_work/03_qc_trim/multiqc"
mkdir -p "$OUT_DIR"

if ! compgen -G "${IN_DIR}/*_fastqc.zip" > /dev/null; then
  echo "[ERROR] No FastQC zip files found in ${IN_DIR}."
  echo "[ERROR] Run 03_scripts/03_qc_trim/01_fastqc_trim_array.sbatch first."
  deactivate
  exit 1
fi

echo "[INFO] TMPDIR: $TMPDIR"
echo "[INFO] Running MultiQC..."
multiqc "$IN_DIR" -o "$OUT_DIR" --filename multiqc_report.html

echo "[INFO] Done: $OUT_DIR/multiqc_report.html"
deactivate
