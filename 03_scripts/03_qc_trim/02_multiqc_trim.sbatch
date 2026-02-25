#!/usr/bin/env bash
#SBATCH --job-name=multiqc_trim
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:30:00
#SBATCH --output=04_logs/slurm/%x_%j.out
#SBATCH --error=04_logs/slurm/%x_%j.err

set -euo pipefail
umask 007

module purge
module load stack/2024-06
module load python/3.11.6

cd "$SLURM_SUBMIT_DIR"

IN_DIR="02_work/03_qc_trim/fastqc"
OUT_BASE="02_work/03_qc_trim"
OUT_DIR="${OUT_BASE}/multiqc"
TMPDIR="$SLURM_SUBMIT_DIR/90_scratch/tmp/${SLURM_JOB_ID}"

mkdir -p "$OUT_BASE" "$OUT_DIR" "$TMPDIR" "04_logs/slurm"
chmod -R u+rwX,g+rwX "$OUT_BASE" "$OUT_DIR" "$TMPDIR" "04_logs/slurm" 2>/dev/null || true

source 00_admin/env/multiqc_venv/bin/activate

export TMPDIR
export TEMP="$TMPDIR"
export TMP="$TMPDIR"

echo "[INFO] TMPDIR: $TMPDIR"
echo "[INFO] Running MultiQC..."
multiqc "$IN_DIR" -o "$OUT_DIR" --filename multiqc_report.html

echo "[INFO] Done: $OUT_DIR/multiqc_report.html"
deactivate
