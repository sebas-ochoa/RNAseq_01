#!/usr/bin/env bash
#SBATCH --job-name=multiqc_raw
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:20:00
#SBATCH --output=04_logs/slurm/%x_%j.out
#SBATCH --error=04_logs/slurm/%x_%j.err

set -euo pipefail
source 00_admin/env/modules.sh

OUT_FASTQC="02_work/01_qc_raw/fastqc"
OUT_MULTI="02_work/01_qc_raw/multiqc"
mkdir -p "$OUT_MULTI"

multiqc "$OUT_FASTQC" -o "$OUT_MULTI"
echo "[INFO] Wrote: $OUT_MULTI/multiqc_report.html"
