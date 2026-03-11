#!/usr/bin/env bash
# Submit Qualimap BAM QC (array) + global MultiQC with SLURM dependency.
# Run from the project root: bash 03_scripts/00_helpers/02_submit_qualimap_multiqc.sh
set -euo pipefail
umask 007

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$PROJECT_ROOT"

source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

mkdir -p "$SLURM_LOG_DIR" "$TMP_ROOT"
chmod -R u+rwX,g+rwX "$SLURM_LOG_DIR" "$TMP_ROOT" 2>/dev/null || true

N_SAMPLES=$(awk 'NR>1 && NF>=1 {n++} END{print n+0}' "$SAMPLES_TSV")
[[ "$N_SAMPLES" -gt 0 ]] || { echo "[ERROR] No samples found in $SAMPLES_TSV"; exit 1; }

J_QUALIMAP=$(sbatch --parsable --array="1-${N_SAMPLES}" \
  "${PROJECT_ROOT}/03_scripts/05_align/04_qualimap_bamqc_array.sbatch")

J_MULTIQC_ALL=$(sbatch --parsable \
  --dependency="afterok:${J_QUALIMAP}" \
  "${PROJECT_ROOT}/03_scripts/05_align/05_multiqc_full.sbatch")

printf "qualimap_bamqc_array\t%s\n" "$J_QUALIMAP"
printf "multiqc_full\t%s\n" "$J_MULTIQC_ALL"
echo "squeue -u $USER"
echo "sacct -j ${J_MULTIQC_ALL} --format=JobID,JobName%30,State,ExitCode,Elapsed,MaxRSS,NCPUS"
