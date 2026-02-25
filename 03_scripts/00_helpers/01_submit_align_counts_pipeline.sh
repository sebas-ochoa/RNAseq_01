#!/usr/bin/env bash
set -euo pipefail
umask 007

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$PROJECT_ROOT"

source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

mkdir -p "$SLURM_LOG_DIR" "$TMP_ROOT" "$ALIGN_DIR" "$ALIGN_QC_DIR" "$COUNTS_DIR"
chmod -R u+rwX,g+rwX "$SLURM_LOG_DIR" "$TMP_ROOT" "$ALIGN_DIR" "$ALIGN_QC_DIR" "$COUNTS_DIR" 2>/dev/null || true

N_SAMPLES=$(awk 'NR>1 && NF>=1 {n++} END{print n+0}' "$SAMPLES_TSV")

J_ALIGN=$(sbatch --parsable --array="1-${N_SAMPLES}" "${PROJECT_ROOT}/03_scripts/05_align/01_star_align_trim_array.sbatch")
J_ALIGN_QC=$(sbatch --parsable --dependency="afterok:${J_ALIGN}" "${PROJECT_ROOT}/03_scripts/05_align/02_align_qc_multiqc.sbatch")
J_FINAL=$(sbatch --parsable --dependency="afterok:${J_ALIGN}" "${PROJECT_ROOT}/03_scripts/06_counts/02_featurecounts_final.sbatch")

printf "star_align_array\t%s\n" "$J_ALIGN"
printf "align_multiqc\t%s\n" "$J_ALIGN_QC"
printf "final_featurecounts\t%s\n" "$J_FINAL"
echo "squeue -u $USER"
echo "sacct -j ${J_FINAL} --format=JobID,JobName%30,State,ExitCode,Elapsed,MaxRSS,NCPUS"
