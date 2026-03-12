#!/usr/bin/env bash
# Submit the full 1*_ pipeline with SLURM dependencies.
# Run from project root: bash 03_scripts/00_helpers/10_submit_pipeline.sh
set -euo pipefail
umask 007

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$PROJECT_ROOT"

source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

N=$(awk 'NR>1 && NF>=1 {n++} END{print n+0}' "$SAMPLES_TSV")
[[ "$N" -gt 0 ]] || { echo "[ERROR] No samples in $SAMPLES_TSV"; exit 1; }
echo "[INFO] Submitting pipeline for $N samples"

S="${PROJECT_ROOT}/03_scripts"

J_QC=$(sbatch  --parsable --array="1-${N}" "${S}/11_qc_raw/11_fastqc_array.sbatch")
J_RNA=$(sbatch --parsable --array="1-${N}" "${S}/12_sortmerna/12_sortmerna_array.sbatch")
J_FQS=$(sbatch --parsable --array="1-${N}" --dependency="afterok:${J_RNA}" "${S}/13_fastq_screen/13_fastq_screen_array.sbatch")
J_TRM=$(sbatch --parsable --array="1-${N}" --dependency="afterok:${J_RNA}" "${S}/14_trim/14_fastp_array.sbatch")
J_ALN=$(sbatch --parsable --array="1-${N}" --dependency="afterok:${J_TRM}" "${S}/15_star/15_star_align_array.sbatch")
J_CNT=$(sbatch --parsable                  --dependency="afterok:${J_ALN}" "${S}/16_featurecounts/16_featurecounts.sbatch")
J_QMP=$(sbatch --parsable --array="1-${N}" --dependency="afterok:${J_ALN}" "${S}/17_qualimap/17_qualimap_array.sbatch")
J_MQC=$(sbatch --parsable                  --dependency="afterok:${J_CNT}:${J_QMP}" "${S}/00_multiqc/00_multiqc.sbatch")

printf "11_fastqc         %s\n" "$J_QC"
printf "12_sortmerna      %s\n" "$J_RNA"
printf "13_fastq_screen   %s  (dep: %s)\n" "$J_FQS" "$J_RNA"
printf "14_fastp          %s  (dep: %s)\n" "$J_TRM" "$J_RNA"
printf "15_star           %s  (dep: %s)\n" "$J_ALN" "$J_TRM"
printf "16_featurecounts  %s  (dep: %s)\n" "$J_CNT" "$J_ALN"
printf "17_qualimap       %s  (dep: %s)\n" "$J_QMP" "$J_ALN"
printf "00_multiqc        %s  (dep: %s:%s)\n" "$J_MQC" "$J_CNT" "$J_QMP"
echo ""
echo "squeue -u $USER"
