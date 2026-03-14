#!/usr/bin/env bash
# Resume pipeline from pair-sync step.
# Use when SortMeRNA already ran but pairs need syncing before fastp.
# HSF4_2H_2 (array 5) already re-ran with sync built-in — excluded from J_SYNC.
# Run from project root: bash 03_scripts/00_helpers/11_resume_from_sync.sh
set -euo pipefail
umask 007

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$PROJECT_ROOT"

source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

N=$(awk 'NR>1 && NF>=1 {n++} END{print n+0}' "$SAMPLES_TSV")
S="${PROJECT_ROOT}/03_scripts"

# Sync pairs for all 12 samples
J_SYNC=$(sbatch --parsable --array="1-${N}" "${S}/12_sortmerna/12b_sync_pairs_array.sbatch")

# FastQ Screen and fastp both wait for sync to finish
J_FQS=$(sbatch --parsable --array="1-${N}" --dependency="afterok:${J_SYNC}" "${S}/13_fastq_screen/13_fastq_screen_array.sbatch")
J_TRM=$(sbatch --parsable --array="1-${N}" --dependency="afterok:${J_SYNC}" "${S}/14_trim/14_fastp_array.sbatch")
J_SUB=$(sbatch --parsable --array="1-${N}" --dependency="afterok:${J_TRM}"  "${S}/15_subsample/15_subsample_array.sbatch")
J_QC2=$(sbatch --parsable --array="1-${N}" --dependency="afterok:${J_SUB}"  "${S}/16_fastqc/16_fastqc_array.sbatch")
J_ALN=$(sbatch --parsable --array="1-${N}" --dependency="afterok:${J_SUB}"  "${S}/17_star/17_star_align_array.sbatch")
J_CNT=$(sbatch --parsable                  --dependency="afterok:${J_ALN}"  "${S}/18_featurecounts/18_featurecounts.sbatch")
J_QMP=$(sbatch --parsable --array="1-${N}" --dependency="afterok:${J_ALN}"  "${S}/19_qualimap/19_qualimap_array.sbatch")
J_MQC=$(sbatch --parsable                  --dependency="afterok:${J_CNT}:${J_QMP}" "${S}/00_multiqc/00_multiqc.sbatch")

printf "12b_sync          %s  (arrays 1-%s)\n" "$J_SYNC" "$N"
printf "13_fastq_screen   %s  (dep: %s)\n" "$J_FQS"  "$J_SYNC"
printf "14_fastp          %s  (dep: %s)\n" "$J_TRM"  "$J_SYNC"
printf "15_subsample      %s  (dep: %s)\n" "$J_SUB"  "$J_TRM"
printf "16_fastqc_sub     %s  (dep: %s)\n" "$J_QC2"  "$J_SUB"
printf "17_star           %s  (dep: %s)\n" "$J_ALN"  "$J_SUB"
printf "18_featurecounts  %s  (dep: %s)\n" "$J_CNT"  "$J_ALN"
printf "19_qualimap       %s  (dep: %s)\n" "$J_QMP"  "$J_ALN"
printf "00_multiqc        %s  (dep: %s:%s)\n" "$J_MQC" "$J_CNT" "$J_QMP"
echo ""
echo "squeue -u $USER"
