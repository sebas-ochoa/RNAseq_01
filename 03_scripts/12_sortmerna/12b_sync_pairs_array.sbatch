#!/usr/bin/env bash
#SBATCH --job-name=12b_sync
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=01:00:00
#SBATCH --output=04_logs/slurm/%x_%A_%a.out
#SBATCH --error=04_logs/slurm/%x_%A_%a.err
#SBATCH --array=1-12

set -euo pipefail
umask 007

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"
source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

module load "${SEQTK_MODULE}"

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES_TSV")
[[ -n "$LINE" ]] || exit 0
SAMPLE=$(echo "$LINE" | cut -f1)

FWD="${SORTMERNA_DIR}/${SAMPLE}_fwd.fq.gz"
REV="${SORTMERNA_DIR}/${SAMPLE}_rev.fq.gz"

echo "[INFO] Sample: $SAMPLE"

SYNCDIR="${TMP_ROOT}/sync_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$SYNCDIR"

# Extract read names from rev (the reference — smaller due to sortmerna v4 bug)
zcat "$REV" | awk 'NR%4==1{print substr($1,2)}' > "${SYNCDIR}/names.txt"
N_REV=$(wc -l < "${SYNCDIR}/names.txt")
echo "[INFO] Rev reads: $N_REV"

# Filter fwd to keep only reads present in rev (maintains input order)
seqtk subseq "$FWD" "${SYNCDIR}/names.txt" | gzip -c > "${FWD}.tmp"
mv "${FWD}.tmp" "$FWD"
chmod u+rw,g+rw "$FWD"

N_FWD=$(zcat "$FWD" | awk 'END{print NR/4}')
echo "[INFO] After sync — fwd: $N_FWD  rev: $N_REV"
[[ "$N_FWD" -eq "$N_REV" ]] && echo "[INFO] OK: synchronized" || echo "[WARN] Still unequal"

rm -rf "$SYNCDIR"
