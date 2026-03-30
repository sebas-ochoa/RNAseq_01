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

module load stack/2024-06
module load "${SEQTK_MODULE}"

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES_TSV")
[[ -n "$LINE" ]] || exit 0
SAMPLE=$(echo "$LINE" | cut -f1)

FWD="${SORTMERNA_DIR}/${SAMPLE}_fwd.fq.gz"
REV="${SORTMERNA_DIR}/${SAMPLE}_rev.fq.gz"

echo "[INFO] Sample: $SAMPLE"

SYNCDIR="${TMP_ROOT}/sync_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$SYNCDIR"

# Extract read names from both files
zcat "$FWD" | awk 'NR%4==1{print substr($1,2)}' > "${SYNCDIR}/fwd_names.txt"
zcat "$REV" | awk 'NR%4==1{print substr($1,2)}' > "${SYNCDIR}/rev_names.txt"

N_FWD_RAW=$(wc -l < "${SYNCDIR}/fwd_names.txt")
N_REV_RAW=$(wc -l < "${SYNCDIR}/rev_names.txt")
echo "[INFO] Pre-sync — fwd: $N_FWD_RAW  rev: $N_REV_RAW"

if [[ "$N_FWD_RAW" -eq "$N_REV_RAW" ]]; then
  echo "[INFO] Already equal — skipping sync"
else
  # Bidirectional: keep only names present in BOTH files
  awk 'NR==FNR{a[$0]; next} $0 in a' \
    "${SYNCDIR}/fwd_names.txt" "${SYNCDIR}/rev_names.txt" > "${SYNCDIR}/common.txt"
  N_COMMON=$(wc -l < "${SYNCDIR}/common.txt")
  echo "[INFO] Common reads: $N_COMMON"

  seqtk subseq "$FWD" "${SYNCDIR}/common.txt" | gzip -c > "${FWD}.tmp"
  seqtk subseq "$REV" "${SYNCDIR}/common.txt" | gzip -c > "${REV}.tmp"
  mv "${FWD}.tmp" "$FWD"
  mv "${REV}.tmp" "$REV"
  chmod u+rw,g+rw "$FWD" "$REV"

  # Verify
  N_FWD_FINAL=$(zcat "$FWD" | awk 'END{print NR/4}')
  N_REV_FINAL=$(zcat "$REV" | awk 'END{print NR/4}')
  echo "[INFO] After sync — fwd: $N_FWD_FINAL  rev: $N_REV_FINAL"
  [[ "$N_FWD_FINAL" -eq "$N_REV_FINAL" ]] && echo "[INFO] OK: synchronized" || echo "[WARN] Still unequal"
fi

rm -rf "$SYNCDIR"
