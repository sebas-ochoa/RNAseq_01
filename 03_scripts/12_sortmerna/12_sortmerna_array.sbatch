#!/usr/bin/env bash
#SBATCH --job-name=12_sortmerna
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=08:00:00
#SBATCH --output=04_logs/slurm/%x_%A_%a.out
#SBATCH --error=04_logs/slurm/%x_%A_%a.err
#SBATCH --array=1-12

set -euo pipefail
umask 007

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"
source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

mkdir -p "$SORTMERNA_DIR" "$TMP_ROOT" "04_logs/slurm"
chmod -R u+rwX,g+rwX "$SORTMERNA_DIR" 2>/dev/null || true

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES_TSV")
[[ -n "$LINE" ]] || exit 0
SAMPLE=$(echo "$LINE" | cut -f1)
R1=$(echo "$LINE" | cut -f2)
R2=$(echo "$LINE" | cut -f3)

WORKDIR="${TMP_ROOT}/sortmerna_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$WORKDIR"

echo "[INFO] Sample: $SAMPLE"

"${SORTMERNA_HOME}/bin/sortmerna" \
  --ref "${SORTMERNA_DB_DIR}/silva-arc-16s-id95.fasta" \
  --ref "${SORTMERNA_DB_DIR}/silva-arc-23s-id98.fasta" \
  --ref "${SORTMERNA_DB_DIR}/silva-bac-16s-id90.fasta" \
  --ref "${SORTMERNA_DB_DIR}/silva-bac-23s-id98.fasta" \
  --ref "${SORTMERNA_DB_DIR}/silva-euk-18s-id95.fasta" \
  --ref "${SORTMERNA_DB_DIR}/silva-euk-28s-id98.fasta" \
  --ref "${SORTMERNA_DB_DIR}/rfam-5s-database-id98.fasta" \
  --ref "${SORTMERNA_DB_DIR}/rfam-5.8s-database-id98.fasta" \
  --reads "$R1" \
  --reads "$R2" \
  --fastx \
  --paired_out \
  --out2 \
  --aligned "${WORKDIR}/${SAMPLE}.rRNA" \
  --other   "${SORTMERNA_DIR}/${SAMPLE}" \
  --threads "$SLURM_CPUS_PER_TASK" \
  --workdir "$WORKDIR"

rm -rf "$WORKDIR"

chmod u+rw,g+rw "${SORTMERNA_DIR}/${SAMPLE}_fwd.fq.gz" \
                "${SORTMERNA_DIR}/${SAMPLE}_rev.fq.gz"

# Bidirectional sync: keep only read names present in BOTH fwd and rev
module load stack/2024-06
module load "${SEQTK_MODULE}"
SYNCDIR="${TMP_ROOT}/sync_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$SYNCDIR"

FWD="${SORTMERNA_DIR}/${SAMPLE}_fwd.fq.gz"
REV="${SORTMERNA_DIR}/${SAMPLE}_rev.fq.gz"

zcat "$FWD" | awk 'NR%4==1{print substr($1,2)}' > "${SYNCDIR}/fwd_names.txt"
zcat "$REV" | awk 'NR%4==1{print substr($1,2)}' > "${SYNCDIR}/rev_names.txt"

N_FWD_RAW=$(wc -l < "${SYNCDIR}/fwd_names.txt")
N_REV_RAW=$(wc -l < "${SYNCDIR}/rev_names.txt")
echo "[INFO] Pre-sync — fwd: $N_FWD_RAW  rev: $N_REV_RAW"

if [[ "$N_FWD_RAW" -eq "$N_REV_RAW" ]]; then
  echo "[INFO] Already equal — skipping sync"
else
  # Intersection: names present in both files
  awk 'NR==FNR{a[$0]; next} $0 in a' \
    "${SYNCDIR}/fwd_names.txt" "${SYNCDIR}/rev_names.txt" > "${SYNCDIR}/common.txt"
  echo "[INFO] Common reads: $(wc -l < "${SYNCDIR}/common.txt")"

  seqtk subseq "$FWD" "${SYNCDIR}/common.txt" | gzip -c > "${FWD}.tmp"
  seqtk subseq "$REV" "${SYNCDIR}/common.txt" | gzip -c > "${REV}.tmp"
  mv "${FWD}.tmp" "$FWD"
  mv "${REV}.tmp" "$REV"
  chmod u+rw,g+rw "$FWD" "$REV"
fi

rm -rf "$SYNCDIR"

echo "[INFO] Done: ${SORTMERNA_DIR}/${SAMPLE}_fwd.fq.gz"
