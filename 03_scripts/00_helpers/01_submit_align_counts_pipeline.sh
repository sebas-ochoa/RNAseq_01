#!/usr/bin/env bash
set -euo pipefail
umask 007

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$PROJECT_ROOT"

source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

ensure_writable_dir() {
  local dir="$1"
  mkdir -p "$dir"
  chmod 2770 "$dir" 2>/dev/null || true

  if [[ ! -d "$dir" || ! -w "$dir" ]]; then
    echo "[ERROR] Directory missing or not writable: $dir" >&2
    ls -ld "$dir" >&2 || true
    return 1
  fi

  local probe="${dir}/.write_test_submit_$$"
  if ! touch "$probe" 2>/dev/null; then
    echo "[ERROR] Could not create file in: $dir" >&2
    ls -ld "$dir" >&2 || true
    return 1
  fi
  rm -f "$probe"
}

if [[ ! -s "$SAMPLES_TSV" ]]; then
  echo "[ERROR] Missing samples table: $SAMPLES_TSV" >&2
  exit 1
fi

N_SAMPLES=$(awk 'NR>1 && NF>=3 {n++} END{print n+0}' "$SAMPLES_TSV")
if [[ "$N_SAMPLES" -le 0 ]]; then
  echo "[ERROR] No samples found in $SAMPLES_TSV" >&2
  exit 1
fi

echo "[INFO] Running preflight checks..."
ensure_writable_dir "$SLURM_LOG_DIR"
ensure_writable_dir "$TMP_ROOT"
ensure_writable_dir "$ALIGN_DIR"
ensure_writable_dir "$ALIGN_QC_DIR"
ensure_writable_dir "$STRAND_SWEEP_DIR"
ensure_writable_dir "$COUNTS_DIR"
ensure_writable_dir "${PROJECT_ROOT}/00_admin/metadata"

if [[ ! -d "$STAR_INDEX_DIR" || ! -s "${STAR_INDEX_DIR}/SA" ]]; then
  echo "[ERROR] STAR index not found/complete in: $STAR_INDEX_DIR" >&2
  echo "[ERROR] Build index first (e.g., 03_scripts/04_ref/02_star_index.sbatch)." >&2
  exit 1
fi

if [[ ! -s "$GTF" ]]; then
  echo "[ERROR] GTF not found: $GTF" >&2
  echo "[ERROR] Set GTF path in 00_admin/config/pipeline.env or place it under REFERENCE_DIR." >&2
  exit 1
fi

if ! compgen -G "${TRIM_DIR}/*_1.trim.fq.gz" > /dev/null; then
  echo "[ERROR] No trimmed FASTQ files found in: $TRIM_DIR" >&2
  echo "[ERROR] Run trimming first." >&2
  exit 1
fi

echo "[INFO] Samples detected: $N_SAMPLES"
echo "[INFO] Starting from alignment stage (download/index skipped)."

J_ALIGN=$(sbatch --parsable --array="1-${N_SAMPLES}" "${PROJECT_ROOT}/03_scripts/05_align/01_star_align_trim_array.sbatch")
J_ALIGN_QC=$(sbatch --parsable --dependency="afterok:${J_ALIGN}" "${PROJECT_ROOT}/03_scripts/05_align/02_align_qc_multiqc.sbatch")
J_SWEEP=$(sbatch --parsable --dependency="afterok:${J_ALIGN}" "${PROJECT_ROOT}/03_scripts/06_counts/01_featurecounts_strand_sweep.sbatch")
J_FINAL=$(sbatch --parsable --dependency="afterok:${J_SWEEP}" "${PROJECT_ROOT}/03_scripts/06_counts/02_featurecounts_final.sbatch")

echo "[INFO] Submitted pipeline with dependencies:"
printf "  %-24s %s\n" "star_align_array" "$J_ALIGN"
printf "  %-24s %s\n" "align_multiqc" "$J_ALIGN_QC"
printf "  %-24s %s\n" "strand_sweep" "$J_SWEEP"
printf "  %-24s %s\n" "final_featurecounts" "$J_FINAL"

echo "[INFO] Track jobs with:"
echo "  squeue -u $USER"
echo "  sacct -j ${J_FINAL} --format=JobID,JobName%30,State,ExitCode,Elapsed,MaxRSS,NCPUS"
