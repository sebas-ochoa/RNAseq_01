#!/usr/bin/env bash
#SBATCH --job-name=fc_final
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --output=04_logs/slurm/%x_%j.out
#SBATCH --error=04_logs/slurm/%x_%j.err

set -euo pipefail
umask 007

on_error() {
  local line="$1"
  local code="$2"
  echo "[ERROR] Job ${SLURM_JOB_ID:-NA} failed at line ${line}: ${BASH_COMMAND}" >&2
  exit "$code"
}
trap 'on_error $LINENO $?' ERR

ensure_writable_dir() {
  local dir="$1"
  mkdir -p "$dir"
  chmod 2770 "$dir" 2>/dev/null || true
  if [[ ! -d "$dir" || ! -w "$dir" ]]; then
    echo "[ERROR] Directory missing or not writable: $dir" >&2
    ls -ld "$dir" >&2 || true
    exit 1
  fi
  local probe="${dir}/.write_test_${SLURM_JOB_ID:-manual}"
  if ! touch "$probe" 2>/dev/null; then
    echo "[ERROR] Could not create file in: $dir" >&2
    ls -ld "$dir" >&2 || true
    exit 1
  fi
  rm -f "$probe"
}

PROJECT_ROOT="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is empty}"
cd "$PROJECT_ROOT"
source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

module purge
module load "$STACK_MODULE"
module load "$GCC_MODULE"
module load "$SUBREAD_MODULE"

ensure_writable_dir "$SLURM_LOG_DIR"
ensure_writable_dir "$COUNTS_DIR"

DECISION_TSV="${STRAND_SWEEP_DIR}/strand_decision.tsv"
if [[ ! -s "$DECISION_TSV" ]]; then
  echo "[ERROR] Missing strand decision file: $DECISION_TSV" >&2
  echo "[ERROR] Run 03_scripts/06_counts/01_featurecounts_strand_sweep.sbatch first." >&2
  exit 1
fi

BEST_S=$(awk 'NR==2{print $1}' "$DECISION_TSV")
if [[ ! "$BEST_S" =~ ^[012]$ ]]; then
  echo "[ERROR] Invalid selected_s in $DECISION_TSV: $BEST_S" >&2
  exit 1
fi

if [[ ! -s "$GTF" ]]; then
  echo "[ERROR] GTF not found: $GTF" >&2
  exit 1
fi

mapfile -t BAMS < <(find "$ALIGN_DIR" -maxdepth 1 -type f -name '*.sorted.bam' | sort)
if [[ "${#BAMS[@]}" -eq 0 ]]; then
  echo "[ERROR] No BAM files found in $ALIGN_DIR" >&2
  exit 1
fi

OUT_COUNTS="${COUNTS_DIR}/gene_counts.tsv"

echo "[INFO] Running final featureCounts with -s ${BEST_S}"
featureCounts \
  -T "$SLURM_CPUS_PER_TASK" \
  -p -B -C \
  -t exon -g gene_id \
  -s "$BEST_S" \
  -a "$GTF" \
  -o "$OUT_COUNTS" \
  "${BAMS[@]}"

if [[ ! -s "$OUT_COUNTS" || ! -s "${OUT_COUNTS}.summary" ]]; then
  echo "[ERROR] Final featureCounts output missing." >&2
  exit 1
fi

MATRIX_OUT="${COUNTS_DIR}/gene_counts.matrix.tsv"
awk 'BEGIN{OFS="\t"}
  NR==2 {
    printf "Geneid"
    for (i=7; i<=NF; i++) {
      sample=$i
      sub(/^.*\//, "", sample)
      sub(/\.sorted\.bam$/, "", sample)
      printf OFS sample
    }
    printf "\n"
    next
  }
  NR>2 {
    printf "%s", $1
    for (i=7; i<=NF; i++) {
      printf OFS $i
    }
    printf "\n"
  }
' "$OUT_COUNTS" > "$MATRIX_OUT"

echo "[INFO] Final counts done: $OUT_COUNTS"
echo "[INFO] Matrix done: $MATRIX_OUT"
