#!/usr/bin/env bash
#SBATCH --job-name=fc_strand_sweep
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
ensure_writable_dir "$STRAND_SWEEP_DIR"

if [[ ! -s "$GTF" ]]; then
  echo "[ERROR] GTF not found: $GTF" >&2
  exit 1
fi

mapfile -t BAMS < <(find "$ALIGN_DIR" -maxdepth 1 -type f -name '*.sorted.bam' | sort)
if [[ "${#BAMS[@]}" -eq 0 ]]; then
  echo "[ERROR] No BAM files found in $ALIGN_DIR" >&2
  exit 1
fi

METRICS_TSV="${STRAND_SWEEP_DIR}/strand_sweep_metrics.tsv"
{
  echo -e "strand_mode\tassigned_reads\ttotal_reads\tassigned_pct"
} > "$METRICS_TSV"

for s in 0 1 2; do
  OUT="${STRAND_SWEEP_DIR}/gene_counts.s${s}.txt"
  echo "[INFO] Running featureCounts with -s ${s}"

  featureCounts \
    -T "$SLURM_CPUS_PER_TASK" \
    -p -B -C \
    -t exon -g gene_id \
    -s "$s" \
    -a "$GTF" \
    -o "$OUT" \
    "${BAMS[@]}"

  SUMMARY="${OUT}.summary"
  if [[ ! -s "$SUMMARY" ]]; then
    echo "[ERROR] featureCounts summary missing: $SUMMARY" >&2
    exit 1
  fi

  ASSIGNED=$(awk 'NR>1 && $1=="Assigned" {for(i=2;i<=NF;i++) sum+=$i} END{print sum+0}' "$SUMMARY")
  TOTAL=$(awk 'NR>1 {for(i=2;i<=NF;i++) sum+=$i} END{print sum+0}' "$SUMMARY")
  PCT=$(awk -v a="$ASSIGNED" -v t="$TOTAL" 'BEGIN{if(t==0){printf "0.00"}else{printf "%.2f", (100*a)/t}}')

  printf "%s\t%s\t%s\t%s\n" "$s" "$ASSIGNED" "$TOTAL" "$PCT" >> "$METRICS_TSV"
done

BEST_S=$(awk 'NR==1{next} {if ($4+0 > best+0) {best=$4; mode=$1}} END{print mode}' "$METRICS_TSV")
if [[ -z "$BEST_S" ]]; then
  echo "[ERROR] Could not determine best strandedness from $METRICS_TSV" >&2
  exit 1
fi

case "$BEST_S" in
  0) LABEL="unstranded" ;;
  1) LABEL="stranded_forward" ;;
  2) LABEL="stranded_reverse" ;;
  *)
    echo "[ERROR] Invalid strandedness mode selected: $BEST_S" >&2
    exit 1
    ;;
esac

{
  echo -e "selected_s\tlabel\tbasis"
  echo -e "${BEST_S}\t${LABEL}\thighest_assigned_pct"
} > "${STRAND_SWEEP_DIR}/strand_decision.tsv"

echo "[INFO] Strand sweep done. Best mode: -s ${BEST_S} (${LABEL})"
echo "[INFO] Metrics: ${METRICS_TSV}"
