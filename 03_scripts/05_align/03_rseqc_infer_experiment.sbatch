#!/usr/bin/env bash
#SBATCH --job-name=rseqc_infer_exp
#SBATCH --cpus-per-task=4
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

RSEQC_MODULE="${RSEQC_MODULE:-rseqc}"
RSEQC_DIR="${RSEQC_DIR:-${PROJECT_ROOT}/02_work/05_align/rseqc}"
RSEQC_BED="${RSEQC_BED:-${REFERENCE_DIR}/rseqc_reference.bed}"

module purge
module load "$STACK_MODULE"
module load "$GCC_MODULE"
module load "$PYTHON_MODULE"
module load "$SAMTOOLS_MODULE"

if ! module load "$RSEQC_MODULE" 2>/dev/null; then
  echo "[WARN] Could not load module '$RSEQC_MODULE'. Will try infer_experiment.py from PATH."
fi

ensure_writable_dir "$SLURM_LOG_DIR"
ensure_writable_dir "$RSEQC_DIR"
ensure_writable_dir "$REFERENCE_DIR"

if ! command -v infer_experiment.py >/dev/null 2>&1; then
  echo "[ERROR] infer_experiment.py not found in PATH." >&2
  echo "[ERROR] Load an RSeQC module or set PATH to include RSeQC binaries." >&2
  exit 1
fi

if [[ ! -s "$GTF" ]]; then
  echo "[ERROR] GTF not found: $GTF" >&2
  exit 1
fi

if [[ ! -s "$RSEQC_BED" ]]; then
  if ! command -v gtf2bed >/dev/null 2>&1; then
    echo "[ERROR] RSeQC BED missing: $RSEQC_BED" >&2
    echo "[ERROR] Could not auto-generate because gtf2bed is not available." >&2
    echo "[ERROR] Provide a BED12 gene model via RSEQC_BED and rerun." >&2
    exit 1
  fi

  echo "[INFO] Generating RSeQC BED from GTF: $RSEQC_BED"
  TMP_BED="${RSEQC_BED}.tmp.${SLURM_JOB_ID:-manual}"
  if [[ "$GTF" == *.gz ]]; then
    zcat "$GTF" | gtf2bed > "$TMP_BED"
  else
    gtf2bed < "$GTF" > "$TMP_BED"
  fi
  mv -f "$TMP_BED" "$RSEQC_BED"
fi

mapfile -t BAMS < <(find "$ALIGN_DIR" -maxdepth 1 -type f -name '*.sorted.bam' | sort)
if [[ "${#BAMS[@]}" -eq 0 ]]; then
  echo "[ERROR] No BAM files found in $ALIGN_DIR" >&2
  exit 1
fi

SUMMARY_TSV="${RSEQC_DIR}/infer_experiment_summary.tsv"
{
  echo -e "sample\tfraction_failed\tfraction_1pp_1mm_2pm_2mp\tfraction_1pm_1mp_2pp_2mm\tinferred_layout\tsuggested_featurecounts_s"
} > "$SUMMARY_TSV"

for BAM in "${BAMS[@]}"; do
  SAMPLE="$(basename "$BAM" .sorted.bam)"
  BAI="${BAM}.bai"
  OUT_TXT="${RSEQC_DIR}/${SAMPLE}.infer_experiment.txt"

  if [[ ! -s "$BAI" ]]; then
    samtools index -@ "$SLURM_CPUS_PER_TASK" "$BAM"
  fi

  infer_experiment.py -r "$RSEQC_BED" -i "$BAM" > "$OUT_TXT" 2>&1

  FAILED="$(awk -F': ' '/Fraction of reads failed to determine/ {gsub(/ /, "", $2); print $2; exit}' "$OUT_TXT")"
  FR="$(awk -F': ' '/1\+\+,1--,2\+-,2-\+/ {gsub(/ /, "", $2); gsub(/"/, "", $2); print $2; exit}' "$OUT_TXT")"
  RF="$(awk -F': ' '/1\+-,1-\+,2\+\+,2--/ {gsub(/ /, "", $2); gsub(/"/, "", $2); print $2; exit}' "$OUT_TXT")"

  if [[ -z "$FAILED" || -z "$FR" || -z "$RF" ]]; then
    echo "[WARN] Could not parse infer_experiment output for $SAMPLE" >&2
    printf "%s\tNA\tNA\tNA\tparse_failed\tNA\n" "$SAMPLE" >> "$SUMMARY_TSV"
    continue
  fi

  read -r LAYOUT SUGGESTED_S <<<"$(awk -v f="$FAILED" -v fr="$FR" -v rf="$RF" 'BEGIN {
    if (f+0 >= 0.50) {
      print "ambiguous_high_failed", 0
    } else if ((rf-fr) >= 0.10) {
      print "stranded_reverse_like", 2
    } else if ((fr-rf) >= 0.10) {
      print "stranded_forward_like", 1
    } else {
      print "ambiguous_s1_s2_close", 0
    }
  }')"

  printf "%s\t%.6f\t%.6f\t%.6f\t%s\t%s\n" \
    "$SAMPLE" "$FAILED" "$FR" "$RF" "$LAYOUT" "$SUGGESTED_S" >> "$SUMMARY_TSV"
done

RECO_TSV="${RSEQC_DIR}/infer_experiment_recommendation.tsv"
read -r MEAN_FAILED MEAN_FR MEAN_RF N_OK <<<"$(awk 'BEGIN{n=0;f=0;fr=0;rf=0}
  NR>1 && $2!="NA" {n++; f+=$2; fr+=$3; rf+=$4}
  END {
    if (n==0) {
      printf "NA NA NA 0"
    } else {
      printf "%.6f %.6f %.6f %d", f/n, fr/n, rf/n, n
    }
  }' "$SUMMARY_TSV")"

if [[ "$N_OK" -eq 0 ]]; then
  FINAL_LAYOUT="no_valid_samples"
  FINAL_S="NA"
else
  read -r FINAL_LAYOUT FINAL_S <<<"$(awk -v f="$MEAN_FAILED" -v fr="$MEAN_FR" -v rf="$MEAN_RF" 'BEGIN {
    if (f+0 >= 0.50) {
      print "ambiguous_high_failed", 0
    } else if ((rf-fr) >= 0.10) {
      print "stranded_reverse_like", 2
    } else if ((fr-rf) >= 0.10) {
      print "stranded_forward_like", 1
    } else {
      print "ambiguous_s1_s2_close", 0
    }
  }')"
fi

{
  echo -e "metric\tvalue"
  echo -e "n_samples_evaluated\t${N_OK}"
  echo -e "mean_fraction_failed\t${MEAN_FAILED}"
  echo -e "mean_fraction_1pp_1mm_2pm_2mp\t${MEAN_FR}"
  echo -e "mean_fraction_1pm_1mp_2pp_2mm\t${MEAN_RF}"
  echo -e "final_inferred_layout\t${FINAL_LAYOUT}"
  echo -e "recommended_featurecounts_s\t${FINAL_S}"
  echo -e "decision_rule\tthreshold_delta_0.10_and_failed_lt_0.50"
} > "$RECO_TSV"

echo "[INFO] infer_experiment summary: $SUMMARY_TSV"
echo "[INFO] infer_experiment recommendation: $RECO_TSV"
