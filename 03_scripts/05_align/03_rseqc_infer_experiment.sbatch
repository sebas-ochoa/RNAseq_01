#!/usr/bin/env bash
#SBATCH --job-name=rseqc_infer_exp
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --output=04_logs/slurm/%x_%j.out
#SBATCH --error=04_logs/slurm/%x_%j.err

set -euo pipefail
umask 007

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"
source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

RSEQC_DIR="${RSEQC_DIR:-${PROJECT_ROOT}/02_work/05_align/rseqc}"
RSEQC_BED="${RSEQC_BED:-${REFERENCE_DIR}/rseqc_reference.bed}"
RSEQC_VENV="${RSEQC_VENV:-${PROJECT_ROOT}/00_admin/env/rseqc_venv}"

module purge
module load "$STACK_MODULE"
module load "$GCC_MODULE"
module load "$PYTHON_MODULE"
module load "$SAMTOOLS_MODULE"

mkdir -p "$SLURM_LOG_DIR" "$RSEQC_DIR" "$REFERENCE_DIR"
chmod -R u+rwX,g+rwX "$SLURM_LOG_DIR" "$RSEQC_DIR" "$REFERENCE_DIR" 2>/dev/null || true

if [[ -f "${RSEQC_VENV}/bin/activate" ]]; then
  source "${RSEQC_VENV}/bin/activate"
fi

if [[ ! -s "$RSEQC_BED" ]]; then
  if [[ "$GTF" == *.gz ]]; then
    zcat "$GTF" | gtf2bed > "$RSEQC_BED"
  else
    gtf2bed < "$GTF" > "$RSEQC_BED"
  fi
fi

mapfile -t BAMS < <(find "$ALIGN_DIR" -maxdepth 1 -type f -name '*.sorted.bam' | sort)

SUMMARY_TSV="${RSEQC_DIR}/infer_experiment_summary.tsv"
{
  echo -e "sample\tfraction_failed\tfraction_1pp_1mm_2pm_2mp\tfraction_1pm_1mp_2pp_2mm\tinferred_layout\tsuggested_featurecounts_s"
} > "$SUMMARY_TSV"

for BAM in "${BAMS[@]}"; do
  SAMPLE="$(basename "$BAM" .sorted.bam)"
  OUT_TXT="${RSEQC_DIR}/${SAMPLE}.infer_experiment.txt"

  [[ -s "${BAM}.bai" ]] || samtools index -@ "$SLURM_CPUS_PER_TASK" "$BAM"

  infer_experiment.py -r "$RSEQC_BED" -i "$BAM" > "$OUT_TXT" 2>&1

  FAILED="$(awk -F': ' '/Fraction of reads failed to determine/ {gsub(/ /, "", $2); print $2; exit}' "$OUT_TXT")"
  FR="$(awk -F': ' '/1\+\+,1--,2\+-,2-\+/ {gsub(/ /, "", $2); gsub(/"/, "", $2); print $2; exit}' "$OUT_TXT")"
  RF="$(awk -F': ' '/1\+-,1-\+,2\+\+,2--/ {gsub(/ /, "", $2); gsub(/"/, "", $2); print $2; exit}' "$OUT_TXT")"

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
  NR>1 {n++; f+=$2; fr+=$3; rf+=$4}
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
