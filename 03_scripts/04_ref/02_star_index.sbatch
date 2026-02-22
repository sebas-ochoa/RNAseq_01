#!/usr/bin/env bash
#SBATCH --job-name=star_index
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4G
#SBATCH --time=06:00:00
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
module load "$STAR_MODULE"

ensure_writable_dir "$SLURM_LOG_DIR"
ensure_writable_dir "$STAR_INDEX_DIR"

if [[ ! -s "$FASTA" || ! -s "$GTF" ]]; then
  echo "[ERROR] Missing reference FASTA/GTF. Run 03_scripts/04_ref/01_download_reference.sbatch first." >&2
  exit 1
fi

FIRST_SAMPLE=$(awk 'NR==2{print $1}' "$SAMPLES_TSV")
if [[ -z "$FIRST_SAMPLE" ]]; then
  echo "[ERROR] Could not read first sample_id from $SAMPLES_TSV" >&2
  exit 1
fi

FIRST_R1="${TRIM_DIR}/${FIRST_SAMPLE}_1.trim.fq.gz"
if [[ ! -s "$FIRST_R1" ]]; then
  FIRST_R1=$(find "$TRIM_DIR" -maxdepth 1 -type f -name '*_1.trim.fq.gz' | sort | head -n 1 || true)
fi

if [[ -z "$FIRST_R1" || ! -s "$FIRST_R1" ]]; then
  echo "[ERROR] No trimmed R1 FASTQ found to infer read length in $TRIM_DIR" >&2
  exit 1
fi

READ_LEN=$(gzip -dc "$FIRST_R1" | awk 'NR==2{print length($0); exit}')
if [[ -z "$READ_LEN" || "$READ_LEN" -le 1 ]]; then
  echo "[ERROR] Could not infer read length from: $FIRST_R1" >&2
  exit 1
fi

SJDB_OVERHANG=$((READ_LEN - 1))

echo "[INFO] Building STAR index"
echo "[INFO] Read length: $READ_LEN"
echo "[INFO] sjdbOverhang: $SJDB_OVERHANG"

STAR \
  --runThreadN "$SLURM_CPUS_PER_TASK" \
  --runMode genomeGenerate \
  --genomeDir "$STAR_INDEX_DIR" \
  --genomeFastaFiles "$FASTA" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang "$SJDB_OVERHANG" \
  --outFileNamePrefix "${STAR_INDEX_DIR}/log_"

printf "metric\tvalue\nread_length\t%s\nsjdbOverhang\t%s\n" "$READ_LEN" "$SJDB_OVERHANG" > "${STAR_INDEX_DIR}/index_build_info.tsv"

echo "[INFO] STAR index ready: $STAR_INDEX_DIR"
