#!/usr/bin/env bash
#SBATCH --job-name=star_index
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4G
#SBATCH --time=06:00:00
#SBATCH --output=04_logs/slurm/%x_%j.out
#SBATCH --error=04_logs/slurm/%x_%j.err

set -euo pipefail
umask 007

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"
source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

module purge
module load "$STACK_MODULE"
module load "$GCC_MODULE"
module load "$STAR_MODULE"

mkdir -p "$SLURM_LOG_DIR" "$STAR_INDEX_DIR"
chmod -R u+rwX,g+rwX "$SLURM_LOG_DIR" "$STAR_INDEX_DIR" 2>/dev/null || true

mapfile -t R1_FASTQ < <(find "$TRIM_DIR" -maxdepth 1 -type f -name '*_1.trim.fq.gz' | sort)
FIRST_R1="${R1_FASTQ[0]:-}"
[[ -n "$FIRST_R1" ]] || { echo "[ERROR] No trimmed R1 found in $TRIM_DIR"; exit 1; }

READ_LEN="$(gzip -dc "$FIRST_R1" | awk 'NR==2{print length($0); exit}')"
SJDB_OVERHANG=$((READ_LEN - 1))

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
