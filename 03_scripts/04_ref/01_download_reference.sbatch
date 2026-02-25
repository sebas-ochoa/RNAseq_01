#!/usr/bin/env bash
#SBATCH --job-name=ref_download
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00
#SBATCH --output=04_logs/slurm/%x_%j.out
#SBATCH --error=04_logs/slurm/%x_%j.err

set -euo pipefail
umask 007

download_file() {
  local url="$1"
  local out="$2"
  [[ -s "$out" ]] && return 0

  if command -v curl >/dev/null 2>&1; then
    curl --fail --location --retry 3 --retry-delay 5 --output "$out" "$url"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "$out" "$url"
  else
    echo "[ERROR] curl/wget not available on this node" >&2
    exit 1
  fi
}

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"
source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

module purge
module load "$STACK_MODULE"

mkdir -p "$SLURM_LOG_DIR" "$REFERENCE_DIR"
chmod -R u+rwX,g+rwX "$SLURM_LOG_DIR" "$REFERENCE_DIR" 2>/dev/null || true

download_file "$FASTA_URL" "$FASTA_GZ"
download_file "$GTF_URL" "$GTF_GZ"

[[ -s "$FASTA" ]] || gzip -dc "$FASTA_GZ" > "$FASTA"
[[ -s "$GTF" ]] || gzip -dc "$GTF_GZ" > "$GTF"

echo "[INFO] FASTA: $FASTA"
echo "[INFO] GTF:   $GTF"
