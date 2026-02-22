#!/usr/bin/env bash
#SBATCH --job-name=ref_download
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00
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

download_file() {
  local url="$1"
  local out="$2"

  if [[ -s "$out" ]]; then
    echo "[INFO] Reusing existing file: $out"
    return 0
  fi

  echo "[INFO] Downloading: $url"
  if command -v curl >/dev/null 2>&1; then
    curl --fail --location --retry 3 --retry-delay 5 --output "$out" "$url"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "$out" "$url"
  else
    echo "[ERROR] Neither curl nor wget is available on this node." >&2
    exit 1
  fi
}

PROJECT_ROOT="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is empty}"
cd "$PROJECT_ROOT"
source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

module purge
module load "$STACK_MODULE"

ensure_writable_dir "$SLURM_LOG_DIR"
ensure_writable_dir "$REFERENCE_DIR"

download_file "$FASTA_URL" "$FASTA_GZ"
download_file "$GTF_URL" "$GTF_GZ"

if [[ ! -s "$FASTA_GZ" || ! -s "$GTF_GZ" ]]; then
  echo "[ERROR] Downloaded reference files are missing or empty." >&2
  exit 1
fi

if [[ ! -s "$FASTA" ]]; then
  echo "[INFO] Decompressing FASTA..."
  gzip -dc "$FASTA_GZ" > "$FASTA"
fi

if [[ ! -s "$GTF" ]]; then
  echo "[INFO] Decompressing GTF..."
  gzip -dc "$GTF_GZ" > "$GTF"
fi

if ! head -n 1 "$FASTA" | grep -q '^>'; then
  echo "[ERROR] FASTA validation failed (no header line starting with '>')." >&2
  exit 1
fi

if ! awk 'BEGIN{ok=0} $0 !~ /^#/ && NF>=9 {ok=1; exit} END{exit ok?0:1}' "$GTF"; then
  echo "[ERROR] GTF validation failed (no feature rows with 9+ fields)." >&2
  exit 1
fi

echo "[INFO] Reference ready:"
echo "[INFO] FASTA: $FASTA"
echo "[INFO] GTF:   $GTF"
