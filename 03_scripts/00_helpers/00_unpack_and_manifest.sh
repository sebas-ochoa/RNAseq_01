#!/usr/bin/env bash
set -euo pipefail

# Run from project root
RAW_DIR="01_raw"
OUT_FASTQ="${RAW_DIR}/fastq"
META_DIR="00_admin/metadata"

TAR=$(ls -1 ${RAW_DIR}/*.tar 2>/dev/null | head -n 1 || true)
if [[ -z "${TAR}" ]]; then
  echo "[ERROR] No .tar found in ${RAW_DIR}/"
  exit 1
fi

mkdir -p "${OUT_FASTQ}" "${META_DIR}"

echo "[INFO] TAR: ${TAR}"
echo "[INFO] Listing FASTQ entries inside TAR..."
tar -tf "${TAR}" | grep -E '\.fq\.gz$' | sort > "${RAW_DIR}/tar_fastq_paths.txt"

N_IN_TAR=$(wc -l < "${RAW_DIR}/tar_fastq_paths.txt" | tr -d ' ')
echo "[INFO] FASTQ(.fq.gz) in TAR: ${N_IN_TAR}"
if [[ "${N_IN_TAR}" -eq 0 ]]; then
  echo "[ERROR] No .fq.gz found inside TAR."
  exit 1
fi

echo "[INFO] Extracting ONLY .fq.gz into ${OUT_FASTQ}/ (flattened)..."
# -i : ignore zeros
# --wildcards : allow pattern matching
# We extract to OUT_FASTQ, then flatten by moving basename into OUT_FASTQ
tar -xf "${TAR}" -C "${OUT_FASTQ}" --wildcards --no-anchored "*.fq.gz"

# Flatten: move any nested .fq.gz up to OUT_FASTQ/
echo "[INFO] Flattening extracted FASTQs..."
find "${OUT_FASTQ}" -type f -name "*.fq.gz" -print0 | while IFS= read -r -d '' f; do
  bn=$(basename "$f")
  if [[ "$f" != "${OUT_FASTQ}/${bn}" ]]; then
    mv -f "$f" "${OUT_FASTQ}/${bn}"
  fi
done

# Clean now-empty dirs (optional, keeps it tidy)
find "${OUT_FASTQ}" -type d -empty -delete || true

echo "[INFO] Writing manifest..."
find "${OUT_FASTQ}" -maxdepth 1 -type f -name "*.fq.gz" | sort > "${META_DIR}/fastq_manifest.txt"
N=$(wc -l < "${META_DIR}/fastq_manifest.txt" | tr -d ' ')
echo "[INFO] FASTQ files on disk: ${N}"
if [[ "${N}" -eq 0 ]]; then
  echo "[ERROR] Extraction produced 0 FASTQ files."
  exit 1
fi

echo "[INFO] Building samples.tsv (pairs _1/_2)..."
echo -e "sample_id\tR1\tR2" > "${META_DIR}/samples.tsv"

# Pair by "_1.fq.gz" -> "_2.fq.gz"
# sample_id from filename prefix BEFORE first "_" (e.g., HSF4_2H_1_..._1.fq.gz -> HSF4_2H_1)
find "${OUT_FASTQ}" -maxdepth 1 -type f -name "*_1.fq.gz" | sort | while read -r r1; do
  r2="${r1/_1.fq.gz/_2.fq.gz}"
  if [[ ! -f "$r2" ]]; then
    echo "[WARN] Missing pair for: $(basename "$r1")" >&2
    continue
  fi

  bn=$(basename "$r1")
  # sample_id is the first 3 underscore-separated fields: HSF4_0H_1
  sample_id=$(echo "$bn" | cut -d'_' -f1-3)

  echo -e "${sample_id}\t${r1}\t${r2}" >> "${META_DIR}/samples.tsv"
done

echo "[INFO] samples.tsv lines:"
wc -l "${META_DIR}/samples.tsv" | awk '{print "[INFO] " $0}'

echo "[INFO] Preview:"
column -t -s $'\t' "${META_DIR}/samples.tsv" | head -n 20

echo "[INFO] Done."
