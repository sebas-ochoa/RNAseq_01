#!/usr/bin/env bash
set -euo pipefail
umask 007

RAW_DIR="01_raw"
OUT_FASTQ="${RAW_DIR}/fastq"
META_DIR="00_admin/metadata"

TAR="$(ls -1 "${RAW_DIR}"/*.tar 2>/dev/null | head -n 1 || true)"
[[ -n "$TAR" ]] || { echo "[ERROR] No .tar found in ${RAW_DIR}/"; exit 1; }

mkdir -p "${OUT_FASTQ}" "${META_DIR}"
chmod -R u+rwX,g+rwX "${OUT_FASTQ}" "${META_DIR}" 2>/dev/null || true

echo "[INFO] TAR: ${TAR}"
tar -tf "${TAR}" | grep -E '\.fq\.gz$' | sort > "${RAW_DIR}/tar_fastq_paths.txt"

tar -xf "${TAR}" -C "${OUT_FASTQ}" --wildcards --no-anchored "*.fq.gz"

find "${OUT_FASTQ}" -type f -name "*.fq.gz" -print0 | while IFS= read -r -d '' f; do
  bn="$(basename "$f")"
  [[ "$f" == "${OUT_FASTQ}/${bn}" ]] || mv -f "$f" "${OUT_FASTQ}/${bn}"
done

find "${OUT_FASTQ}" -type d -empty -delete || true

find "${OUT_FASTQ}" -maxdepth 1 -type f -name "*.fq.gz" | sort > "${META_DIR}/fastq_manifest.txt"

{
  echo -e "sample_id\tR1\tR2"
  find "${OUT_FASTQ}" -maxdepth 1 -type f -name "*_1.fq.gz" | sort | while read -r r1; do
    r2="${r1/_1.fq.gz/_2.fq.gz}"
    [[ -f "$r2" ]] || continue
    sample_id="$(basename "$r1" | cut -d'_' -f1-3)"
    echo -e "${sample_id}\t${r1}\t${r2}"
  done
} > "${META_DIR}/samples.tsv"

echo "[INFO] FASTQ files: $(wc -l < "${META_DIR}/fastq_manifest.txt" | tr -d ' ')"
echo "[INFO] Samples: $(( $(wc -l < "${META_DIR}/samples.tsv" | tr -d ' ') - 1 ))"
column -t -s $'\t' "${META_DIR}/samples.tsv" | head -n 20
