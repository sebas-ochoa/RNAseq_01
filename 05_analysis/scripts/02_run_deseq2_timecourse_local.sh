#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"

if [[ ! -s "${PROJECT_ROOT}/02_work/06_counts/final/gene_counts.matrix.tsv" ]]; then
  echo "[ERROR] Counts matrix not found under PROJECT_ROOT=${PROJECT_ROOT}" >&2
  echo "[ERROR] Run this from project root (013_RNA_Seq) or set PROJECT_ROOT." >&2
  exit 1
fi

if [[ ! -s "${PROJECT_ROOT}/00_admin/metadata/samples.tsv" ]]; then
  echo "[ERROR] samples.tsv not found under PROJECT_ROOT=${PROJECT_ROOT}" >&2
  exit 1
fi

Rscript "${PROJECT_ROOT}/05_analysis/scripts/01_deseq2_timecourse.R" \
  --counts="${PROJECT_ROOT}/02_work/06_counts/final/gene_counts.matrix.tsv" \
  --samples="${PROJECT_ROOT}/00_admin/metadata/samples.tsv" \
  --outdir="${PROJECT_ROOT}/05_analysis/de_timecourse"
