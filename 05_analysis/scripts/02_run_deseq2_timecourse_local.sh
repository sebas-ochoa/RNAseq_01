#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"

Rscript "${PROJECT_ROOT}/05_analysis/scripts/01_deseq2_timecourse.R" \
  --counts="${PROJECT_ROOT}/02_work/06_counts/final/gene_counts.matrix.tsv" \
  --samples="${PROJECT_ROOT}/00_admin/metadata/samples.tsv" \
  --outdir="${PROJECT_ROOT}/05_analysis/de_timecourse"
