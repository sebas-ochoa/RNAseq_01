#!/usr/bin/env bash
#SBATCH --job-name=setup_pipeline
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --output=04_logs/slurm/%x_%j.out
#SBATCH --error=04_logs/slurm/%x_%j.err
#
# Install SortMeRNA and FastQ Screen, then build Bowtie2 indices.
# Run from project root: sbatch 03_scripts/00_helpers/01_setup_cluster.sh
set -euo pipefail
umask 007

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"
source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

STAGING="${PROJECT_ROOT}/00_admin/env/staging"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "$SORTMERNA_HOME" "$SORTMERNA_DB_DIR" "$FASTQ_SCREEN_HOME" "04_logs/slurm"
chmod -R u+rwX,g+rwX "$SORTMERNA_HOME" "$FASTQ_SCREEN_HOME" 2>/dev/null || true

# ---------------------------------------------------------------------------
# SortMeRNA
# ---------------------------------------------------------------------------
echo "[INFO] Installing SortMeRNA..."
SMRNA_VER="4.3.7"
SMRNA_INSTALLER="${STAGING}/sortmerna/sortmerna-${SMRNA_VER}-Linux.sh"
bash "$SMRNA_INSTALLER" --skip-license --prefix="$SORTMERNA_HOME"

cp "${STAGING}/sortmerna/db/"*.fasta "$SORTMERNA_DB_DIR/"
echo "[INFO] SortMeRNA installed: ${SORTMERNA_HOME}/bin/sortmerna"

# ---------------------------------------------------------------------------
# FastQ Screen
# ---------------------------------------------------------------------------
echo "[INFO] Installing FastQ Screen..."
FQS_VER="0.16.0"
tar -xzf "${STAGING}/fastq_screen/fastq_screen_v${FQS_VER}.tar.gz" \
  -C "$FASTQ_SCREEN_HOME" --strip-components=1
chmod +x "${FASTQ_SCREEN_HOME}/fastq_screen"
echo "[INFO] FastQ Screen installed: ${FASTQ_SCREEN_HOME}/fastq_screen"

# ---------------------------------------------------------------------------
# Build Bowtie2 indices for FastQ Screen
# ---------------------------------------------------------------------------
module purge
module load "$STACK_MODULE"
module load "$BOWTIE2_MODULE"

IDX_DIR="${FASTQ_SCREEN_HOME}/indices"
mkdir -p \
  "${IDX_DIR}/drosophila" \
  "${IDX_DIR}/human" \
  "${IDX_DIR}/ecoli" \
  "${IDX_DIR}/phix" \
  "${IDX_DIR}/vectors" \
  "${IDX_DIR}/adapters"

echo "[INFO] Building Bowtie2 index: Drosophila..."
if [[ -r "$FASTA" ]]; then
  DROSO_REF="$FASTA"
elif [[ -r "$FASTA_GZ" ]]; then
  DROSO_REF="$FASTA_GZ"
else
  echo "[ERROR] Drosophila FASTA not found or not readable at $FASTA or $FASTA_GZ" >&2
  exit 1
fi
echo "[INFO] Using: $DROSO_REF"
bowtie2-build --threads "$THREADS" \
  "$DROSO_REF" "${IDX_DIR}/drosophila/drosophila"

echo "[INFO] Building Bowtie2 index: Human GRCh38..."
bowtie2-build --threads "$THREADS" \
  "${STAGING}/fastq_screen/refs/human/GRCh38.primary.fa.gz" \
  "${IDX_DIR}/human/human"

echo "[INFO] Building Bowtie2 index: E. coli..."
bowtie2-build --threads "$THREADS" \
  "${STAGING}/fastq_screen/refs/ecoli/ecoli_k12_mg1655.fa.gz" \
  "${IDX_DIR}/ecoli/ecoli"

echo "[INFO] Building Bowtie2 index: PhiX..."
bowtie2-build --threads "$THREADS" \
  "${STAGING}/fastq_screen/refs/phix/phix174.fa.gz" \
  "${IDX_DIR}/phix/phix"

echo "[INFO] Building Bowtie2 index: Adapters..."
bowtie2-build --threads "$THREADS" \
  "${STAGING}/fastq_screen/refs/adapters.fa" \
  "${IDX_DIR}/adapters/adapters"

echo "[INFO] Building Bowtie2 index: Vectors (UniVec_Core)..."
bowtie2-build --threads "$THREADS" \
  "${STAGING}/fastq_screen/refs/vectors/UniVec_Core.fa" \
  "${IDX_DIR}/vectors/vectors"

# ---------------------------------------------------------------------------
# Generate fastq_screen.conf
# ---------------------------------------------------------------------------
cat > "$FASTQ_SCREEN_CONF" <<EOF
BOWTIE2 $(which bowtie2)
THREADS $THREADS

DATABASE Drosophila    ${IDX_DIR}/drosophila/drosophila
DATABASE Human         ${IDX_DIR}/human/human
DATABASE E_coli        ${IDX_DIR}/ecoli/ecoli
DATABASE PhiX          ${IDX_DIR}/phix/phix
DATABASE Adapters      ${IDX_DIR}/adapters/adapters
DATABASE Vectors       ${IDX_DIR}/vectors/vectors
EOF

echo "[INFO] fastq_screen.conf written: $FASTQ_SCREEN_CONF"

# ---------------------------------------------------------------------------
# Cleanup staging
# ---------------------------------------------------------------------------
rm -rf "$STAGING"
echo "[INFO] Setup complete. Staging directory removed."
