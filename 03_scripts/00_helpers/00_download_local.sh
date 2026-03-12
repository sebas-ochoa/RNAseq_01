#!/usr/bin/env bash
# Run locally (Mac) to download SortMeRNA and FastQ Screen references.
# After completion, copy staging/ to the cluster with the scp command shown.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
STAGING="${PROJECT_ROOT}/staging"

CLUSTER_USER="jochoa"
CLUSTER_HOST="euler.ethz.ch"
CLUSTER_PROJECT="/nfs/nas22/fs2202/biol_bc_jagannathan_comp_1/Sebastian/013.0_RNAseq_hsf4_14d29_0-6h25"

echo "[INFO] Staging dir: $STAGING"
mkdir -p \
  "${STAGING}/sortmerna/db" \
  "${STAGING}/fastq_screen" \
  "${STAGING}/fastq_screen/refs/drosophila" \
  "${STAGING}/fastq_screen/refs/human" \
  "${STAGING}/fastq_screen/refs/ecoli" \
  "${STAGING}/fastq_screen/refs/phix" \
  "${STAGING}/fastq_screen/refs/vectors"

# ---------------------------------------------------------------------------
# SortMeRNA binary (Linux x86_64)
# ---------------------------------------------------------------------------
echo "[INFO] Downloading SortMeRNA..."
SMRNA_VER="4.3.7"
curl -fL \
  "https://github.com/biocore/sortmerna/releases/download/v${SMRNA_VER}/sortmerna-${SMRNA_VER}-Linux.sh" \
  -o "${STAGING}/sortmerna/sortmerna-${SMRNA_VER}-Linux.sh"

# ---------------------------------------------------------------------------
# SortMeRNA rRNA databases
# ---------------------------------------------------------------------------
echo "[INFO] Downloading rRNA databases..."
DB_BASE="https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases"
for DB in \
  silva-arc-16s-id95.fasta \
  silva-arc-23s-id98.fasta \
  silva-bac-16s-id90.fasta \
  silva-bac-23s-id98.fasta \
  silva-euk-18s-id95.fasta \
  silva-euk-28s-id98.fasta \
  rfam-5s-database-id98.fasta \
  rfam-5.8s-database-id98.fasta
do
  echo "  -> $DB"
  curl -fL "${DB_BASE}/${DB}" -o "${STAGING}/sortmerna/db/${DB}"
done

# ---------------------------------------------------------------------------
# FastQ Screen (Perl script, GitHub source archive)
# ---------------------------------------------------------------------------
echo "[INFO] Downloading FastQ Screen..."
FQS_VER="0.16.0"
curl -fL \
  "https://github.com/StevenWingett/FastQ-Screen/archive/refs/tags/v${FQS_VER}.tar.gz" \
  -o "${STAGING}/fastq_screen/fastq_screen_v${FQS_VER}.tar.gz"

# ---------------------------------------------------------------------------
# FastQ Screen reference genomes
# ---------------------------------------------------------------------------

# Drosophila melanogaster BDGP6.46 (~47MB) — needed for bowtie2 index
echo "[INFO] Downloading Drosophila melanogaster FASTA..."
curl -fL \
  "https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz" \
  -o "${STAGING}/fastq_screen/refs/drosophila/drosophila.fa.gz"

# Human GRCh38 primary assembly (no alt, ~900MB)
echo "[INFO] Downloading Human GRCh38 (primary assembly)..."
curl -fL \
  "https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
  -o "${STAGING}/fastq_screen/refs/human/GRCh38.primary.fa.gz"

# E. coli K-12 MG1655 (~5MB)
echo "[INFO] Downloading E. coli K-12 MG1655..."
curl -fL \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz" \
  -o "${STAGING}/fastq_screen/refs/ecoli/ecoli_k12_mg1655.fa.gz"

# PhiX174 (~5KB)
echo "[INFO] Downloading PhiX174..."
curl -fL \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz" \
  -o "${STAGING}/fastq_screen/refs/phix/phix174.fa.gz"

# UniVec_Core — common cloning vectors (~1MB)
echo "[INFO] Downloading UniVec_Core..."
curl -fL \
  "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core" \
  -o "${STAGING}/fastq_screen/refs/vectors/UniVec_Core.fa"

# Illumina adapters (inlined — TruSeq + Nextera)
cat > "${STAGING}/fastq_screen/refs/adapters.fa" <<'ADAPTERS'
>TruSeq_Universal
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>TruSeq_R1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>TruSeq_R2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>Nextera_R1
CTGTCTCTTATACACATCT
>Nextera_R2
CTGTCTCTTATACACATCT
>PolyA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>PolyC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
ADAPTERS

echo ""
echo "==================================================================="
echo " Download complete. Transfer to cluster with:"
echo "==================================================================="
echo ""
echo "  scp -r '${STAGING}/' \\"
echo "    ${CLUSTER_USER}@${CLUSTER_HOST}:${CLUSTER_PROJECT}/00_admin/env/staging/"
echo ""
echo " Then on the cluster run:"
echo "  sbatch ${CLUSTER_PROJECT}/03_scripts/00_helpers/01_setup_cluster.sh"
echo "==================================================================="
