#!/usr/bin/env bash
# Download and install Qualimap 2.3 locally under 00_admin/env/qualimap/
# Run once on the cluster from the project root:
#   bash 03_scripts/00_helpers/00_install_qualimap.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$PROJECT_ROOT"

source "${PROJECT_ROOT}/00_admin/config/pipeline.env"

QUALIMAP_HOME="${QUALIMAP_HOME:-${PROJECT_ROOT}/00_admin/env/qualimap}"
QUALIMAP_VERSION="2.3"
QUALIMAP_ZIP="qualimap_v${QUALIMAP_VERSION}.zip"
QUALIMAP_URL="https://bitbucket.org/kokonech/qualimap/downloads/${QUALIMAP_ZIP}"

if [[ -x "${QUALIMAP_HOME}/qualimap" ]]; then
  echo "[INFO] Qualimap already installed at ${QUALIMAP_HOME}/qualimap"
  "${QUALIMAP_HOME}/qualimap" --version 2>&1 || true
  exit 0
fi

mkdir -p "$QUALIMAP_HOME" "$TMP_ROOT"

TMP_ZIP="${TMP_ROOT}/${QUALIMAP_ZIP}"
echo "[INFO] Downloading qualimap v${QUALIMAP_VERSION}..."
wget -q -O "$TMP_ZIP" "$QUALIMAP_URL"

echo "[INFO] Extracting..."
EXTRACT_DIR="${TMP_ROOT}/qualimap_extract_$$"
mkdir -p "$EXTRACT_DIR"
unzip -q "$TMP_ZIP" -d "$EXTRACT_DIR"

# The zip contains a single top-level directory (qualimap_v2.3/)
INNER_DIR="$(find "$EXTRACT_DIR" -maxdepth 1 -mindepth 1 -type d | head -1)"
cp -r "${INNER_DIR}/." "$QUALIMAP_HOME/"

chmod +x "${QUALIMAP_HOME}/qualimap"

# Cleanup
rm -rf "$TMP_ZIP" "$EXTRACT_DIR"

echo "[INFO] Testing installation..."
"${QUALIMAP_HOME}/qualimap" --version 2>&1 || true
echo "[INFO] Qualimap installed at: ${QUALIMAP_HOME}"
echo "[INFO] Next step: verify JAVA_MODULE in pipeline.env (run: module avail java)"
