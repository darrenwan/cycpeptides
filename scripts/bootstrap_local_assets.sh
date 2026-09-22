#!/usr/bin/env bash
# Copy campaign-owned model weights and ProteinMPNN into this repo.
# Safe to re-run (skips existing files unless FORCE=1).
#
# Usage:
#   bash scripts/bootstrap_local_assets.sh
#   FORCE=1 bash scripts/bootstrap_local_assets.sh
#
# Optional source overrides (defaults are historical sibling paths):
#   SRC_RFD_MODELS   SRC_AFCYC_PARAMS   SRC_PROTEINMPNN
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FORCE="${FORCE:-0}"

SRC_RFD_MODELS="${SRC_RFD_MODELS:-/mnt/data4t/aidd/models/rfmodels}"
SRC_AFCYC_PARAMS="${SRC_AFCYC_PARAMS:-/mnt/data4t/cursor_project/RFdiffusion/projects/pdl1_macrocycle/models/afcyc}"
SRC_PROTEINMPNN="${SRC_PROTEINMPNN:-/mnt/data4t/kg_project/ProteinMPNN}"

need() {
  local f="$1"
  if [[ "${FORCE}" == "1" ]]; then return 0; fi
  [[ ! -e "$f" ]]
}

mkdir -p "${ROOT}/models/rfdiffusion" "${ROOT}/models/afcyc" "${ROOT}/third_party"

# --- RFdiffusion Complex_base ---
DST_CKPT="${ROOT}/models/rfdiffusion/Complex_base_ckpt.pt"
if need "${DST_CKPT}"; then
  if [[ ! -f "${SRC_RFD_MODELS}/Complex_base_ckpt.pt" ]]; then
    echo "Missing source: ${SRC_RFD_MODELS}/Complex_base_ckpt.pt" >&2
    exit 1
  fi
  echo "→ models/rfdiffusion/Complex_base_ckpt.pt"
  cp -a "${SRC_RFD_MODELS}/Complex_base_ckpt.pt" "${DST_CKPT}"
  [[ -f "${SRC_RFD_MODELS}/url.txt" ]] && cp -a "${SRC_RFD_MODELS}/url.txt" "${ROOT}/models/rfdiffusion/" || true
else
  echo "skip ${DST_CKPT}"
fi

# --- AfCyc AF2 params (params/ only) ---
DST_PARAMS="${ROOT}/models/afcyc/params"
if need "${DST_PARAMS}"; then
  if [[ ! -d "${SRC_AFCYC_PARAMS}/params" ]]; then
    echo "Missing source: ${SRC_AFCYC_PARAMS}/params" >&2
    exit 1
  fi
  echo "→ models/afcyc/params/  (large ~5GB)"
  rsync -a --info=progress2 "${SRC_AFCYC_PARAMS}/params" "${ROOT}/models/afcyc/"
else
  echo "skip ${DST_PARAMS}"
fi

# --- ProteinMPNN (minimal vendor) ---
DST_MPNN="${ROOT}/third_party/ProteinMPNN"
if need "${DST_MPNN}/protein_mpnn_run.py"; then
  if [[ ! -f "${SRC_PROTEINMPNN}/protein_mpnn_run.py" ]]; then
    echo "Missing source: ${SRC_PROTEINMPNN}/protein_mpnn_run.py" >&2
    exit 1
  fi
  echo "→ third_party/ProteinMPNN/"
  mkdir -p "${DST_MPNN}/helper_scripts" "${DST_MPNN}/vanilla_model_weights"
  cp -a "${SRC_PROTEINMPNN}/protein_mpnn_run.py" "${DST_MPNN}/"
  cp -a "${SRC_PROTEINMPNN}/protein_mpnn_utils.py" "${DST_MPNN}/"
  cp -a "${SRC_PROTEINMPNN}/LICENSE" "${DST_MPNN}/" 2>/dev/null || true
  cp -a "${SRC_PROTEINMPNN}/README.md" "${DST_MPNN}/" 2>/dev/null || true
  cp -a "${SRC_PROTEINMPNN}/helper_scripts/parse_multiple_chains.py" "${DST_MPNN}/helper_scripts/"
  cp -a "${SRC_PROTEINMPNN}/helper_scripts/assign_fixed_chains.py" "${DST_MPNN}/helper_scripts/"
  cp -a "${SRC_PROTEINMPNN}/vanilla_model_weights/"*.pt "${DST_MPNN}/vanilla_model_weights/"
else
  echo "skip ${DST_MPNN}"
fi

# --- ColabDesign (AfCyc Docker build context) ---
bash "${ROOT}/scripts/fetch_colabdesign.sh"

echo
echo "OK. Layout:"
du -sh "${ROOT}/models/rfdiffusion" "${ROOT}/models/afcyc" "${ROOT}/third_party/ProteinMPNN" "${ROOT}/third_party/ColabDesign" 2>/dev/null || true
echo "Provenance: ${ROOT}/models/SOURCES.md"
