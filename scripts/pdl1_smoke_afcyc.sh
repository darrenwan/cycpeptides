#!/usr/bin/env bash
# S2c: AfCycDesign structure check on MPNN-designed complexes (required for S2).
# Uses Docker image pdl1-afcyc; AF2 params mounted from host.
# Does NOT use system/conda python for campaign logic.
#
# Usage (after MPNN smoke):
#   bash scripts/pdl1_smoke_afcyc.sh
#   AFCYC_LIMIT=2 bash scripts/pdl1_smoke_afcyc.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DES_DIR="${ROOT}/results/PDL1/smoke/mpnn/designed_pdbs"
OUT_DIR="${ROOT}/results/PDL1/smoke/afcyc"
LOG_DIR="${ROOT}/results/PDL1/smoke/logs"
IMAGE="${AFCYC_IMAGE:-pdl1-afcyc:latest}"
PARAMS_HOST="${AFCYC_PARAMS_DIR:-${ROOT}/models/afcyc}"
PARAMS_CTR="/models/afcyc"
NUM_RECYCLES="${AFCYC_NUM_RECYCLES:-3}"
LIMIT="${AFCYC_LIMIT:-}"
BINDER_CHAIN="${DESIGN_CHAIN:-A}"
TARGET_CHAIN="${FIXED_CHAINS:-B}"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"
LOG="${LOG_DIR}/afcyc_smoke_$(date +%Y%m%d_%H%M%S).log"
SUMMARY="${OUT_DIR}/afcyc_summary.json"

shopt -s nullglob
PDBS=("${DES_DIR}"/*.pdb)
if [[ ${#PDBS[@]} -lt 1 ]]; then
  echo "No designed PDBs in ${DES_DIR}. Run bash scripts/pdl1_smoke_mpnn.sh first." >&2
  exit 1
fi
if [[ -n "${LIMIT}" ]]; then
  PDBS=("${PDBS[@]:0:${LIMIT}}")
fi
if [[ ! -d "${PARAMS_HOST}/params" ]]; then
  echo "Missing AF2 params under ${PARAMS_HOST}/params" >&2
  echo "Run: bash scripts/bootstrap_local_assets.sh" >&2
  echo "Or set AFCYC_PARAMS_DIR to a directory that contains params/." >&2
  exit 1
fi
if ! docker image inspect "${IMAGE}" >/dev/null 2>&1; then
  echo "Docker image ${IMAGE} not found. Build with:" >&2
  echo "  bash scripts/build_academic_images.sh afcyc" >&2
  exit 1
fi

# Build argv lists (container paths under /work)
PDB_ARGS=()
TAG_ARGS=()
for p in "${PDBS[@]}"; do
  tag="$(basename "${p}" .pdb)"
  # designed pdb mounted at /work/results/...
  rel="$(realpath --relative-to="${ROOT}" "${p}")"
  PDB_ARGS+=(--pdb "/work/${rel}")
  TAG_ARGS+=(--tag "${tag}")
done

{
  echo "=== S2c AfCycDesign smoke ==="
  echo "image=${IMAGE}"
  echo "params=${PARAMS_HOST} -> ${PARAMS_CTR}"
  echo "n=${#PDBS[@]} recycles=${NUM_RECYCLES} binder=${BINDER_CHAIN} target=${TARGET_CHAIN}"
} | tee "${LOG}"

set +e
docker run --rm --gpus all \
  -v "${ROOT}:/work" \
  -v "${PARAMS_HOST}:${PARAMS_CTR}:ro" \
  -e AFCYC_PARAMS_DIR="${PARAMS_CTR}" \
  -w /work \
  --entrypoint /opt/venv/bin/python \
  "${IMAGE}" \
  /work/scripts/pdl1_afcyc_predict.py \
  "${PDB_ARGS[@]}" \
  "${TAG_ARGS[@]}" \
  --params_dir "${PARAMS_CTR}" \
  --out_dir /work/results/PDL1/smoke/afcyc \
  --summary_json /work/results/PDL1/smoke/afcyc/afcyc_summary.json \
  --binder_chain "${BINDER_CHAIN}" \
  --target_chain "${TARGET_CHAIN}" \
  --num_recycles "${NUM_RECYCLES}" \
  --use_binder_template 1 \
  2>&1 | tee -a "${LOG}"
rc=${PIPESTATUS[0]}
set -e

n_pred="$(find "${OUT_DIR}" -maxdepth 1 -name '*_afcyc_pred.pdb' | wc -l)"
echo "exit=${rc} n_pred=${n_pred}" | tee -a "${LOG}"
if [[ "${rc}" -ne 0 || "${n_pred}" -lt 1 ]]; then
  echo "AfCyc smoke FAILED" >&2
  exit 1
fi
echo "AfCyc smoke OK → ${OUT_DIR}"
