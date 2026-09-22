#!/usr/bin/env bash
# S2 smoke: RFpeptides cyclic binder generation via Docker (not the uv venv).
# Campaign filter/scoring scripts still use: uv run python ...
#
# Usage:
#   bash scripts/pdl1_smoke_rfpeptides.sh
#   NUM_DESIGNS=8 bash scripts/pdl1_smoke_rfpeptides.sh
#
# Requires: docker + GPU, image rosettacommons/rfdiffusion:latest
# Models:   ${ROOT}/models/rfdiffusion/Complex_base_ckpt.pt (override RFD_MODELS_DIR)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PDB_HOST="${ROOT}/data/targets/PDL1/pdb/PDL1_IgV_5O45_A18-134.pdb"
OUT_DIR="${ROOT}/results/PDL1/smoke/rfpeptides"
LOG_DIR="${ROOT}/results/PDL1/smoke/logs"
MODELS_DIR="${RFD_MODELS_DIR:-${ROOT}/models/rfdiffusion}"
IMAGE="${RFD_IMAGE:-rosettacommons/rfdiffusion:latest}"
NUM_DESIGNS="${NUM_DESIGNS:-8}"
DIFFUSER_T="${DIFFUSER_T:-50}"
SHM="${RFD_SHM_SIZE:-8g}"
HOTSPOTS="${HOTSPOTS:-A56,A123,A115,A113,A66,A121}"
CONTIGS="${CONTIGS:-[12-16 A18-134/0]}"
PREFIX_HOST="${OUT_DIR}/pdl1_cyc"
PREFIX_CTR="/outputs/pdl1_cyc"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"
LOG="${LOG_DIR}/rfpeptides_smoke_$(date +%Y%m%d_%H%M%S).log"

if [[ ! -f "${PDB_HOST}" ]]; then
  echo "Missing cropped receptor. Run: uv run python scripts/pdl1_prepare_receptor.py" >&2
  exit 1
fi
if [[ ! -f "${MODELS_DIR}/Complex_base_ckpt.pt" ]]; then
  echo "Missing Complex_base_ckpt.pt under ${MODELS_DIR}" >&2
  echo "Run: bash scripts/bootstrap_local_assets.sh" >&2
  exit 1
fi

# Hydra list: ppi.hotspot_res=['A56','A123',...]
hotspot_arg() {
  local joined
  joined="$(printf "'%s'," ${HOTSPOTS//,/ })"
  echo "ppi.hotspot_res=[${joined%,}]"
}

{
  echo "=== S2 RFpeptides smoke ==="
  echo "image=${IMAGE}"
  echo "pdb=${PDB_HOST}"
  echo "models=${MODELS_DIR}"
  echo "N=${NUM_DESIGNS} contigs=${CONTIGS} T=${DIFFUSER_T}"
  echo "hotspots=${HOTSPOTS}"
  echo "out=${OUT_DIR}"
} | tee "${LOG}"

# Mount only the cropped PDB as /inputs/receptor.pdb
INPUTS_TMP="${OUT_DIR}/_inputs"
mkdir -p "${INPUTS_TMP}"
cp -f "${PDB_HOST}" "${INPUTS_TMP}/receptor.pdb"

set +e
docker run --rm --gpus all \
  --shm-size "${SHM}" \
  -v "${MODELS_DIR}:/models:ro" \
  -v "${INPUTS_TMP}:/inputs:ro" \
  -v "${OUT_DIR}:/outputs" \
  "${IMAGE}" \
  "inference.input_pdb=/inputs/receptor.pdb" \
  "inference.output_prefix=${PREFIX_CTR}" \
  "inference.model_directory_path=/models" \
  "inference.num_designs=${NUM_DESIGNS}" \
  "inference.design_startnum=0" \
  "inference.cyclic=True" \
  "inference.cyc_chains=a" \
  "contigmap.contigs=${CONTIGS}" \
  "diffuser.T=${DIFFUSER_T}" \
  "denoiser.noise_scale_ca=0" \
  "denoiser.noise_scale_frame=0" \
  "$(hotspot_arg)" \
  2>&1 | tee -a "${LOG}"
rc=${PIPESTATUS[0]}
set -e

n_pdb="$(find "${OUT_DIR}" -maxdepth 1 -name 'pdl1_cyc_*.pdb' | wc -l)"
echo "exit=${rc} n_pdb=${n_pdb}" | tee -a "${LOG}"
if [[ "${rc}" -ne 0 || "${n_pdb}" -lt 1 ]]; then
  echo "RFpeptides smoke FAILED" >&2
  exit 1
fi
echo "RFpeptides smoke OK → ${OUT_DIR}"
