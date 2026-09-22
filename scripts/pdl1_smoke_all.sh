#!/usr/bin/env bash
# S2 end-to-end smoke orchestrator (three-stage required).
# - Campaign Python: always `uv run`
# - GPU: Docker RFpeptides → ProteinMPNN → AfCycDesign
#
# Usage:
#   bash scripts/pdl1_smoke_all.sh
#   NUM_DESIGNS=4 SEQS_PER_BB=2 AFCYC_LIMIT=2 bash scripts/pdl1_smoke_all.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${ROOT}"

export UV_LINK_MODE="${UV_LINK_MODE:-copy}"
export NUM_DESIGNS="${NUM_DESIGNS:-8}"
export SEQS_PER_BB="${SEQS_PER_BB:-2}"

echo "[0] uv sync"
uv sync

echo "[1] S1 gate check (thresholds + cropped receptor)"
uv run python scripts/pdl1_prepare_receptor.py
test -f data/targets/PDL1/thresholds_v1.yaml || uv run python scripts/pdl1_calibrate_filters.py

echo "[2] S2a RFpeptides Docker smoke N=${NUM_DESIGNS}"
bash scripts/pdl1_smoke_rfpeptides.sh

echo "[3] S2b ProteinMPNN Docker smoke"
bash scripts/pdl1_smoke_mpnn.sh

echo "[4] S2c AfCycDesign Docker smoke (required)"
bash scripts/pdl1_smoke_afcyc.sh

echo "[5] S2d Score with uv (prefers AfCyc preds)"
uv run python scripts/pdl1_smoke_score.py

echo "S2 three-stage smoke finished → results/PDL1/smoke/"
