#!/usr/bin/env bash
# Example RFpeptides inference for PD-L1 IgV (track A).
# Prefer: bash scripts/pdl1_smoke_rfpeptides.sh  (Docker + cropped A18-134 receptor)
# Campaign scoring/filters: uv run python scripts/...
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
echo "Use: bash ${ROOT}/scripts/pdl1_smoke_all.sh"
echo "Or:  NUM_DESIGNS=100 bash ${ROOT}/scripts/pdl1_smoke_rfpeptides.sh"
echo "Receptor: ${ROOT}/data/targets/PDL1/pdb/PDL1_IgV_5O45_A18-134.pdb"
echo "Python: uv run (see README.md) — do not use system/conda python for scripts/"
