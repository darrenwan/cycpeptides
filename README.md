# cycpeptides

Pure **in silico** AI cyclic peptide design (PD-L1 campaign).

## Project-local assets (independence)

Runtime defaults point **inside this repo** (no sibling-project paths):

| Asset | Path |
|-------|------|
| RFdiffusion `Complex_base_ckpt.pt` | `models/rfdiffusion/` |
| AfCyc / AF2 params | `models/afcyc/params/` |
| ProteinMPNN (code + vanilla weights) | `third_party/ProteinMPNN/` |

Provenance: `models/SOURCES.md`. Re-copy if missing:

```bash
bash scripts/bootstrap_local_assets.sh
```

Large `*.pt` / AF2 `params/` are gitignored; keep them on the machine for GPU runs.

## Python environment (uv only)

Do **not** use system/`conda` Python for project scripts. Use [uv](https://github.com/astral-sh/uv):

```bash
cd /mnt/data4t/aidd/cycpeptides
uv sync
uv run python scripts/pdl1_prepare_receptor.py
uv run python scripts/pdl1_blockade_metrics.py
uv run python scripts/pdl1_calibrate_filters.py
bash scripts/pdl1_smoke_all.sh          # Docker: RFD → MPNN → AfCyc; score via uv
# NUM_DESIGNS=100 bash scripts/pdl1_smoke_all.sh
```

GPU stacks (RFdiffusion / ProteinMPNN / **AfCycDesign**) run in Docker images; campaign filter/score scripts always use this uv `.venv`.

**S2 three-stage (required):** RFpeptides → ProteinMPNN (omit C) → **AfCycDesign** → `uv run` score.

Override paths only if needed: `RFD_MODELS_DIR`, `PROTEINMPNN_DIR`, `AFCYC_PARAMS_DIR`.

## Docs

- `docs/PDL1_大环肽_战役规格.md`
- `docs/PDL1_设计方案科学评估与工程方案.md`
- `docs/PDL1_结构审计与修复优先级.html` ← S0–S2 结构审计（含编号契约 P0 与修复顺序）
- `data/targets/PDL1/roadmap.yaml`
