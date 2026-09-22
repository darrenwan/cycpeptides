# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A **pure in silico** (dry-lab) AI design campaign for macrocyclic peptide binders of the
PD-L1 IgV domain (CD274). There is no wet-lab scope: the campaign terminus is a ranked
candidate pool (24–48 sequences) plus metrics tables, predicted complexes, and provenance.

Machine-readable campaign state lives in [data/targets/PDL1/roadmap.yaml](data/targets/PDL1/roadmap.yaml)
(stages S0–S7 with `status` and per-task DoD). Read it first — it is the authoritative
"where are we / what's next" record. Narrative docs are in Chinese:
[docs/PDL1_大环肽_战役规格.md](docs/PDL1_大环肽_战役规格.md) (spec) and
[docs/PDL1_设计方案科学评估与工程方案.md](docs/PDL1_设计方案科学评估与工程方案.md) (science review + engineering contract).

## Runtime rules (these are hard constraints, not preferences)

| Layer | Tool | Notes |
|-------|------|-------|
| Campaign scripts (`scripts/pdl1_*.py`) | **`uv` only** | `uv sync`, `uv run python scripts/...`. System/conda `python` is **forbidden** for these. |
| GPU stages (RFpeptides, ProteinMPNN, AfCycDesign) | **Docker** | Never install torch/JAX into the project venv. |
| Model weights / vendored code | in-repo | `models/`, `third_party/` — no sibling-project paths at runtime. |

The venv only holds the light campaign deps (biopython, numpy, pandas, pyyaml; ruff in the
`dev` group). Override asset locations only if needed: `RFD_MODELS_DIR`, `PROTEINMPNN_DIR`,
`AFCYC_PARAMS_DIR`.

## Common commands

```bash
# Environment
uv sync
bash scripts/bootstrap_local_assets.sh     # copy weights + ProteinMPNN into repo (idempotent; FORCE=1 to redo)

# S1: receptor prep → metrics → threshold calibration (all uv, all CPU)
uv run python scripts/pdl1_prepare_receptor.py
uv run python scripts/pdl1_blockade_metrics.py
uv run python scripts/pdl1_calibrate_filters.py

# S2: full three-stage smoke (Docker GPU + uv scoring)
bash scripts/pdl1_smoke_all.sh
NUM_DESIGNS=100 SEQS_PER_BB=2 bash scripts/pdl1_smoke_all.sh

# Individual stages (each depends on the previous stage's outputs)
bash scripts/pdl1_smoke_rfpeptides.sh      # NUM_DESIGNS, HOTSPOTS, CONTIGS, DIFFUSER_T
bash scripts/pdl1_smoke_mpnn.sh            # SEQS_PER_BB, DESIGN_CHAIN=A, FIXED_CHAINS=B
AFCYC_LIMIT=2 bash scripts/pdl1_smoke_afcyc.sh
uv run python scripts/pdl1_smoke_score.py

# Docker images for the academic stacks
bash scripts/fetch_colabdesign.sh          # prefetch ColabDesign @ v1.1.3 (needed before afcyc build)
bash scripts/build_academic_images.sh afcyc
bash scripts/build_academic_images.sh pyrosetta

uv run ruff check scripts/                 # lint (no ruff config in pyproject; defaults apply)
```

There is **no test suite and no CI**. Verification is by the DoD checks built into each
script: `pdl1_prepare_receptor.py` self-validates hotspots and exits non-zero on failure,
the smoke scripts count their expected PDB outputs, and `pdl1_smoke_score.py` refuses to
report S2 as complete unless AfCyc predictions exist. Re-running the relevant stage is the
regression test.

## The pipeline

```
S0 spec/assets → S1 filter scaffold → S2 smoke (3-stage) → S3 main gen (~10k)
   → S4 structure/physics funnel → S5 blockade + diversity selection → S6 deliver
   → S7 track B (optional/deferred)
```

**S2 is three-stage and mandatory:** RFpeptides (cyclic backbone) → ProteinMPNN (omit C) →
**AfCycDesign** (cyclic complex re-prediction) → uv scoring. A run without AfCyc output is
explicitly *not* a completed S2 — the scripts print a warning and then fail.

Two tracks share the pipeline but never share evaluation criteria:
- **Track A (default):** de novo RFpeptides binders against the 5O45 receptor.
- **Track B:** redesign around the 8ALX macrocycle pharmacophore (pAC65); internal benchmark only.

### Gate policy (encode this, don't reinvent it)

- **Primary gate:** `pd1_jaccard` — overlap between peptide-contact residues and the 4ZQK
  PD-1 footprint. Threshold in [data/targets/PDL1/thresholds_v1.yaml](data/targets/PDL1/thresholds_v1.yaml)
  (currently `T_pd1_jaccard: 0.15`, calibrated so all 5 crystal macro positives pass and all
  negatives fail).
- **Rim residues (A122/A124/A125):** a *soft weighted score* for ranking and quota, **never**
  a `require(D122 and K124 and R125)` conjunction. All 5 calibration positives would fail that
  AND, which is exactly why `forbid_require_all_three: true` is in the thresholds file.
- **Diversity:** `rmsd_vs_8alx` binning with a quota for high-RMSD (opposite-ring / eccentric)
  designs. Selection mix 40% composite / 30% high-rim / 30% high-RMSD.

### Chain and numbering conventions

- RFdiffusion/RFpeptides outputs put the **cyclic binder in chain A and PD-L1 in chain B**.
  The MPNN stage designs chain A and fixes chain B.
- 5O4Y is the exception: PD-L1 is chain B and the peptide is chain D (the receptor cut was
  rewritten B→A for the `PDL1_IgV_5O4Y_chainA.pdb` derivative).
- Residue numbering follows UniProt IgV numbering as in 5O45/7OUN chain A (18–134), so
  hotspot strings like `A56,A123,A115,A113,A66,A121` mean the same thing across structures.
- `PDL1_IgV_5O45_A18-134.pdb` is **generated** — never hand-edit it; it is also the
  repro source for the receptor prep report. Note the `.WRONG_peptide.bak` file next to it in
  the PDB dir if you go looking for history.

## Script roles

- [scripts/pdl1_hotspot_consensus.py](scripts/pdl1_hotspot_consensus.py) — cross-structure
  hotspot consensus (5 macrocycles + 4ZQK + antibody epitopes; includes ncAA HETATM).
  Note: this one script has an **absolute hardcoded `ROOT`** rather than deriving it from
  `__file__` like every other script.
- [scripts/pdl1_prepare_receptor.py](scripts/pdl1_prepare_receptor.py) — crop 5O45 to A18–134,
  validate that all hotspots are present with the expected residues.
- [scripts/pdl1_blockade_metrics.py](scripts/pdl1_blockade_metrics.py) — the shared metrics
  library. `score_complex()` is imported by the smoke scorer, so metric definitions
  (`pd1_jaccard`, `rim_*`, `rim_score`, `rmsd_vs_8alx`) are single-sourced here.
- [scripts/pdl1_calibrate_filters.py](scripts/pdl1_calibrate_filters.py) — positive/negative
  crystal controls; **writes `thresholds_v1.yaml`**. Calibrate before S3, not after.
- [scripts/pdl1_afcyc_predict.py](scripts/pdl1_afcyc_predict.py) — runs *inside* the
  `pdl1-afcyc` image (ColabDesign + JAX). Applies the cyclic relative-positional encoding at
  runtime (`add_cyclic_offset`, `offset_type=2`); this is not baked into the image.
- `pdl1_smoke_*.sh` — Docker orchestration. They mount only what each stage needs (cropped
  receptor as `/inputs/receptor.pdb`, weights as read-only `/models`) and write outputs under
  `results/PDL1/smoke/`. Logs land in `results/PDL1/smoke/logs/`.

## Docker images

Both academic images pin heavily on purpose; the rationale is recorded in
[models/SOURCES.md](models/SOURCES.md) and mirrors comments in the Dockerfiles.
- `pdl1-afcyc` ([docker/afcyc/Dockerfile](docker/afcyc/Dockerfile)): CUDA 12.1 cudnn8 base,
  Python 3.10.15, `jax[cuda12]==0.4.30`, `dm-haiku==0.0.12`, `numpy==1.26.4`. JAX is capped
  at 0.4.30 because dm-haiku ≥0.0.17 needs `jax.interpreters.xla.xe`, which jax ≥0.5 removed.
  AF2 params and ColabDesign are not baked in / not cloned in-Docker respectively — params are
  mounted at runtime and ColabDesign must be host-prefetched (`fetch_colabdesign.sh`).
- `pdl1-pyrosetta` ([docker/pyrosetta/Dockerfile](docker/pyrosetta/Dockerfile)): not needed for
  the S2 smoke path; only for the S4 Rosetta physics stage. Uses a host-prefetched ~1.5 GiB
  wheel via a BuildKit bind-mount when present, else downloads from quarterly find-links.

Rebuilding an image invalidates prior smoke results — re-run
`AFCYC_LIMIT=2 bash scripts/pdl1_smoke_afcyc.sh` as the regression check.

## Assets and artifacts

Large binaries (`*.pt`, AF2 `params/`, PyRosetta `.whl`, `third_party/ColabDesign/`) are
gitignored but expected to be **on disk** for GPU runs — `bootstrap_local_assets.sh` restores
them from documented source paths. This directory is not currently a git repository.

`results/PDL1/` holds the generated evidence and is the reference for what "done" looks like
per stage (`receptor_prep_report.json`, `filter_calibration/`, `smoke/SMOKE_REPORT.md`).
`data/targets/PDL1/` holds the inputs of record. Do not hand-write results there; regenerate
them by re-running the owning script.
