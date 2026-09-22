# Local model / vendored assets (project-owned)

Copied or fetched into this repo for campaign independence. Do not rely on sibling
project paths at runtime.

| Asset | Local path | Source (provenance) | Notes |
|-------|------------|---------------------|-------|
| RFdiffusion Complex_base | `models/rfdiffusion/Complex_base_ckpt.pt` | `/mnt/data4t/aidd/models/rfmodels/` (RosettaCommons RFdiffusion binder ckpt) | Track A binder + hotspots |
| AfCyc / AF2 params | `models/afcyc/params/` | cursor_project `.../models/afcyc/params` (**alphafold_params_2022-12-06**) | Mount into Docker `pdl1-afcyc` as `/models/afcyc` |
| ColabDesign | `third_party/ColabDesign/` | `https://github.com/sokrypton/ColabDesign.git` **tag `v1.1.3`** | `bash scripts/fetch_colabdesign.sh`; pin file `.colabdesign_pin` |
| ProteinMPNN | `third_party/ProteinMPNN/` | `/mnt/data4t/kg_project/ProteinMPNN` (Dauparas et al.) | Scripts + `vanilla_model_weights` |

## AfCyc image pins (`pdl1-afcyc`)

Built by `bash scripts/build_academic_images.sh afcyc` from `docker/afcyc/Dockerfile`:

| Component | Pin |
|-----------|-----|
| Base | `nvidia/cuda:12.1.1-cudnn8-runtime-ubuntu22.04` (DaoCloud mirror OK) |
| CPython | 3.10.15 (uv-managed) |
| jax / jaxlib | `jax[cuda12]==0.4.30` (GCS; pulls cudnn9 pip wheels onto cudnn8 base) |
| dm-haiku | 0.0.12 (`--no-deps`, avoid ≥0.0.17 + jax≥0.5 break) |
| flax | 0.8.5 |
| numpy | 1.26.4 (`<2` for AF2/ColabDesign stability) |
| biopython | 1.84 |
| ColabDesign | host `third_party/ColabDesign` @ **v1.1.3** (`uv pip install --no-deps`; no git-in-Docker) |
| uv HTTP timeout | 7200s (large CUDA wheels) |

Cyclic encoding for AfCycDesign is applied at runtime in `scripts/pdl1_afcyc_predict.py`
(`add_cyclic_offset`, `offset_type=2`), not baked into the image.

Large binaries are gitignored; keep them on disk for local runs.

```bash
bash scripts/bootstrap_local_assets.sh   # weights + ProteinMPNN
bash scripts/fetch_colabdesign.sh        # ColabDesign source for Docker build
bash scripts/build_academic_images.sh afcyc
```
