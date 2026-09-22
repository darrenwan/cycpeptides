# S2 冒烟报告（uv + Docker 三阶段）

> 日期：2026-09-22  
> 战役 Python：**仅 uv**（`.venv`，`uv run`）  
> GPU：**Docker** RFpeptides → ProteinMPNN → **AfCycDesign**（`pdl1-afcyc`）

## 结论

| 阶段 | 结果 |
|------|------|
| S1a–c | PASS（uv） |
| S2a RFpeptides | PASS（N=4 cyclic binder） |
| S2b ProteinMPNN | PASS（8 designed PDBs，omit C） |
| S2c **AfCycDesign** | PASS（本轮 `AFCYC_LIMIT=2`；`n_ok=2/2`，pep RMSD≈0.3 Å） |
| S2d uv 评分 | PASS（对 AfCyc 预测结构打分；`smoke_metrics.csv`） |

**S2 三阶段 DoD 已满足**（RFD → MPNN → AfCyc → uv 评分）。  
扩规模：`NUM_DESIGNS=100 bash scripts/pdl1_smoke_all.sh`；AfCyc 全量可不设 `AFCYC_LIMIT`。

## 环境约定

```bash
cd /mnt/data4t/aidd/cycpeptides
uv sync
uv run python scripts/...          # 战役脚本唯一入口
bash scripts/pdl1_smoke_all.sh     # 三阶段 Docker + uv 评分
# 仅补 AfCyc：AFCYC_LIMIT=2 bash scripts/pdl1_smoke_afcyc.sh
```

禁止用系统 Python / conda base 跑 `scripts/pdl1_*.py`。

## 产物路径

- `results/PDL1/smoke/rfpeptides/pdl1_cyc_*.pdb`
- `results/PDL1/smoke/mpnn/designed_pdbs/*.pdb`
- `results/PDL1/smoke/afcyc/*_afcyc_pred.pdb` + `afcyc_summary.json`
- `results/PDL1/smoke/smoke_metrics.csv` / `smoke_report.json`
- 日志：`results/PDL1/smoke/logs/`

## 参数（本轮）

- Receptor：`PDL1_IgV_5O45_A18-134.pdb`
- Contig：`[12-16 A18-134/0]`，`inference.cyclic=True`，`cyc_chains=a`
- Hotspots：`A56,A123,A115,A113,A66,A121`
- MPNN：`SEQS_PER_BB=2`，`omit_AAs=C`
- AfCyc：镜像 `pdl1-afcyc:latest`；params 挂载自 sibling `models/afcyc`；recycles=3；`use_binder_template=1`
