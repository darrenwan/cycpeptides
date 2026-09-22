# PD-L1 过滤阈值校准报告（S1c）

> 日期：2026-09-22  ·  接触阈值：4.5 Å  ·  阈值文件：`data/targets/PDL1/thresholds_v1.yaml`

## 结论

- **主门控**：`pd1_jaccard >= 0.15`
- 阳性通过率：100%（min jaccard=0.2917）
- 阴性落选率：100%（max jaccard=0.0）
- 可分性：是
- **rim AND 诊断**：5/5 个阳性宏环 **不过** `D122∧K124∧R125` → 禁止作默认硬门

## 阳性宏环

| tag | pd1_jaccard | rim_D122 | rim_K124 | rim_R125 | rim_score | rim_all_three |
| --- | --- | --- | --- | --- | --- | --- |
| 7OUN_peptide104 | 0.6667 | 12 | 3 | 0 | 0.2437 | 0 |
| 5O45_peptide57 | 0.4583 | 9 | 0 | 0 | 0.1406 | 0 |
| 8ALX_pAC65 | 0.44 | 0 | 0 | 11 | 0.1719 | 0 |
| 6PV9_macro | 0.375 | 0 | 0 | 0 | 0.0 | 0 |
| 5O4Y_peptide71_BD | 0.2917 | 0 | 0 | 0 | 0.0 | 0 |

## 阴性对照

| tag | pd1_jaccard | rim_D122 | rim_K124 | rim_R125 | rim_score |
| --- | --- | --- | --- | --- | --- |
| decoy_7OUN_pep_plus50A | 0.0 | 0 | 0 | 0 | 0.0 |
| decoy_5O45_pep_plus50A | 0.0 | 0 | 0 | 0 | 0.0 |
| decoy_8ALX_pep_plus50A | 0.0 | 0 | 0 | 0 | 0.0 |
| 5O45_no_ligand | 0.0 | 0 | 0 | 0 | 0.0 |

## 门控策略（冻结 v1）

1. P0：足迹 Jaccard ≥ T
2. rim：分项软分 + 配额（非 AND）
3. 结构/物理阈值待 S2 冒烟后补校准

原始表：`results/PDL1/filter_calibration/calibration_metrics.csv`
