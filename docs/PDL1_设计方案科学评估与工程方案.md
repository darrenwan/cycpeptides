# PD-L1 大环肽干实验方案：科学评估与可落地工程方案

> **日期**：2026-09-22  
> **范围**：正式写流水线代码前的方案审查  
> **依据**：`docs/PDL1_大环肽_战役规格.md`、`data/targets/PDL1/design_spec.yaml`、`results/PDL1/hotspot_consensus.*`、本地 PDB 复核  
> **项目边界**：纯干实验；止于计算候选池

---

## 0. 总判（Verdict）

| 维度 | 判定 | 一句话 |
|------|------|--------|
| 靶点与机制 | **合理** | PD-L1 IgV 面阻断有多枚宏环共晶 + 公开 BMS 分子支撑 |
| 受体/轨道分工 | **合理** | 5O45 生成 / 8ALX 药效团 / 7OUN 解读 已拆开，避免混用 |
| Hotspot 生成集 | **合理** | 6 残基共识 ∩ PD-1 足迹；本地编号复核通过 |
| 工具主线 | **合理** | RFpeptides → MPNN → AfCyc/HighFold → Rosetta 符合 2025 证据强度 |
| 阻断「硬 AND」过滤 | **需修正** | 文献阳性宏环几乎都不同时接触 D122∧K124∧R125 |
| 差异化叙事 | **方向对、执行需软化** | 极性边缘应作**加分/配额**，不宜作一票否决 |
| 工程就绪度 | **可开工，但先 S1 校准** | 缺过滤实现、缺阳性校准、缺受体裁剪与冒烟门控 |

**结论**：科学大方向可进入工程实现；**唯一必须在写大规模生成代码前改掉的设计点**是把 D122/K124/R125 从「三者同时硬门控」改为「分项遮挡分 + 校准后的软门控/配额」。其余为优先级与实现顺序问题。

---

## 1. 科学评估

### 1.1 已确认合理之处（有本地/文献证据）

1. **靶点选择**  
   Holak/BMS 等多枚 PD-L1–宏环共晶（5O45、7OUN、8ALX、5O4Y、6PV9）表明该平坦 CC′FG 面可用大环肽遮挡；与抗体/PD-1 同面收敛的判断与结构一致。

2. **单体 IgV 面阻断 vs 小分子二聚**  
   明确排除 5N2F 等二聚口袋，机制上正确；避免假差异化。

3. **轨道 A/B 拆分**  
   - 轨道 A：de novo canonical，IP 友好  
   - 轨道 B：8ALX 药效团内部对标  
   避免把 BMS 几何烙进主生成受体，逻辑清楚。

4. **受体分辨率 vs 结合模式**  
   本地复核（修正 Kabsch）：5O45 vs 7OUN / 8ALX 界面 Cα RMSD ≈ **0.30 / 0.32 Å**；全共享 CA ≈ **0.71 / 0.83 Å**。战役规格中「~0.4 Å」量级成立。生成用最高分辨率受体、过滤阶段交叉构象，合理。

5. **Hotspot 编号**  
   切链受体上 Y56/Q66/R113/M115/A121/Y123 与 D122/K124/R125 均存在且氨基酸正确（5O45/7OUN/8ALX/5O4Y）。

6. **Hotspot ≤6、差异化放过滤**  
   符合 RFdiffusion hotspot 用法：条件生成宜少而准，扩展接触勿塞进 `hotspot_res`。

7. **代理读出层级（足迹 > 贴面）**  
   对纯干实验正确：iPAE/ddG 不能代替「是否挡住 PD-1 足迹」。

8. **先校准再放行**  
   用公开宏环作阳性对照校准分位数，是平坦 PPI 上避免自欺的必要步骤。

### 1.2 必须修正：阻断硬过滤与阳性对照冲突

对文献复合物，重原子 ≤4.5 Å 接触原子对数（本地计算，2026-09-22）：

| 复合物 | D122 | K124 | R125 | 若要求三者皆接触 |
|--------|------|------|------|------------------|
| 7OUN peptide-104 | 12 | 3（弱） | **0** | **不通过** |
| 5O45 peptide-57 | 9 | **0** | **0** | **不通过** |
| 8ALX pAC65 | **0** | **0** | 11 | **不通过** |
| 6PV9 | 0 | 0 | 0 | **不通过** |
| 5O4Y peptide-71 | 0 | 0 | 0 | **不通过** |
| 4ZQK PD-1（足迹定义） | 19 | 27 | 16 | （配体是 PD-1，非宏环） |

**含义**：当前规格把 D122/K124/R125 写成「阻断硬条件 / 强制遮挡」，若实现为 **AND 门控**，则：

- 战役规定的**阳性骨架一律落选** → 校准逻辑自相矛盾；  
- 会系统性偏向「比已知活性宏环更贪心的极性边缘覆盖」，假阴性极高；  
- 「差异化」被误写成「活性必要条件」。

**修正建议（写入规格后再编码）**

| 层级 | 规则 |
|------|------|
| **P0 门控** | PD-1 足迹重叠（残基 Jaccard 或溶剂可及体积遮挡）过校准阈值 |
| **P0.5 差异化评分** | `rim_score = f(接触 D122, K124, R125)`，分项计分，**不要求三者同时** |
| **配额（轨道 A）** | 优先池中强制保留高 `rim_score` 与高 vs-8ALX RMSD 的簇，各占一定比例 |
| **阳性期望** | 7OUN/5O45 应过 P0；8ALX 在 R125 分项上应高；**不要求**阳性过「三残基全中」 |

差异化报告中的 Lys124「宏环空白」仍成立，但应解读为 **机会位点**，不是活性充分必要条件。

### 1.3 中等风险 / 设计缺口

| 问题 | 严重度 | 说明 | 对策 |
|------|--------|------|------|
| Canonical 头尾环 vs 文献 ncAA 活性宏环 | 中 | 轨道 A 化学空间更难；成功率预期低于 ncAA 对标 | 接受；冒烟与对照看「能否过过滤」而非「能否复现 nM」 |
| 平坦界面假阳性 | 高（已知） | 文献亦称 PD-L1 难；置信度≠亲和 | 双模型 + 足迹 + 对照；报告禁止写 Kd |
| Group II 仅靠过滤保留 | 中 | 高 vs-8ALX RMSD ≠ 5O4Y 反向族 | S3 可选小批量 5O4Y 受体生成，或 motif 反向种子 |
| `contig A18-134` vs 5O45 链 17–145 | 中 | 受体多出 17 与 135–145 | **生成前裁剪**到 18–134（或与 contig 严格一致） |
| 双模型全量过贵 | 中 | 10k×4 全跑 AfCyc+HighFold+Rosetta 不现实 | 漏斗：粗筛 → 精筛 → 交叉构象 |
| vs-8ALX RMSD 多样性 | 低–中 | 高 RMSD 有利 IP，不利药效团迁移 | 分箱配额，勿只留高 RMSD |
| 净电荷 0–+2、omit Cys | 低 | 启发式合理 | 作软约束/标注 |
| 鼠 PD-L1 分簇 | 低 | 仅标注即可 | S6 可选列 |
| RFpeptides 未在公开材料中以 PD-L1 为成功案例主打 | 中 | 方法外推 | 用本靶点冒烟 + 阳性宏环校准界定「计算成功」 |

### 1.4 能力边界（必须写进交付报告）

1. 本战役产出是**排序代理候选池**，不是实测阻断剂清单。  
2. 打分与 Kd 相关性弱（对接/能量文献共识）；校准只保证「已知结构不至于被自己的过滤器杀掉」。  
3. 轨道 A canonical 不能复现 pAC65 酰化-Trp 药效团；差异化是设计选择，不是效力保证。

### 1.5 科学结论摘要

- **可以开发**：靶点、结构资产、hotspot、工具栈、双轨与纯干实验边界均成立。  
- **开发前改规格**：阻断残基从硬 AND → 分项评分 + 足迹主门控 + 多样性/差异化配额。  
- **开发顺序**：S1 过滤校准 → S2 冒烟 → 再 S3 大规模；禁止颠倒。

---

## 2. 可落地工程方案

### 2.1 目标与非目标

| 要做 | 不做 |
|------|------|
| 可复现流水线：规格 YAML → 生成 → 过滤 → 候选包 | 合成、assay、给药、临床 |
| 过滤阈值版本化与阳性/阴性校准报告 | 把分位数写成「预测 IC50」 |
| 冒烟门控后再扩到 ~10k | 未校准就全量 GPU 烧预算 |

### 2.2 推荐仓库布局（增量）

```text
cycpeptides/
├── data/targets/PDL1/
│   ├── design_spec.yaml      # 科学规格（已有）
│   ├── roadmap.yaml          # 阶段状态机（新增）
│   ├── hotspots*.json
│   └── pdb/                  # 原始 + 切链；另增 cropped 受体
├── scripts/
│   ├── pdl1_hotspot_consensus.py          # 已有
│   ├── extract_pdl1_5o4y_receptor.py      # 已有
│   ├── pdl1_prepare_receptor.py           # 新增：裁剪 18–134、校验 hotspot
│   ├── pdl1_blockade_metrics.py           # 新增：足迹 / rim_score / vs-8ALX
│   ├── pdl1_calibrate_filters.py          # 新增：阳/阴性对照跑分
│   ├── pdl1_select_candidates.py          # 新增：聚类 + 配额挑选
│   └── pipelines/                         # 可选：Nextflow/Snakemake
├── results/PDL1/
│   ├── filter_calibration/
│   ├── smoke/
│   ├── rfpeptides/  mpnn/  structure_filter/
│   └── candidates_round1/
└── docs/
    ├── PDL1_大环肽_战役规格.md
    └── PDL1_设计方案科学评估与工程方案.md   # 本文
```

### 2.3 过滤逻辑（编码契约）

对每个「肽–PD-L1」复合物模型输出一行指标：

| 字段 | 定义 | 用途 |
|------|------|------|
| `pd1_jaccard` | 肽接触残基 ∩ 4ZQK PD-1 足迹 / 并集 | **主门控**（阈值由校准定） |
| `rim_D122/K124/R125` | 各残基接触原子对或距离 proxy | 分项分 |
| `rim_score` | 分项归一化加权和（权重可配） | 排序 + 配额 |
| `ipae` / `iptm` / `iface_plddt` | AfCyc 或 HighFold | 结构置信 |
| `pep_ca_rmsd` | 设计 vs 预测 | 自洽 |
| `ddg` / `cms` / `sap` | Rosetta | 物理软过滤 |
| `rmsd_vs_8alx` | 肽骨架对齐后 Cα RMSD | 多样性分箱 |
| `charge` / `cycpeptmp` | 可选 | 弱权重 |

**默认门控伪代码**（阈值数字须经校准改写）：

```text
pass_P0 = pd1_jaccard >= T_jaccard          # 校准所得
pass_struct = pep_ca_rmsd <= 2.0 and ipae in top_q
pass_physics = ddg <= median_pool and sap not extreme

eligible = pass_P0 and pass_struct and pass_physics

# 挑选 24–48：在 eligible 内
#   40% 按 (pd1_jaccard, rim_score, -ipae) 综合排序
#   30% 高 rim_score 配额（哪怕综合分略低）
#   30% 高 rmsd_vs_8alx 配额（反向/偏心）
```

**禁止**：`require(D122 and K124 and R125)` 作为默认。

### 2.4 阶段工程计划（与 §10 对齐，含门控）

| 阶段 | 工程任务 | 完成定义（DoD） | 预估 |
|------|----------|-----------------|------|
| **S1a** | `pdl1_prepare_receptor.py`：裁剪 5O45→18–134；校验 hotspot；写 `PDL1_IgV_5O45_A18-134.pdb` | contig 与文件残基集合一致；hotspot 全在 | 0.5 d |
| **S1b** | `pdl1_blockade_metrics.py`：足迹 Jaccard、rim 分项、vs-8ALX RMSD | 对 7OUN/8ALX/5O45 晶体复合物跑通并落盘 | 1–2 d |
| **S1c** | `pdl1_calibrate_filters.py`：阳/阴对照 | 报告：阳性过 P0 率、阴性落选率；写出 `thresholds_v1.yaml` | 1 d |
| **S1 门控** | 人工读校准报告 | 阳性宏环 **不得**被默认门控杀光；否则降 T 或改 rim | — |
| **S2** | RFpeptides → MPNN → **AfCycDesign** 冒烟 | 三阶段全通；AfCyc 预测 PDB + 指标表非空；路径文档化 | 2–5 d（含环境） |
| **S3** | 扩至 ~10k 骨架 ×4 seq | 仅在 S2 通过后；断点续跑；parquet 清单 | GPU 墙钟视集群 |
| **S4** | 漏斗：Top 结构置信 → Rosetta → 交叉构象 | 明确每级保留比例（如 5%→1%→交叉 200） | 与 S3 重叠 |
| **S5** | 套用 `thresholds_v1` + 配额挑选 | 分箱直方图 + 入池名单 | 1 d |
| **S6** | `candidates_round1/` 打包 | 序列、PDB、指标、溯源 JSON、一页报告 | 0.5–1 d |
| **S7** | 轨道 B（可选） | 8ALX motif；与轨道 A 指标表隔离 | 后置 |

### 2.5 计算漏斗（建议默认）

```text
10k backbones × 4 seq
    → 廉价预筛（可选：快速 AF/开放式 iPAE 或仅 MPNN score）保留 ~10–20%
    → AfCyc 或 HighFold 单模型 → 保留 Top ~5%
    → 第二模型共识 + Rosetta → ~200–500
    → 足迹 P0 + rim/多样性配额 → 24–48
    → 仅对最终池做 7OUN/8ALX 交叉构象复评
```

全量双模型从第一天起跑 **不推荐**。

### 2.6 环境与依赖（工程前置）

#### 战役 Python：**仅 uv**

禁止使用系统 Python / conda base 运行本仓库 `scripts/pdl1_*.py`。

| 项 | 约定 |
|----|------|
| 工具 | [uv](https://github.com/astral-sh/uv) |
| 清单 | `pyproject.toml`（`requires-python >=3.11`） |
| 锁文件 | `uv.lock` |
| 虚拟环境 | 项目根目录 `.venv`（`uv sync` 创建） |
| 入口 | `uv run python scripts/...` |
| 核心依赖 | biopython、numpy、pandas、pyyaml |

```bash
cd /mnt/data4t/aidd/cycpeptides
uv sync
uv run python scripts/pdl1_prepare_receptor.py
uv run python scripts/pdl1_blockade_metrics.py
uv run python scripts/pdl1_calibrate_filters.py
bash scripts/pdl1_smoke_all.sh
```

#### GPU / 重依赖栈：Docker（与 uv 分工）

| 组件 | 备注 |
|------|------|
| RFdiffusion + cyclic（RFpeptides） | Docker `rosettacommons/rfdiffusion:latest`；权重默认 `models/rfdiffusion/`；`inference.cyclic=True` |
| ProteinMPNN | 在 RFD 镜像内用 torch 跑；代码/权重默认 `third_party/ProteinMPNN/`；战役侧 omit C |
| AfCycDesign | **S2 必跑**（Docker `pdl1-afcyc`）；AF2 params 默认 `models/afcyc/`；环化复合物复核；不进 uv `.venv` |
| HighFold | 可选交叉验证；S4 可并列 |
| Rosetta | FastRelax + 界面指标；学术许可；后续阶段 |
| （可选）Nextflow/Snakemake | S3 以后再上 |

机器可读：`data/targets/PDL1/roadmap.yaml` → `runtime`；`design_spec.yaml` → `structures.runtime`。

**S1 不依赖** RFdiffusion / AfCyc GPU；可在 uv CPU 环境完成过滤脚手架与晶体对照校准。  
**S2 依赖** RFpeptides + MPNN + **AfCyc** 三阶段；缺 AfCyc 不得宣称 S2 完成。

### 2.7 立即开工的最小增量

1. ~~改规格措辞：硬 AND → rim 评分~~（已完成）  
2. ~~`uv sync` + 受体裁剪 / blockade / 校准~~（S1 已完成）  
3. ~~Docker RFpeptides + MPNN + **AfCycDesign**~~（S2 三阶段已通；可扩 `NUM_DESIGNS` / 取消 `AFCYC_LIMIT`）  
4. 下一优先：S3 大规模生成（仅在 S2 通过后）。

---

## 3. 对现行文档的修订清单

| 文档位置 | 修订 |
|----------|------|
| 战役规格 §3.3 / 过滤表 | 「阻断硬条件」→「差异化 rim 分项 + 配额；主门控为 PD-1 足迹」 |
| 战役规格 §6.3 | P0 明确为足迹；rim 为 P0.5/配额 |
| `design_spec.yaml` `biology_proxy` | `blockade_hard_residues` → `blockade_rim_score_residues` + `gate: pd1_footprint` |
| §10 S1 | 增加受体裁剪；校准失败则禁止进 S3 |
| 差异化报告 | Lys124 等改为「机会位点 / 配额」，删除「强制遮挡三者」语气 |

（修订可与 S1 编码同一 PR/同一会话完成。）

---

## 4. 溯源（本评估用到的本地证据）

- 界面 RMSD：5O45 vs 7OUN / 8ALX ≈ 0.30 / 0.32 Å（共享界面 CA，Kabsch）。  
- 硬过滤冲突表：见 §1.2（≤4.5 Å 重原子接触对）。  
- Hotspot 存在性：切链 PDB 残基复核。  
- 5O45 超出 contig 的残基：17，135–145。  
- 方法学主线：RFpeptides *Nat Chem Biol* 2025 doi:10.1038/s41589-025-01929-w。

---

## 5. 一页决策

| 决策 | 选择 |
|------|------|
| 是否按现方案开发？ | **是，但先改 rim 门控定义** |
| 是否上 Kiro 全套？ | **否**（见前次讨论）；用 `roadmap.yaml` 即可 |
| 第一行代码写什么？ | 受体裁剪 + `blockade_metrics` + 晶体校准，**不是** 10k RFdiffusion |
| 战役成功标准（干实验） | 校准合理 + 冒烟通路通 + 交付 24–48 条带完整溯源的候选包 |
