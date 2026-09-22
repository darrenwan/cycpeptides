#!/usr/bin/env python3
"""
PD-L1 hotspot consensus + differentiation opportunity analysis.

Inputs: co-crystal PDBs under data/targets/PDL1/pdb/
Outputs:
  results/PDL1/hotspot_consensus.csv
  results/PDL1/hotspot_consensus.json
  results/PDL1/differentiation_opportunities.md
  data/targets/PDL1/hotspots_ranked.json

Method: heavy-atom contacts ≤4.5 Å; ligand includes non-water HETATM (ncAA).
Date: 2026-09-22.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

ROOT = Path("/mnt/data4t/aidd/cycpeptides")
PDB_DIR = ROOT / "data/targets/PDL1/pdb"
OUT_DIR = ROOT / "results/PDL1"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CUTOFF = 4.5
PARSER = PDBParser(QUIET=True)
WATER = {"HOH", "WAT", "H2O", "DOD"}

# Fixed complexes (5O4Y paired automatically)
COMPLEXES = [
    ("7OUN", "A", ["B"], "macrocycle", "peptide-104"),
    ("5O45", "A", ["B"], "macrocycle", "peptide-57"),
    ("8ALX", "A", ["B"], "macrocycle", "pAC65"),
    ("6PV9", "A", ["B"], "macrocycle", "macrocycle_6PV9"),
    ("4ZQK", "A", ["B"], "pd1", "PD-1"),
]


def polymer_residues(chain):
    return [r for r in chain if r.id[0] == " "]


def ligand_residues(chain):
    return [r for r in chain if r.resname not in WATER]


def heavy_atoms(residue):
    return [a for a in residue if a.element != "H"]


def contact_counts(structure, pdl1_chain_id, partner_ids, cutoff=CUTOFF):
    model = structure[0]
    pdl1 = model[pdl1_chain_id]
    partners = []
    for cid in partner_ids:
        if cid in model:
            partners.extend(ligand_residues(model[cid]))
    if not partners:
        return Counter()
    partner_coords = np.array(
        [a.coord for r in partners for a in heavy_atoms(r)]
    )
    if len(partner_coords) == 0:
        return Counter()
    counts = Counter()
    for res in polymer_residues(pdl1):
        atoms = heavy_atoms(res)
        if not atoms:
            continue
        coords = np.array([a.coord for a in atoms])
        d2 = ((coords[:, None, :] - partner_coords[None, :, :]) ** 2).sum(axis=2)
        n_pairs = int((d2 <= cutoff**2).sum())
        if n_pairs > 0:
            counts[res.id[1]] += n_pairs
    return counts


def aa_map(structure, chain_id):
    return {
        res.id[1]: res.resname
        for res in polymer_residues(structure[0][chain_id])
    }


def auto_pair_5o4y(structure):
    model = structure[0]
    pdl1_ids = [c for c in ("B", "C", "E") if c in model]
    pep_ids = [c for c in ("A", "D", "F") if c in model]

    def com(chain_id):
        coords = []
        for r in ligand_residues(model[chain_id]):
            if "CA" in r:
                coords.append(r["CA"].coord)
            else:
                atoms = heavy_atoms(r)
                if atoms:
                    coords.append(np.mean([a.coord for a in atoms], axis=0))
        return np.mean(coords, axis=0) if coords else None

    pairs, used = [], set()
    for p in pdl1_ids:
        cp = com(p)
        best, best_d = None, 1e9
        for q in pep_ids:
            if q in used:
                continue
            cq = com(q)
            if cp is None or cq is None:
                continue
            d = float(np.linalg.norm(cp - cq))
            if d < best_d:
                best, best_d = q, d
        if best is not None:
            used.add(best)
            pairs.append((p, best, best_d))
            print(f"5O4Y pair PD-L1 {p} ↔ peptide {best} (COM {best_d:.1f} Å)")
    return pairs


def probe_antibody_chains():
    specs = []
    for pdb, name in [
        ("5XXY", "atezolizumab"),
        ("5GRJ", "avelumab"),
        ("5X8M", "durvalumab"),
    ]:
        path = PDB_DIR / f"{pdb}.pdb"
        if not path.exists():
            continue
        s = PARSER.get_structure(pdb, str(path))
        model = s[0]
        pdl1_id = None
        candidates = []
        for chain in model:
            am = aa_map(s, chain.id)
            n = len(polymer_residues(chain))
            if (
                am.get(56) == "TYR"
                and am.get(115) == "MET"
                and am.get(123) == "TYR"
            ):
                candidates.append((chain.id, n))
        # Prefer IgV-sized chain; else accept full ectodomain
        for cid, n in candidates:
            if 90 <= n <= 150:
                pdl1_id = cid
                break
        if pdl1_id is None and candidates:
            pdl1_id = candidates[0][0]
        if pdl1_id is None:
            print(f"WARN: no PD-L1 chain in {pdb}")
            continue
        ab_chains = [
            c.id
            for c in model
            if c.id != pdl1_id and len(polymer_residues(c)) >= 50
        ]
        specs.append((pdb, pdl1_id, ab_chains, "antibody", name))
        print(f"Ab {pdb}/{name}: PD-L1={pdl1_id} partners={ab_chains}")
    return specs


def write_diff_report(rows, hotspots_ranked, n_macro, n_ab):
    core = [r for r in rows if r["tier"] == "core_required"]
    strong = [r for r in rows if r["tier"] == "core_strong"]
    extended = [r for r in rows if r["tier"] == "extended"]
    rf_nums = [
        int(x[1:]) for x in hotspots_ranked.get("recommended_rf_string_tier1", [])
    ]
    rf_rows = sorted(
        [r for r in rows if r["resnum"] in rf_nums],
        key=lambda r: rf_nums.index(r["resnum"]),
    )
    pd1_rich = sorted(
        [r for r in rows if r["pd1_contact"]], key=lambda x: -x["pd1_atom_pairs"]
    )
    macro_under_pd1 = [
        r
        for r in pd1_rich
        if r["macro_n_structures"] <= 2 and r["pd1_atom_pairs"] >= 5
    ]
    ab_rim = [
        r
        for r in rows
        if r["antibody_n_structures"] >= max(n_ab - 1, 1)
        and r["macro_n_structures"] <= 2
        and r["pd1_contact"]
    ]
    r125 = next((r for r in rows if r["resnum"] == 125), None)
    r124 = next((r for r in rows if r["resnum"] == 124), None)
    r122 = next((r for r in rows if r["resnum"] == 122), None)
    role_map = {
        56: "中央芳香锚；宏环/抗体/PD-1 共有",
        123: "PD-1 最密接触；第二芳香锚",
        115: "疏水沟 + S/π",
        113: "静电/堆叠边缘",
        66: "极性锚定；PD-1 强；pAC65 最密宏环接触之一",
        121: "浅沟底；PD-1 中等密度",
        54: "沟缘疏水；宏环高共识",
        58: "氧–π / 沟缘",
    }
    rf = hotspots_ranked.get("recommended_rf_string_tier1", [])
    lines = [
        "# PD-L1 大环肽：Hotspot 共识统计与差异化机会",
        "",
        f"> 日期：{date.today()}  ·  接触阈值：重原子 ≤ {CUTOFF} Å  ·  "
        f"宏环：{n_macro}（7OUN/5O45/8ALX/5O4Y共识/6PV9）  ·  "
        f"抗体：{n_ab}（5XXY/5GRJ/5X8M）  ·  PD-1：4ZQK",
        "",
        "配体原子含 **ncAA HETATM**（否则漏计 pAC65 TrpNAc–Arg125 等）。",
        "",
        "原始表：`results/PDL1/hotspot_consensus.csv`  ·  "
        "`results/PDL1/hotspot_consensus.json`  ·  "
        "`data/targets/PDL1/hotspots_ranked.json`",
        "",
        "---",
        "",
        "## 1. 统计结论：最佳 Hotspot",
        "",
        "评分 = `3×宏环分数 + 2×PD-1命中 + min(密度/20,1.5) + 0.5×抗体分数`；"
        "宏环≥3 且 PD-1 接触对≥15 → 升入 Tier-1（保护 Tyr123）。",
        "",
        "### 推荐 `ppi.hotspot_res`（6 残基）",
        "",
        "| 残基 | 宏环命中 | 均接触对数 | PD-1对数 | Ab命中 | 角色 |",
        "|------|----------|------------|----------|--------|------|",
    ]
    for r in rf_rows:
        lines.append(
            f"| **{r['label']}** | {r['macro_n_structures']}/{n_macro} | "
            f"{r['macro_mean_pairs_when_hit']} | {r['pd1_atom_pairs']} | "
            f"{r['antibody_n_structures']}/{n_ab} | "
            f"{role_map.get(r['resnum'], '')} |"
        )
    lines += [
        "",
        "```text",
        f"ppi.hotspot_res={rf}",
        "```",
        "",
        "### 全量 Tier-1 / Tier-2 / Extended",
        "",
        "**core_required：** "
        + ", ".join(f"{r['label']}({r['macro_n_structures']}/{n_macro})" for r in core),
        "",
        "**core_strong：** " + (", ".join(r["label"] for r in strong) or "—"),
        "",
        "**extended：** " + (", ".join(r["label"] for r in extended[:12]) or "—"),
        "",
        "### 与原战役规格对照",
        "",
        "| 原规格 | 本轮 | 建议 |",
        "|--------|------|------|",
        "| Y56/M115/Y123/R113/A121/Q66 | 高度重现 | **保留为默认 RF 集** |",
        "| — | I54 宏环共识升至 Tier-1 | 过滤加分或替补 A121 |",
        "| — | E58/N63 宏环常见、PD-1 弱 | 勿进 hotspot_res |",
        "",
        "**原则**：`hotspot_res` 不超过 6–7；差异化放在阻断足迹过滤与构象/化学族。",
        "",
        "---",
        "",
        "## 2. 差异化机会（优先级）",
        "",
        "抗体 / BMS 宏环 / 部分小分子在同一 CC′FG 面收敛。"
        "再盖 Tyr56–Tyr123 难形成科学或 IP 差异；应差在极性边缘化学、环取向、偏心遮挡与格式。",
        "",
        "### ★★★ A. 极性边缘「方向柄」重构（相对 BMS TrpNAc–Arg 盐桥）",
        "",
        "pAC65 用酰化 Trp (V6W) 伸向 Arg113/Arg125——文献标志药效团。"
        "修正 HETATM 后 Arg125 确被 pAC65 接触；Lys124 仍是宏环相对空白。",
        "",
    ]
    if r125:
        lines.append(
            f"- **Arg125**：宏环 {r125['macro_n_structures']}/{n_macro}，"
            f"PD-1 {r125['pd1_atom_pairs']}，Ab {r125['antibody_n_structures']}/{n_ab}"
        )
    if r124:
        lines.append(
            f"- **Lys124**：宏环 {r124['macro_n_structures']}/{n_macro}，"
            f"PD-1 {r124['pd1_atom_pairs']}（PD-1 第二密）——**最佳宏环空白**"
        )
    if r122:
        lines.append(
            f"- **Asp122**：宏环 {r122['macro_n_structures']}/{n_macro}，"
            f"PD-1 {r122['pd1_atom_pairs']} —— canonical Asp/Arg 锁适宜位点"
        )
    lines += [
        "",
        "PD-1 强 / 宏环弱（≤2/5）：",
        "",
        "| 残基 | PD-1 接触对 | 宏环命中 |",
        "|------|-------------|----------|",
    ]
    for r in macro_under_pd1[:10]:
        lines.append(
            f"| {r['label']} | {r['pd1_atom_pairs']} | "
            f"{r['macro_n_structures']}/{n_macro} |"
        )
    lines += [
        "",
        "**动作**：生成仍用 6 残基 hotspot；过滤以 **PD-1 足迹为主门控**，"
        "对 **K124/R125/D122 做 rim 分项分 + 配额**（禁止三者 AND 硬门）；"
        "轨道 A 用 Asp/Glu/Arg 溶剂侧做非酰化-Trp 极性锁（对标 104 而非 BMS）；"
        "canonical → 与 BMS ncAA **化学空间天然分流**。",
        "",
        "### ★★★ B. Group II / 反向环走向（5O4Y 族）",
        "",
        "peptide-71 相对 peptide-57 反向环走向；公开专利/文献多在 Group I / BMS hairpin。",
        "",
        "**动作**：按 vs 8ALX/7OUN 肽 Cα RMSD 分箱，**配额保留**高 RMSD 簇；"
        "可选小批量 5O4Y 受体交叉生成。",
        "",
        "### ★★ C. 偏心遮挡 + 可缀合溶剂面",
        "",
        "相对 pAC65「最小抗体」：偏心斑块优先挡 Y123/D122/K124/Q66，"
        "允许 M115 半占有，**以 PD-1 足迹主门控 + rim 配额作干实验放行**；"
        "溶剂面预留 Lys 做缀合格式标注。",
        "",
    ]
    if ab_rim:
        lines.append(
            "抗体重、宏环偏轻 rim： "
            + ", ".join(r["label"] for r in ab_rim[:10])
        )
        lines.append("")
    lines += [
        "### ★★ D. 人特异 vs 人鼠保守性两条序列簇（in silico 标注）",
        "",
        "鼠界面非全保守；按保守性分簇标注，避免全池只服务一种假设；"
        "**非动物实验计划**。",
        "",
        "### ★ E. 假差异化（不做）",
        "",
        "| 陷阱 | 原因 |",
        "|------|------|",
        "| 小分子二聚口袋 | 机制不同 |",
        "| 复制 pAC65 双 Trp + 酰化盐桥 | IP 撞车 |",
        "| 无极性补偿纯疏水斑 | 聚集 / 可制造性启发式失败 |",
        "| hotspot_res ≥10 | 多样性坍缩 |",
        "",
        "---",
        "",
        "## 3. 战役级落地（纯干实验）",
        "",
        "| 层级 | 内容 |",
        "|------|------|",
        f"| 生成 | `{rf}` |",
        "| 过滤 | PD-1 足迹主门控；**rim 分项分 + 配额**（非 AND） |",
        "| 多样性 | vs 8ALX RMSD 分箱配额；保留反向/偏心 |",
        "| 化学 | 轨道 A = canonical；ncAA → 轨道 B |",
        "| 代理读出 | 足迹遮挡 > 单纯贴面接触 |",
        "| 交付 | 优先候选池 + 指标表 + 复合物 PDB（无湿测） |",
        "",
        "## 4. 溯源",
        "",
        "- 脚本：`scripts/pdl1_hotspot_consensus.py`",
        "- 结构：`data/targets/PDL1/pdb/`",
        "- 评估：`docs/PDL1_设计方案科学评估与工程方案.md`",
        "- 文献：doi:10.1186/s12943-023-01853-4；"
        "doi:10.3390/molecules26164848；"
        "EATRIS/ACS “One Face, Three Solutions”",
        "",
    ]
    (OUT_DIR / "differentiation_opportunities.md").write_text("\n".join(lines))


def main():
    ab_specs = probe_antibody_chains()
    s5 = PARSER.get_structure("5O4Y_pair", str(PDB_DIR / "5O4Y.pdb"))
    o4y_specs = [
        ("5O4Y", p, [q], "macrocycle", f"peptide-71_{p}{q}")
        for p, q, _ in auto_pair_5o4y(s5)
    ]
    all_specs = COMPLEXES + o4y_specs + ab_specs

    by_class = defaultdict(Counter)
    presence = defaultdict(Counter)
    per_complex = {}
    aa_ref = {}
    macro_copy_contacts = []

    for pdb, pdl1_ch, partners, cls, label in all_specs:
        path = PDB_DIR / f"{pdb}.pdb"
        s = PARSER.get_structure(f"{pdb}_{label}", str(path))
        counts = contact_counts(s, pdl1_ch, partners)
        am = aa_map(s, pdl1_ch)
        aa_ref.update({k: v for k, v in am.items() if k not in aa_ref})
        key = f"{pdb}:{label}"
        per_complex[key] = {
            "class": cls,
            "label": label,
            "pdl1_chain": pdl1_ch,
            "partners": partners,
            "contacts": dict(sorted(counts.items())),
            "n_contact_residues": len(counts),
            "total_atom_pairs": int(sum(counts.values())),
        }
        if "peptide-71" in label:
            macro_copy_contacts.append(counts)
            continue
        for r, n in counts.items():
            by_class[cls][r] += n
            presence[cls][r] += 1

    if macro_copy_contacts:
        all_res = set().union(*[set(c) for c in macro_copy_contacts])
        collapsed = Counter()
        for r in all_res:
            vals = [c.get(r, 0) for c in macro_copy_contacts]
            if sum(1 for v in vals if v > 0) >= 2:
                collapsed[r] = int(round(np.mean(vals)))
        per_complex["5O4Y:peptide-71_consensus"] = {
            "class": "macrocycle",
            "label": "peptide-71_consensus",
            "contacts": dict(sorted(collapsed.items())),
            "n_contact_residues": len(collapsed),
            "total_atom_pairs": int(sum(collapsed.values())),
            "note": "auto-paired ASU; ≥2/3 copies; includes ncAA HETATM",
        }
        for r, n in collapsed.items():
            by_class["macrocycle"][r] += n
            presence["macrocycle"][r] += 1

    n_macro = 5
    n_ab = sum(1 for x in all_specs if x[3] == "antibody")
    residues = sorted(
        set(by_class["macrocycle"])
        | set(by_class["pd1"])
        | set(by_class["antibody"])
    )

    rows = []
    for r in residues:
        aa = aa_ref.get(r, "UNK")
        mac_pres = presence["macrocycle"][r]
        mac_pairs = by_class["macrocycle"][r]
        pd1_pairs = by_class["pd1"][r]
        ab_pres = presence["antibody"][r]
        ab_pairs = by_class["antibody"][r]
        consensus_macro = mac_pres / n_macro
        mac_density = mac_pairs / mac_pres if mac_pres else 0.0
        pd1_hit = 1 if pd1_pairs > 0 else 0
        core_score = (
            3.0 * consensus_macro
            + 2.0 * pd1_hit
            + 1.0 * min(mac_density / 20.0, 1.5)
            + 0.5 * (ab_pres / max(n_ab, 1))
        )
        rows.append(
            {
                "resnum": r,
                "aa": aa,
                "label": f"{aa}{r}",
                "macro_n_structures": mac_pres,
                "macro_fraction": round(consensus_macro, 3),
                "macro_atom_pairs_sum": mac_pairs,
                "macro_mean_pairs_when_hit": round(mac_density, 2),
                "pd1_atom_pairs": pd1_pairs,
                "pd1_contact": bool(pd1_pairs),
                "antibody_n_structures": ab_pres,
                "antibody_atom_pairs_sum": ab_pairs,
                "core_score": round(core_score, 3),
                "flag_pd1_gap": pd1_pairs > 0 and mac_pres <= 2,
                "flag_ab_heavy_macro_light": ab_pres >= max(n_ab - 1, 1)
                and mac_pres <= 2,
            }
        )
    rows.sort(
        key=lambda x: (
            -x["core_score"],
            -x["macro_n_structures"],
            -x["macro_atom_pairs_sum"],
        )
    )

    pd1_dense = 15
    for row in rows:
        mac = row["macro_n_structures"]
        pd1 = row["pd1_contact"]
        pd1_n = row["pd1_atom_pairs"]
        if mac >= 4 and pd1:
            tier = "core_required"
        elif mac >= 3 and pd1 and pd1_n >= pd1_dense:
            tier = "core_required"
        elif mac >= 3 and pd1:
            tier = "core_strong"
        elif mac >= 3 or (pd1 and mac >= 2):
            tier = "extended"
        elif pd1 and mac <= 1 and pd1_n >= 5:
            tier = "pd1_gap_opportunity"
        elif row["flag_ab_heavy_macro_light"] and pd1:
            tier = "ab_rim_opportunity"
        else:
            tier = "peripheral"
        row["tier"] = tier

    csv_path = OUT_DIR / "hotspot_consensus.csv"
    cols = [
        "resnum",
        "aa",
        "label",
        "tier",
        "core_score",
        "macro_n_structures",
        "macro_fraction",
        "macro_atom_pairs_sum",
        "macro_mean_pairs_when_hit",
        "pd1_atom_pairs",
        "pd1_contact",
        "antibody_n_structures",
        "antibody_atom_pairs_sum",
        "flag_pd1_gap",
        "flag_ab_heavy_macro_light",
    ]
    with csv_path.open("w") as f:
        f.write(",".join(cols) + "\n")
        for row in rows:
            f.write(",".join(str(row[c]) for c in cols) + "\n")

    core = [r for r in rows if r["tier"] == "core_required"]
    strong = [r for r in rows if r["tier"] == "core_strong"]
    gaps = [
        r
        for r in rows
        if r["tier"] in ("pd1_gap_opportunity", "ab_rim_opportunity")
    ]

    # Curated 6-residue RF set
    priority = [56, 123, 115, 113, 66, 121, 54, 58]
    pool = {r["resnum"]: r for r in core + strong}
    rf_picked = []
    for num in priority:
        if num in pool and len(rf_picked) < 6:
            rf_picked.append(pool[num])
    for r in core:
        if r["resnum"] not in {x["resnum"] for x in rf_picked} and len(rf_picked) < 6:
            rf_picked.append(r)
    recommended_string = [f"A{r['resnum']}" for r in rf_picked]

    ranked = {
        "accessed": str(date.today()),
        "method": {
            "contact_cutoff_A": CUTOFF,
            "macrocycle_structures": n_macro,
            "antibody_structures": n_ab,
            "pd1_structure": "4ZQK",
            "ligand_atoms": "polymer + non-water HETATM (ncAA)",
            "core_score": "3*macro_frac + 2*pd1_hit + min(density/20,1.5) + 0.5*ab_frac",
            "tier_override": "macro>=3 & PD-1 pairs>=15 -> core_required",
        },
        "rfdiffusion_recommended": {
            "tier1_string": recommended_string,
            "rationale": "consensus ∩ PD-1; Tyr123 by density; capped at 6",
        },
        "tiers": {
            "core_required": [
                {
                    "res": r["resnum"],
                    "aa": r["aa"],
                    "score": r["core_score"],
                    "macro_n": r["macro_n_structures"],
                    "pd1_pairs": r["pd1_atom_pairs"],
                }
                for r in rows
                if r["tier"] == "core_required"
            ],
            "core_strong": [
                {
                    "res": r["resnum"],
                    "aa": r["aa"],
                    "score": r["core_score"],
                    "macro_n": r["macro_n_structures"],
                    "pd1_pairs": r["pd1_atom_pairs"],
                }
                for r in rows
                if r["tier"] == "core_strong"
            ],
            "differentiation_candidates": [
                {
                    "res": r["resnum"],
                    "aa": r["aa"],
                    "tier": r["tier"],
                    "macro_n": r["macro_n_structures"],
                    "pd1_pairs": r["pd1_atom_pairs"],
                    "ab_n": r["antibody_n_structures"],
                }
                for r in gaps
            ],
        },
        "top20_by_core_score": rows[:20],
        "per_complex_summary": {
            k: {
                "class": v["class"],
                "n_contact_residues": v["n_contact_residues"],
                "total_atom_pairs": v["total_atom_pairs"],
                "top_contacts": sorted(v["contacts"].items(), key=lambda x: -x[1])[
                    :10
                ],
            }
            for k, v in per_complex.items()
            if "peptide-71_" not in k or k.endswith("consensus")
        },
    }
    (OUT_DIR / "hotspot_consensus.json").write_text(json.dumps(ranked, indent=2))

    hotspots_ranked = {
        "target": "PD-L1",
        "uniprot": "Q9NZQ7",
        "numbering_reference": "UniProt_IgV_as_in_7OUN_chain_A",
        "accessed": str(date.today()),
        "source_analysis": "results/PDL1/hotspot_consensus.json",
        "core_hotspots_rfdiffusion": [
            {
                "res": r["resnum"],
                "aa": r["aa"],
                "rf": f"A{r['resnum']}",
                "macro_n": r["macro_n_structures"],
                "tier": r["tier"],
                "pd1_pairs": r["pd1_atom_pairs"],
            }
            for r in rf_picked
        ],
        "differentiation_growth_hotspots": [
            {
                "res": r["resnum"],
                "aa": r["aa"],
                "reason": r["tier"],
                "pd1_pairs": r["pd1_atom_pairs"],
                "macro_n": r["macro_n_structures"],
                "ab_n": r["antibody_n_structures"],
            }
            for r in gaps
        ],
        "recommended_rf_string_tier1": recommended_string,
        "blockade_filter_residues": [
            {"res": 122, "aa": "ASP", "why": "PD-1 dense; polar rim"},
            {"res": 124, "aa": "LYS", "why": "PD-1 2nd densest; macro-light"},
            {
                "res": 125,
                "aa": "ARG",
                "why": "PD-1 + Abs; BMS salt-bridge — redesign chemically",
            },
        ],
    }
    (ROOT / "data/targets/PDL1/hotspots_ranked.json").write_text(
        json.dumps(hotspots_ranked, indent=2)
    )
    write_diff_report(rows, hotspots_ranked, n_macro, n_ab)

    print("\n=== TOP 15 by core_score ===")
    for r in rows[:15]:
        print(
            f"{r['label']:8} tier={r['tier']:22} score={r['core_score']:5.2f} "
            f"macro={r['macro_n_structures']}/{n_macro} dens={r['macro_mean_pairs_when_hit']:5.1f} "
            f"PD1={r['pd1_atom_pairs']:3} Ab={r['antibody_n_structures']}/{n_ab}"
        )
    print("\nRF tier1:", recommended_string)
    print(f"Wrote {csv_path}")
    print(f"Wrote {OUT_DIR / 'hotspot_consensus.json'}")
    print(f"Wrote {OUT_DIR / 'differentiation_opportunities.md'}")
    print(f"Wrote data/targets/PDL1/hotspots_ranked.json")


if __name__ == "__main__":
    main()
