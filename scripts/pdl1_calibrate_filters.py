#!/usr/bin/env python3
"""Calibrate PD-L1 dry-lab filter thresholds from crystal controls (S1c).

- Positives: crystal macrocycle complexes (must pass P0 footprint gate)
- Negatives: null ligand + translated-peptide decoys (must mostly fail P0)
- rim_all_three is reported but NEVER used as a default hard gate

Writes:
  data/targets/PDL1/thresholds_v1.yaml
  results/PDL1/filter_calibration/calibration_report.md
  results/PDL1/filter_calibration/calibration_summary.json
  results/PDL1/filter_calibration/decoy_*.pdb (几何阴性)

Usage:
  python scripts/pdl1_calibrate_filters.py
"""
from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from Bio.PDB import PDBIO, PDBParser, Structure, Model, Chain

ROOT = Path(__file__).resolve().parents[1]
PDB_DIR = ROOT / "data/targets/PDL1/pdb"
OUT_DIR = ROOT / "results/PDL1/filter_calibration"
THRESH_PATH = ROOT / "data/targets/PDL1/thresholds_v1.yaml"

import sys

sys.path.insert(0, str(ROOT / "scripts"))
from pdl1_blockade_metrics import (  # noqa: E402
    DEFAULT_COMPLEXES,
    CUTOFF,
    RIM_WEIGHTS,
    run as run_metrics,
    score_complex,
    load_footprint,
    polymer_residues,
    ligand_residues,
)

PARSER = PDBParser(QUIET=True)


def write_translated_decoy(
    src_pdb: Path,
    pdl1_chain: str,
    pep_chain: str,
    out_pdb: Path,
    shift: np.ndarray,
) -> None:
    """Copy complex; translate peptide chain by `shift` (Å)."""
    s = PARSER.get_structure("src", str(src_pdb))
    new = Structure.Structure(out_pdb.stem)
    model = Model.Model(0)
    new.add(model)
    for chain in s[0]:
        ch = Chain.Chain(chain.id)
        model.add(ch)
        for res in chain:
            if res.id[0] not in (" ", "H_"):
                # keep polymer + HETATM ligand
                if res.resname in {"HOH", "WAT", "H2O", "DOD"}:
                    continue
            rc = res.copy()
            if chain.id == pep_chain:
                for atom in rc:
                    atom.coord = atom.coord + shift
            ch.add(rc)
    # Bio.PDB residue id filter: include polymer and non-water HETATM
    # Rebuild more carefully
    new = Structure.Structure(out_pdb.stem)
    model = Model.Model(0)
    new.add(model)
    for chain in s[0]:
        ch = Chain.Chain(chain.id)
        model.add(ch)
        for res in chain:
            if res.resname in {"HOH", "WAT", "H2O", "DOD"}:
                continue
            if res.id[0] == " " or res.id[0].startswith("H_"):
                rc = res.copy()
                if chain.id == pep_chain:
                    for atom in rc:
                        atom.coord = np.asarray(atom.coord, dtype=float) + shift
                ch.add(rc)
    io = PDBIO()
    io.set_structure(new)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io.save(str(out_pdb))


def build_decoy_specs() -> list[dict]:
    """Far-translated peptide decoys from 7OUN / 5O45."""
    specs = []
    plans = [
        ("7OUN.pdb", "A", "B", "decoy_7OUN_pep_plus50A.pdb", np.array([50.0, 0.0, 0.0])),
        ("5O45.pdb", "A", "B", "decoy_5O45_pep_plus50A.pdb", np.array([0.0, 50.0, 0.0])),
        ("8ALX.pdb", "A", "B", "decoy_8ALX_pep_plus50A.pdb", np.array([0.0, 0.0, 50.0])),
    ]
    for src, p, pep, name, shift in plans:
        out = OUT_DIR / name
        write_translated_decoy(PDB_DIR / src, p, pep, out, shift)
        specs.append(
            {
                "tag": name.replace(".pdb", ""),
                "role": "negative_decoy",
                "pdb": str(out.relative_to(ROOT)),
                "pdl1": p,
                "partners": [pep],
            }
        )
    # null
    specs.append(
        {
            "tag": "5O45_no_ligand",
            "role": "negative_null",
            "pdb": "PDL1_IgV_5O45_A18-134.pdb",
            "pdl1": "A",
            "partners": [],
        }
    )
    return specs


def choose_threshold(pos_j: np.ndarray, neg_j: np.ndarray) -> dict:
    """Pick T_jaccard so all positives pass and most negatives fail.

    Default: T = min(pos) - small epsilon, floored at 0.15.
    If that fails to separate, use midpoint between min(pos) and max(neg).
    """
    pos_j = np.asarray(pos_j, dtype=float)
    neg_j = np.asarray(neg_j, dtype=float)
    min_pos = float(pos_j.min())
    max_neg = float(neg_j.max()) if len(neg_j) else 0.0
    # Prefer just below min positive, but above max negative when separable
    if min_pos > max_neg:
        t = max(0.15, (min_pos + max_neg) / 2.0)
        # ensure all positives still pass
        t = min(t, min_pos - 1e-4)
        t = max(t, 0.0)
    else:
        # overlapping — set to min_pos * 0.9 and flag
        t = max(0.0, min_pos * 0.9)
    pos_pass = float((pos_j >= t).mean())
    neg_fail = float((neg_j < t).mean()) if len(neg_j) else 1.0
    return {
        "T_pd1_jaccard": round(t, 4),
        "min_positive_jaccard": round(min_pos, 4),
        "max_negative_jaccard": round(max_neg, 4),
        "positive_pass_rate": round(pos_pass, 4),
        "negative_fail_rate": round(neg_fail, 4),
        "separable": bool(min_pos > max_neg),
    }


def df_to_md(df: pd.DataFrame) -> str:
    if df is None or len(df) == 0:
        return "_无_"
    cols = list(df.columns)
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join("---" for _ in cols) + " |",
    ]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
    return "\n".join(lines)


def df_to_md(df: pd.DataFrame) -> str:
    if df is None or len(df) == 0:
        return "_无_"
    cols = list(df.columns)
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join("---" for _ in cols) + " |",
    ]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
    return "\n".join(lines)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # Ensure cropped receptor exists for null negative
    cropped = PDB_DIR / "PDL1_IgV_5O45_A18-134.pdb"
    if not cropped.exists():
        raise SystemExit("Run scripts/pdl1_prepare_receptor.py first")

    decoys = build_decoy_specs()
    positives = [c for c in DEFAULT_COMPLEXES if c["role"] == "positive"]
    footprint_ref = [c for c in DEFAULT_COMPLEXES if c["role"] == "footprint_ref"]
    all_specs = positives + footprint_ref + decoys

    # Write complexes list for metrics
    complexes_json = OUT_DIR / "calibration_complexes.json"
    complexes_json.write_text(json.dumps(all_specs, indent=2) + "\n")

    df = run_metrics(all_specs)
    metrics_csv = OUT_DIR / "calibration_metrics.csv"
    df.to_csv(metrics_csv, index=False)

    pos = df[df["role"] == "positive"].copy()
    neg = df[df["role"].astype(str).str.startswith("negative")].copy()

    cal = choose_threshold(
        pos["pd1_jaccard"].astype(float).values,
        neg["pd1_jaccard"].astype(float).values,
    )

    # Rim diagnostics: show AND would kill positives
    and_kill = int((pos["rim_all_three"] == 0).sum()) if "rim_all_three" in pos else -1

    thresholds = {
        "version": "v1",
        "accessed": date.today().isoformat(),
        "contact_cutoff_A": CUTOFF,
        "primary_gate": {
            "metric": "pd1_jaccard",
            "T_pd1_jaccard": cal["T_pd1_jaccard"],
            "rule": "pass if pd1_jaccard >= T_pd1_jaccard",
        },
        "rim": {
            "mode": "soft_score_plus_quota",
            "residues": ["A122", "A124", "A125"],
            "weights": {f"A{k}": v for k, v in RIM_WEIGHTS.items()},
            "forbid_require_all_three": True,
            "note": "Use rim_score for ranking/quota only",
        },
        "diversity": {
            "rmsd_vs_8alx_bins_A": [0.0, 2.0, 4.0, 6.0, 999.0],
            "high_rmsd_quota_frac": 0.30,
            "high_rim_quota_frac": 0.30,
            "top_composite_frac": 0.40,
        },
        "structure": {
            "peptide_ca_rmsd_max_A": 2.0,
            "note": "Applied after structure prediction (S2+); not from crystal calib",
        },
        "calibration_stats": cal,
        "diagnostics": {
            "n_positives": int(len(pos)),
            "n_negatives": int(len(neg)),
            "positives_failing_rim_all_three": and_kill,
            "positives_rim_all_three_would_kill_frac": round(and_kill / max(len(pos), 1), 4),
        },
        "provenance": {
            "metrics_csv": str(metrics_csv.relative_to(ROOT)),
            "complexes_json": str(complexes_json.relative_to(ROOT)),
            "script": "scripts/pdl1_calibrate_filters.py",
        },
    }
    THRESH_PATH.write_text(yaml.safe_dump(thresholds, sort_keys=False, allow_unicode=True))

    # Markdown report
    lines = [
        "# PD-L1 过滤阈值校准报告（S1c）",
        "",
        f"> 日期：{date.today().isoformat()}  ·  接触阈值：{CUTOFF} Å  ·  阈值文件：`data/targets/PDL1/thresholds_v1.yaml`",
        "",
        "## 结论",
        "",
        f"- **主门控**：`pd1_jaccard >= {cal['T_pd1_jaccard']}`",
        f"- 阳性通过率：{cal['positive_pass_rate']:.0%}（min jaccard={cal['min_positive_jaccard']}）",
        f"- 阴性落选率：{cal['negative_fail_rate']:.0%}（max jaccard={cal['max_negative_jaccard']}）",
        f"- 可分性：{'是' if cal['separable'] else '否（需人工复核）'}",
        f"- **rim AND 诊断**：{and_kill}/{len(pos)} 个阳性宏环 **不过** `D122∧K124∧R125` → 禁止作默认硬门",
        "",
        "## 阳性宏环",
        "",
        df_to_md(
            pos[
                [
                    "tag",
                    "pd1_jaccard",
                    "rim_D122",
                    "rim_K124",
                    "rim_R125",
                    "rim_score",
                    "rim_all_three",
                ]
            ]
        ),
        "",
        "## 阴性对照",
        "",
        df_to_md(
            neg[
                [
                    "tag",
                    "pd1_jaccard",
                    "rim_D122",
                    "rim_K124",
                    "rim_R125",
                    "rim_score",
                ]
            ]
        )
        if len(neg)
        else "_无_",
        "",
        "## 门控策略（冻结 v1）",
        "",
        "1. P0：足迹 Jaccard ≥ T",
        "2. rim：分项软分 + 配额（非 AND）",
        "3. 结构/物理阈值待 S2 冒烟后补校准",
        "",
        f"原始表：`{metrics_csv.relative_to(ROOT)}`",
        "",
    ]
    report_path = OUT_DIR / "calibration_report.md"
    report_path.write_text("\n".join(lines))

    summary = {
        "accessed": date.today().isoformat(),
        "thresholds": thresholds,
        "positive_tags": pos["tag"].tolist(),
        "negative_tags": neg["tag"].tolist(),
    }
    (OUT_DIR / "calibration_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    print(yaml.safe_dump(thresholds["primary_gate"], sort_keys=False))
    print(f"Wrote {THRESH_PATH}")
    print(f"Wrote {report_path}")
    if cal["positive_pass_rate"] < 1.0:
        raise SystemExit("Calibration FAILED: not all positives pass P0")
    if cal["negative_fail_rate"] < 0.8:
        raise SystemExit("Calibration FAILED: negatives not mostly rejected")
    print("S1c calibration PASS")


if __name__ == "__main__":
    main()
