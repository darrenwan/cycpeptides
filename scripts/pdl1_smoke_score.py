#!/usr/bin/env python3
"""Score S2 smoke complexes with blockade metrics + thresholds_v1 (uv run).

Usage:
  uv run python scripts/pdl1_smoke_score.py
"""
from __future__ import annotations

import json
import sys
from datetime import date
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from pdl1_blockade_metrics import score_complex, load_footprint, PARSER  # noqa: E402

SMOKE = ROOT / "results/PDL1/smoke"
THRESH = ROOT / "data/targets/PDL1/thresholds_v1.yaml"
OUT_CSV = SMOKE / "smoke_metrics.csv"
OUT_JSON = SMOKE / "smoke_report.json"


def detect_chains(pdb: Path) -> tuple[str, list[str]]:
    s = PARSER.get_structure(pdb.stem, str(pdb))
    chains = [c.id for c in s[0]]
    # RFpeptides binder design: binder=A, target=B typically
    if "A" in chains and "B" in chains:
        return "B", ["A"]  # PD-L1 is B, peptide A
    if len(chains) >= 2:
        return chains[1], [chains[0]]
    return chains[0], []


def main() -> None:
    thresh = yaml.safe_load(THRESH.read_text())
    t_j = float(thresh["primary_gate"]["T_pd1_jaccard"])
    footprint = load_footprint()
    iface = sorted(footprint)

    afcyc = list((SMOKE / "afcyc").glob("*_afcyc_pred.pdb"))
    designed = list((SMOKE / "mpnn/designed_pdbs").glob("*.pdb"))
    backbones = list((SMOKE / "rfpeptides").glob("pdl1_cyc_*.pdb"))
    if afcyc:
        pool, role = afcyc, "smoke_afcyc"
    elif designed:
        pool, role = designed, "smoke_designed"
    else:
        pool, role = backbones, "smoke_backbone"

    if not pool:
        raise SystemExit(f"No smoke PDBs under {SMOKE}")
    if role != "smoke_afcyc":
        print(
            f"WARN: scoring {role} — AfCyc preds missing; S2 incomplete without S2c",
            flush=True,
        )

    ref_path = ROOT / "data/targets/PDL1/pdb/8ALX.pdb"
    ref_s = PARSER.get_structure("8ALX", str(ref_path))
    ref_8alx = (ref_s, "A", ["B"])

    rows = []
    for pdb in sorted(pool):
        pdl1, partners = detect_chains(pdb)
        row = score_complex(
            pdb,
            pdl1,
            partners,
            footprint,
            pdb.stem,
            role,
            ref_8alx=ref_8alx if partners else None,
            iface_for_align=iface if partners else None,
        )
        row["pass_P0"] = bool(row.get("pd1_jaccard", 0) >= t_j)
        rows.append(row)

    df = pd.DataFrame(rows)
    SMOKE.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    report = {
        "accessed": date.today().isoformat(),
        "n_scored": len(df),
        "source": role,
        "requires_afcyc_for_s2_done": True,
        "afcyc_present": role == "smoke_afcyc",
        "T_pd1_jaccard": t_j,
        "n_pass_P0": int(df["pass_P0"].sum()) if len(df) else 0,
        "metrics_csv": str(OUT_CSV.relative_to(ROOT)),
        "python": sys.executable,
        "note": "S2 score via uv; prefer AfCyc preds when present",
    }
    OUT_JSON.write_text(json.dumps(report, indent=2) + "\n")
    cols = [c for c in ("tag", "pd1_jaccard", "rim_score", "pass_P0") if c in df.columns]
    print(df[cols].to_string(index=False))
    print(json.dumps(report, indent=2))
    if len(df) < 1:
        raise SystemExit("Smoke score FAILED: empty metrics")
    if role != "smoke_afcyc":
        raise SystemExit("Smoke score WARN→FAIL: AfCyc predictions required for S2 complete")
    print("S2 smoke score PASS (AfCyc preds)")


if __name__ == "__main__":
    main()
