#!/usr/bin/env python3
"""PD-L1 blockade / diversity metrics for designed or crystal complexes (S1b).

Metrics (per complex):
  - pd1_jaccard: peptide-contact residues ∩ 4ZQK PD-1 footprint / union
  - rim_D122 / rim_K124 / rim_R125: heavy-atom contact pairs ≤ cutoff
  - rim_score: soft weighted sum of per-residue normalized rim contacts (NOT AND)
  - rmsd_vs_8alx: peptide CA RMSD vs 8ALX peptide after receptor-interface align
    (only when peptide length matches reference; else null)

Usage:
  python scripts/pdl1_blockade_metrics.py
  python scripts/pdl1_blockade_metrics.py --complexes-json path.json
"""
from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parents[1]
PDB_DIR = ROOT / "data/targets/PDL1/pdb"
SPEC = ROOT / "data/targets/PDL1/design_spec.yaml"
OUT_DIR = ROOT / "results/PDL1/filter_calibration"
CUTOFF = 4.5
WATER = {"HOH", "WAT", "H2O", "DOD"}
PARSER = PDBParser(QUIET=True)

RIM_RESIDUES = (122, 124, 125)
RIM_NAMES = {122: "D122", 124: "K124", 125: "R125"}
# Soft weights: K124 gap emphasized for differentiation ranking
RIM_WEIGHTS = {122: 1.0, 124: 1.2, 125: 1.0}
RIM_NORM_CAP = 20.0  # atom-pair cap for soft saturation


def load_footprint() -> set[int]:
    spec = yaml.safe_load(SPEC.read_text())
    return set(spec["hotspots"]["pd1_footprint_4zqk_pdl1"])


def polymer_residues(chain):
    return [r for r in chain if r.id[0] == " "]


def ligand_residues(chain):
    return [r for r in chain if r.resname not in WATER]


def heavy_atoms(residue):
    return [a for a in residue if a.element != "H"]


def get_res(chain, resid: int):
    for r in polymer_residues(chain):
        if r.id[1] == resid:
            return r
    return None


def partner_coords(model, partner_chain_ids: list[str]) -> np.ndarray:
    atoms = []
    for cid in partner_chain_ids:
        if cid not in model:
            continue
        for r in ligand_residues(model[cid]):
            atoms.extend(heavy_atoms(r))
    if not atoms:
        return np.zeros((0, 3))
    return np.array([a.coord for a in atoms], dtype=float)


def contact_residue_set(
    model, pdl1_chain: str, partner_ids: list[str], cutoff: float = CUTOFF
) -> set[int]:
    pc = partner_coords(model, partner_ids)
    if len(pc) == 0:
        return set()
    hit = set()
    for res in polymer_residues(model[pdl1_chain]):
        atoms = heavy_atoms(res)
        if not atoms:
            continue
        coords = np.array([a.coord for a in atoms], dtype=float)
        d2 = ((coords[:, None, :] - pc[None, :, :]) ** 2).sum(axis=2)
        if (d2 <= cutoff**2).any():
            hit.add(res.id[1])
    return hit


def contact_pairs_to_residue(
    model, pdl1_chain: str, partner_ids: list[str], resid: int, cutoff: float = CUTOFF
) -> int:
    res = get_res(model[pdl1_chain], resid)
    if res is None:
        return 0
    atoms = heavy_atoms(res)
    pc = partner_coords(model, partner_ids)
    if not atoms or len(pc) == 0:
        return 0
    coords = np.array([a.coord for a in atoms], dtype=float)
    d2 = ((coords[:, None, :] - pc[None, :, :]) ** 2).sum(axis=2)
    return int((d2 <= cutoff**2).sum())


def kabsch_rmsd(A: np.ndarray, B: np.ndarray) -> float:
    """A, B: (N,3). Return RMSD after optimal rotation of B onto A."""
    if len(A) < 3 or len(A) != len(B):
        return float("nan")
    Ac, Bc = A.mean(0), B.mean(0)
    A0, B0 = A - Ac, B - Bc
    H = A0.T @ B0
    U, _S, Vt = np.linalg.svd(H)
    d = np.linalg.det(U @ Vt)
    D = np.diag([1.0, 1.0, 1.0 if d >= 0 else -1.0])
    R = U @ D @ Vt
    B2 = (R @ B0.T).T + Ac
    return float(np.sqrt(((A - B2) ** 2).sum() / len(A)))


def ca_coords(chain, resids: list[int]) -> np.ndarray | None:
    coords = []
    for rid in resids:
        r = get_res(chain, rid)
        if r is None or "CA" not in r:
            return None
        coords.append(r["CA"].coord.copy())
    return np.array(coords, dtype=float)


def peptide_ca_list(chain) -> list[tuple[int, np.ndarray]]:
    out = []
    for r in polymer_residues(chain):
        if "CA" in r:
            out.append((r.id[1], r["CA"].coord.copy()))
        else:
            # ncAA may lack CA; try any backbone-ish
            for name in ("CA", "C", "N"):
                if name in r:
                    out.append((r.id[1], r[name].coord.copy()))
                    break
    return out


def rmsd_vs_reference_peptide(
    model,
    pdl1_chain: str,
    pep_chain_ids: list[str],
    ref_model,
    ref_pdl1: str,
    ref_pep_ids: list[str],
    iface_resids: list[int],
) -> float | None:
    """Align receptors on interface CA, then RMSD peptide CA (equal-length only)."""
    A = ca_coords(ref_model[ref_pdl1], iface_resids)
    B = ca_coords(model[pdl1_chain], iface_resids)
    if A is None or B is None:
        return None
    Ac, Bc = A.mean(0), B.mean(0)
    A0, B0 = A - Ac, B - Bc
    H = A0.T @ B0
    U, _S, Vt = np.linalg.svd(H)
    d = np.linalg.det(U @ Vt)
    D = np.diag([1.0, 1.0, 1.0 if d >= 0 else -1.0])
    R = U @ D @ Vt

    def transform(coords: np.ndarray) -> np.ndarray:
        return (R @ (coords - Bc).T).T + Ac

    # collect polymer peptide CAs from first partner chain with polymer residues
    def pep_cas(m, cids):
        for cid in cids:
            if cid not in m:
                continue
            polymer = polymer_residues(m[cid])
            if polymer:
                return peptide_ca_list(m[cid])
            # fall back: all ligand residues with CA
            pts = []
            for r in ligand_residues(m[cid]):
                if "CA" in r:
                    pts.append((r.id[1], r["CA"].coord.copy()))
            if pts:
                return pts
        return []

    ref_pep = pep_cas(ref_model, ref_pep_ids)
    mov_pep = pep_cas(model, pep_chain_ids)
    if len(ref_pep) < 3 or len(mov_pep) < 3:
        return None
    # equal-length pairwise by order (crystal macros differ in length → NaN)
    if len(ref_pep) != len(mov_pep):
        return None
    Rref = np.array([c for _, c in ref_pep], dtype=float)
    Rmov = transform(np.array([c for _, c in mov_pep], dtype=float))
    return kabsch_rmsd(Rref, Rmov)


def rim_score_from_pairs(pairs: dict[int, int]) -> float:
    parts = []
    wsum = 0.0
    for rid, w in RIM_WEIGHTS.items():
        n = pairs.get(rid, 0)
        sat = min(n / RIM_NORM_CAP, 1.0)
        parts.append(w * sat)
        wsum += w
    return float(sum(parts) / wsum) if wsum else 0.0


def score_complex(
    pdb_path: Path,
    pdl1_chain: str,
    partner_ids: list[str],
    footprint: set[int],
    tag: str,
    role: str,
    ref_8alx=None,
    iface_for_align: list[int] | None = None,
) -> dict:
    s = PARSER.get_structure(tag, str(pdb_path))
    model = s[0]
    contacts = contact_residue_set(model, pdl1_chain, partner_ids)
    inter = contacts & footprint
    union = contacts | footprint
    jaccard = float(len(inter) / len(union)) if union else 0.0
    rim_pairs = {
        rid: contact_pairs_to_residue(model, pdl1_chain, partner_ids, rid)
        for rid in RIM_RESIDUES
    }
    row = {
        "tag": tag,
        "role": role,
        "pdb": str(pdb_path.relative_to(ROOT)) if pdb_path.is_relative_to(ROOT) else str(pdb_path),
        "pdl1_chain": pdl1_chain,
        "partner_chains": ",".join(partner_ids),
        "n_contact_res": len(contacts),
        "n_footprint_overlap": len(inter),
        "pd1_jaccard": round(jaccard, 4),
        "rim_D122": rim_pairs[122],
        "rim_K124": rim_pairs[124],
        "rim_R125": rim_pairs[125],
        "rim_score": round(rim_score_from_pairs(rim_pairs), 4),
        "rim_all_three": int(all(rim_pairs[r] > 0 for r in RIM_RESIDUES)),
        "rmsd_vs_8alx": None,
    }
    if ref_8alx is not None and iface_for_align is not None:
        ref_s, ref_p, ref_pep = ref_8alx
        rms = rmsd_vs_reference_peptide(
            model,
            pdl1_chain,
            partner_ids,
            ref_s[0],
            ref_p,
            ref_pep,
            iface_for_align,
        )
        row["rmsd_vs_8alx"] = None if rms is None or np.isnan(rms) else round(float(rms), 3)
    return row


DEFAULT_COMPLEXES = [
    # positives (crystal macros)
    {
        "tag": "7OUN_peptide104",
        "role": "positive",
        "pdb": "7OUN.pdb",
        "pdl1": "A",
        "partners": ["B"],
    },
    {
        "tag": "5O45_peptide57",
        "role": "positive",
        "pdb": "5O45.pdb",
        "pdl1": "A",
        "partners": ["B"],
    },
    {
        "tag": "8ALX_pAC65",
        "role": "positive",
        "pdb": "8ALX.pdb",
        "pdl1": "A",
        "partners": ["B"],
    },
    {
        "tag": "6PV9_macro",
        "role": "positive",
        "pdb": "6PV9.pdb",
        "pdl1": "A",
        "partners": ["B"],
    },
    {
        "tag": "5O4Y_peptide71_BD",
        "role": "positive",
        "pdb": "5O4Y.pdb",
        "pdl1": "B",
        "partners": ["D"],
    },
    # PD-1 itself — footprint definition reference (not a macro positive)
    {
        "tag": "4ZQK_PD1",
        "role": "footprint_ref",
        "pdb": "4ZQK.pdb",
        "pdl1": "A",
        "partners": ["B"],
    },
]


def build_negatives_from_positives(footprint: set[int]) -> list[dict]:
    """Geometry negatives: same PD-L1, but score as if peptide only hits non-hotspot face.

    Implemented as synthetic rows by computing contacts then zeroing overlap with
    core hotspots / footprint for a 'hotspot_miss' proxy — see calibrate script
    for true random decoys. Here we add apo-like null: empty partner.
    """
    return [
        {
            "tag": "5O45_no_ligand",
            "role": "negative_null",
            "pdb": "PDL1_IgV_5O45_A18-134.pdb",
            "pdl1": "A",
            "partners": [],  # no ligand → jaccard 0
        }
    ]


def run(complexes: list[dict] | None = None) -> pd.DataFrame:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    footprint = load_footprint()
    iface = sorted(footprint)
    # 8ALX reference for RMSD
    ref_path = PDB_DIR / "8ALX.pdb"
    ref_s = PARSER.get_structure("8ALX", str(ref_path))
    ref_8alx = (ref_s, "A", ["B"])

    rows = []
    specs = complexes if complexes is not None else DEFAULT_COMPLEXES + build_negatives_from_positives(footprint)
    for c in specs:
        pdb_path = PDB_DIR / c["pdb"] if not Path(c["pdb"]).is_absolute() else Path(c["pdb"])
        if not pdb_path.exists():
            # allow relative to ROOT
            alt = ROOT / c["pdb"]
            pdb_path = alt if alt.exists() else pdb_path
        if not pdb_path.exists():
            rows.append(
                {
                    "tag": c["tag"],
                    "role": c["role"],
                    "pdb": c["pdb"],
                    "error": "missing_pdb",
                }
            )
            continue
        partners = c.get("partners") or []
        row = score_complex(
            pdb_path,
            c["pdl1"],
            partners,
            footprint,
            c["tag"],
            c["role"],
            ref_8alx=ref_8alx if partners else None,
            iface_for_align=iface if partners else None,
        )
        rows.append(row)
    df = pd.DataFrame(rows)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complexes-json", type=Path, default=None)
    ap.add_argument(
        "--out-csv",
        type=Path,
        default=OUT_DIR / "crystal_blockade_metrics.csv",
    )
    ap.add_argument(
        "--out-json",
        type=Path,
        default=OUT_DIR / "crystal_blockade_metrics.json",
    )
    args = ap.parse_args()
    complexes = None
    if args.complexes_json:
        complexes = json.loads(args.complexes_json.read_text())
    df = run(complexes)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_csv, index=False)
    payload = {
        "accessed": date.today().isoformat(),
        "cutoff_A": CUTOFF,
        "rim_mode": "soft_score_plus_quota",
        "rim_weights": RIM_WEIGHTS,
        "note": "rim_all_three is diagnostic only; NEVER use as default gate",
        "rows": df.replace({np.nan: None}).to_dict(orient="records"),
    }
    args.out_json.write_text(json.dumps(payload, indent=2) + "\n")
    print(df.to_string(index=False))
    print(f"\nWrote {args.out_csv}")
    print(f"Wrote {args.out_json}")


if __name__ == "__main__":
    main()
