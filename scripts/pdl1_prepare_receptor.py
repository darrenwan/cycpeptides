#!/usr/bin/env python3
"""Crop / validate PD-L1 design receptors for RFpeptides (S1a).

Writes:
  data/targets/PDL1/pdb/PDL1_IgV_5O45_A18-134.pdb
  results/PDL1/receptor_prep_report.json

Usage:
  python scripts/pdl1_prepare_receptor.py
"""
from __future__ import annotations

import json
from datetime import date
from pathlib import Path

from Bio.PDB import PDBIO, PDBParser, Structure, Model, Chain

ROOT = Path(__file__).resolve().parents[1]
PDB_DIR = ROOT / "data/targets/PDL1/pdb"
OUT_DIR = ROOT / "results/PDL1"
HOTSPOTS_JSON = ROOT / "data/targets/PDL1/hotspots.json"

# Contig used by Track A RFpeptides
RES_START, RES_END = 18, 134
SRC = PDB_DIR / "PDL1_IgV_5O45_chainA.pdb"
OUT_PDB = PDB_DIR / "PDL1_IgV_5O45_A18-134.pdb"

EXPECTED_HOTSPOT_AA = {
    56: "TYR",
    66: "GLN",
    113: "ARG",
    115: "MET",
    121: "ALA",
    123: "TYR",
    122: "ASP",
    124: "LYS",
    125: "ARG",
}


def polymer_residues(chain):
    return [r for r in chain if r.id[0] == " "]


def crop_chain(src_chain, out_chain_id: str, start: int, end: int):
    new = Structure.Structure("cropped")
    model = Model.Model(0)
    chain = Chain.Chain(out_chain_id)
    model.add(chain)
    new.add(model)
    kept = []
    for res in polymer_residues(src_chain):
        rid = res.id[1]
        if start <= rid <= end:
            chain.add(res.copy())
            kept.append(rid)
    return new, kept


def write_pdb(structure, out_path: Path, remarks: list[str]) -> None:
    io = PDBIO()
    io.set_structure(structure)
    tmp = out_path.with_suffix(".tmp.pdb")
    io.save(str(tmp))
    header = (
        "HEADER    PD-L1 IGV CROPPED FOR RFPEPTIDES\n"
        "TITLE     PD-L1 IgV 5O45 chain A cropped to design contig\n"
    )
    for r in remarks:
        header += f"REMARK    {r}\n"
    body = "\n".join(
        line for line in tmp.read_text().splitlines() if not line.startswith("END")
    )
    out_path.write_text(header + body + "\nEND\n")
    tmp.unlink(missing_ok=True)


def validate(kept_ids: list[int], structure) -> dict:
    present = set(kept_ids)
    aa_map = {
        r.id[1]: r.resname
        for r in polymer_residues(structure[0]["A"])
    }
    hotspot_ok = {}
    missing = []
    aa_mismatch = []
    for rid, aa in EXPECTED_HOTSPOT_AA.items():
        if rid not in present:
            missing.append(rid)
            hotspot_ok[rid] = {"present": False}
            continue
        got = aa_map.get(rid)
        ok = got == aa
        if not ok:
            aa_mismatch.append({"res": rid, "expected": aa, "got": got})
        hotspot_ok[rid] = {"present": True, "aa": got, "ok": ok}
    return {
        "n_residues": len(kept_ids),
        "range": [min(kept_ids), max(kept_ids)] if kept_ids else None,
        "contig_match": kept_ids == list(range(RES_START, RES_END + 1)),
        "missing_in_contig_span": [
            i for i in range(RES_START, RES_END + 1) if i not in present
        ],
        "hotspot_ok": hotspot_ok,
        "hotspot_missing": missing,
        "aa_mismatch": aa_mismatch,
        "pass": (
            not missing
            and not aa_mismatch
            and min(kept_ids) == RES_START
            and max(kept_ids) == RES_END
        ),
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    parser = PDBParser(QUIET=True)
    if not SRC.exists():
        raise SystemExit(f"Missing source receptor: {SRC}")
    src = parser.get_structure("5O45", str(SRC))
    src_ids = [r.id[1] for r in polymer_residues(src[0]["A"])]
    cropped, kept = crop_chain(src[0]["A"], "A", RES_START, RES_END)
    remarks = [
        f"source: {SRC.name}",
        f"crop: A{RES_START}-{RES_END} (RFpeptides contig)",
        f"extracted: {date.today().isoformat()}",
        f"removed_outside_crop: {[i for i in src_ids if i < RES_START or i > RES_END]}",
    ]
    write_pdb(cropped, OUT_PDB, remarks)
    # re-read for validation
    out_s = parser.get_structure("out", str(OUT_PDB))
    report = {
        "accessed": date.today().isoformat(),
        "source": str(SRC.relative_to(ROOT)),
        "output": str(OUT_PDB.relative_to(ROOT)),
        "contig": f"A{RES_START}-{RES_END}",
        "source_residue_range": [min(src_ids), max(src_ids)],
        "removed": [i for i in src_ids if i < RES_START or i > RES_END],
        "validation": validate(kept, out_s),
        "hotspots_json": str(HOTSPOTS_JSON.relative_to(ROOT)),
    }
    report_path = OUT_DIR / "receptor_prep_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Wrote {OUT_PDB}")
    print(f"Wrote {report_path}")
    print(f"pass={report['validation']['pass']} n={report['validation']['n_residues']}")
    if not report["validation"]["pass"]:
        raise SystemExit("Receptor prep validation FAILED")


if __name__ == "__main__":
    main()
