#!/usr/bin/env python3
"""Extract PD-L1 IgV from 5O4Y (peptide-71 complex) for design use.

In 5O4Y, ASU chains A/D/F are the macrocycle; B/C/E are PD-L1.
This script writes:
  - PDL1_IgV_5O4Y_chainA.pdb  (chain B polymer, rewritten as chain A; preferred — lowest B-factors)
  - PDL1_IgV_5O4Y_chain{B,C,E}.pdb  (original chain IDs, polymer only)

Usage:
  python scripts/extract_pdl1_5o4y_receptor.py
"""
from __future__ import annotations

from datetime import date
from pathlib import Path

from Bio.PDB import PDBIO, PDBParser, Structure, Model, Chain

ROOT = Path(__file__).resolve().parents[1]
PDB_DIR = ROOT / "data/targets/PDL1/pdb"
SRC = PDB_DIR / "5O4Y.pdb"
PREFERRED = "B"  # lowest mean B at hotspots among B/C/E


def extract(chain_src, out_chain_id: str, out_path: Path, remark_extra: str = "") -> None:
    new = Structure.Structure(out_path.stem)
    model = Model.Model(0)
    chain = Chain.Chain(out_chain_id)
    model.add(chain)
    new.add(model)
    for res in chain_src:
        if res.id[0] == " ":
            chain.add(res.copy())
    io = PDBIO()
    io.set_structure(new)
    tmp = out_path.with_suffix(".tmp.pdb")
    io.save(str(tmp))
    header = (
        f"HEADER    PD-L1 IGV FROM 5O4Y\n"
        f"TITLE     PD-L1 IgV from PDB 5O4Y (peptide-71), polymer only\n"
        f"REMARK    source_chain: {chain_src.id} -> output_chain: {out_chain_id}\n"
        f"REMARK    ASU map: protein=B/C/E ; peptide=A/D/F\n"
        f"REMARK    extracted: {date.today().isoformat()}\n"
    )
    if remark_extra:
        header += f"REMARK    {remark_extra}\n"
    body = "\n".join(
        l for l in tmp.read_text().splitlines() if not l.startswith("END")
    )
    out_path.write_text(header + body + "\nEND\n")
    tmp.unlink(missing_ok=True)


def main() -> None:
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("5O4Y", str(SRC))
    # Main design file: B -> A
    out_main = PDB_DIR / "PDL1_IgV_5O4Y_chainA.pdb"
    extract(
        s[0][PREFERRED],
        "A",
        out_main,
        remark_extra="preferred copy (lowest hotspot B-factors); Group II diversity receptor",
    )
    print(f"Wrote {out_main} (from chain {PREFERRED} -> A)")
    for cid in ("B", "C", "E"):
        out = PDB_DIR / f"PDL1_IgV_5O4Y_chain{cid}.pdb"
        extract(s[0][cid], cid, out)
        print(f"Wrote {out}")


if __name__ == "__main__":
    main()
