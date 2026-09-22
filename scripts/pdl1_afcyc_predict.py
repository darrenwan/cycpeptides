#!/usr/bin/env python3
"""AfCycDesign smoke predictor for PD-L1 designed complexes (S2c).

Runs inside Docker image `pdl1-afcyc` (ColabDesign + AF2 params).
Campaign host scoring remains: uv run python scripts/pdl1_smoke_score.py

Usage (usually via scripts/pdl1_smoke_afcyc.sh):
  python scripts/pdl1_afcyc_predict.py \\
    --pdb path1.pdb --tag tag1 \\
    --params_dir /models/afcyc \\
    --out_dir /out/afcyc
"""
from __future__ import annotations

import argparse
import json
import traceback
from datetime import date
from pathlib import Path

import numpy as np
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import PDBIO, PDBParser

AA = {k.upper(): v for k, v in protein_letters_3to1.items()}
PARSER = PDBParser(QUIET=True)


def sequence_from_pdb(pdb: Path, chain_id: str) -> str:
    s = PARSER.get_structure("x", str(pdb))
    seq = []
    for r in s[0][chain_id]:
        if r.id[0] != " ":
            continue
        seq.append(AA.get(r.resname, "X"))
    return "".join(seq)


def chain_lengths(pdb: Path) -> dict[str, int]:
    s = PARSER.get_structure("x", str(pdb))
    out = {}
    for c in s[0]:
        n = sum(1 for r in c if r.id[0] == " ")
        if n:
            out[c.id] = n
    return out


def rewrite_chain_letters(pdb: Path, mapping: dict[str, str]) -> None:
    if all(a == b for a, b in mapping.items()):
        return
    s = PARSER.get_structure("x", str(pdb))
    # Bio.PDB cannot easily rename chains in-place; rewrite ATOM lines.
    text = pdb.read_text()
    lines = []
    for line in text.splitlines():
        if line.startswith(("ATOM", "HETATM", "TER")) and len(line) >= 22:
            cid = line[21]
            if cid in mapping and mapping[cid] != cid:
                line = line[:21] + mapping[cid] + line[22:]
        lines.append(line)
    pdb.write_text("\n".join(lines) + "\n")


def add_cyclic_offset(af_model, offset_type: int = 2) -> None:
    """AfCycDesign cyclic relative positional encoding (Rettie et al.)."""

    def cyclic_offset(L: int) -> np.ndarray:
        i = np.arange(L)
        ij = np.stack([i, i + L], -1)
        offset = i[:, None] - i[None, :]
        c_offset = np.abs(ij[:, None, :, None] - ij[None, :, None, :]).min((2, 3))
        if offset_type >= 2:
            a = c_offset < np.abs(offset)
            c_offset[a] = -c_offset[a]
        return c_offset * np.sign(offset)

    idx = af_model._inputs["residue_index"]
    offset = np.array(idx[:, None] - idx[None, :])
    if af_model.protocol == "binder":
        offset[af_model._target_len :, af_model._target_len :] = cyclic_offset(
            af_model._binder_len
        )
    af_model._inputs["offset"] = offset


def detect_pred_chains(
    pred_pdb: Path, n_binder: int, n_target: int, binder_chain: str, target_chain: str
) -> tuple[str, str]:
    lengths = chain_lengths(pred_pdb)
    binder_cands = [c for c, n in lengths.items() if n == n_binder]
    target_cands = [c for c, n in lengths.items() if n == n_target]
    if len(binder_cands) == 1 and len(target_cands) == 1:
        return binder_cands[0], target_cands[0]
    if lengths.get(binder_chain) == n_binder and lengths.get(target_chain) == n_target:
        return binder_chain, target_chain
    raise ValueError(f"cannot map pred chains: {lengths} vs binder={n_binder} target={n_target}")


def kabsch_rmsd(A: np.ndarray, B: np.ndarray) -> float:
    Ac, Bc = A.mean(0), B.mean(0)
    A0, B0 = A - Ac, B - Bc
    H = A0.T @ B0
    U, _S, Vt = np.linalg.svd(H)
    d = np.linalg.det(U @ Vt)
    D = np.diag([1.0, 1.0, 1.0 if d >= 0 else -1.0])
    R = U @ D @ Vt
    B2 = (R @ B0.T).T + Ac
    return float(np.sqrt(((A - B2) ** 2).sum() / len(A)))


def peptide_ca_rmsd_after_target_align(
    design: Path, pred: Path, binder_chain: str, target_chain: str
) -> float:
    sd = PARSER.get_structure("d", str(design))
    sp = PARSER.get_structure("p", str(pred))

    def cas(struct, cid):
        return np.array(
            [r["CA"].coord for r in struct[0][cid] if r.id[0] == " " and "CA" in r],
            dtype=float,
        )

    Td, Tp = cas(sd, target_chain), cas(sp, target_chain)
    Bd, Bp = cas(sd, binder_chain), cas(sp, binder_chain)
    if len(Td) != len(Tp) or len(Bd) != len(Bp) or len(Td) < 3:
        return float("nan")
    # Align pred target onto design target; apply same transform to binder
    Ac, Bc = Td.mean(0), Tp.mean(0)
    A0, B0 = Td - Ac, Tp - Bc
    H = A0.T @ B0
    U, _S, Vt = np.linalg.svd(H)
    d = np.linalg.det(U @ Vt)
    D = np.diag([1.0, 1.0, 1.0 if d >= 0 else -1.0])
    R = U @ D @ Vt
    Bp2 = (R @ (Bp - Bc).T).T + Ac
    return kabsch_rmsd(Bd, Bp2)


def run_one(
    pdb: Path,
    tag: str,
    out_dir: Path,
    params_dir: Path,
    binder_chain: str,
    target_chain: str,
    num_recycles: int,
    use_binder_template: bool,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    pred_path = out_dir / f"{tag}_afcyc_pred.pdb"
    row: dict = {
        "tag": tag,
        "input_pdb": str(pdb),
        "pred_pdb": str(pred_path),
        "status": "pending",
        "stage": "S2c_afcyc",
    }
    try:
        pep_seq = sequence_from_pdb(pdb, binder_chain)
        tgt_seq = sequence_from_pdb(pdb, target_chain)
        row["peptide_seq"] = pep_seq
        row["peptide_len"] = len(pep_seq)
        row["target_len"] = len(tgt_seq)

        from colabdesign import clear_mem, mk_afdesign_model  # type: ignore

        clear_mem()
        af_model = mk_afdesign_model(
            protocol="binder",
            num_models=1,
            sample_models=False,
            data_dir=str(params_dir),
        )
        af_model.prep_inputs(
            pdb_filename=str(pdb),
            chain=target_chain,
            binder_chain=binder_chain,
            use_binder_template=use_binder_template,
            rm_binder_seq=False,
            rm_binder_sc=False,
        )
        add_cyclic_offset(af_model, offset_type=2)
        af_model.set_seq(seq=pep_seq)
        af_model.predict(num_recycles=num_recycles, verbose=True)
        af_model.save_pdb(str(pred_path))

        src_b, src_t = detect_pred_chains(
            pred_path, len(pep_seq), len(tgt_seq), binder_chain, target_chain
        )
        rewrite_chain_letters(
            pred_path, {src_b: binder_chain, src_t: target_chain}
        )
        row["pred_chain_map"] = {
            "src_binder": src_b,
            "src_target": src_t,
            "dst_binder": binder_chain,
            "dst_target": target_chain,
        }

        aux = getattr(af_model, "aux", {}) or {}
        conf = {}
        for key in ("plddt", "ptm", "i_ptm"):
            if key in aux:
                conf[key] = float(np.mean(aux[key]))
        row["confidence"] = conf
        row["pep_ca_rmsd_vs_design"] = round(
            peptide_ca_rmsd_after_target_align(
                pdb, pred_path, binder_chain, target_chain
            ),
            3,
        )
        row["status"] = "ok"
    except Exception as exc:  # noqa: BLE001
        row["status"] = "failed"
        row["error"] = f"{type(exc).__name__}: {exc}"
        row["traceback"] = traceback.format_exc(limit=20)
    return row


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", action="append", required=True)
    ap.add_argument("--tag", action="append", required=True)
    ap.add_argument("--params_dir", type=Path, required=True)
    ap.add_argument("--out_dir", type=Path, required=True)
    ap.add_argument("--summary_json", type=Path, required=True)
    ap.add_argument("--binder_chain", default="A")
    ap.add_argument("--target_chain", default="B")
    ap.add_argument("--num_recycles", type=int, default=3)
    ap.add_argument("--use_binder_template", type=int, choices=(0, 1), default=1)
    args = ap.parse_args()
    if len(args.pdb) != len(args.tag):
        raise SystemExit("--pdb and --tag counts must match")

    rows = []
    for pdb_s, tag in zip(args.pdb, args.tag):
        print(f"=== AfCyc predict {tag} ===", flush=True)
        rows.append(
            run_one(
                Path(pdb_s),
                tag,
                args.out_dir,
                args.params_dir,
                args.binder_chain,
                args.target_chain,
                args.num_recycles,
                bool(args.use_binder_template),
            )
        )
        print(json.dumps({k: rows[-1].get(k) for k in ("tag", "status", "pep_ca_rmsd_vs_design", "confidence")}, ensure_ascii=False), flush=True)

    n_ok = sum(1 for r in rows if r["status"] == "ok")
    summary = {
        "accessed": date.today().isoformat(),
        "tool": "AfCycDesign",
        "n_input": len(rows),
        "n_ok": n_ok,
        "n_failed": len(rows) - n_ok,
        "rows": rows,
    }
    args.summary_json.parent.mkdir(parents=True, exist_ok=True)
    args.summary_json.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {args.summary_json} n_ok={n_ok}/{len(rows)}", flush=True)
    return 0 if n_ok > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
