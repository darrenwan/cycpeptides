#!/usr/bin/env bash
# S2 smoke: ProteinMPNN on RFpeptides backbones (Docker RFD python + host MPNN code).
# Does NOT use system/conda python for campaign logic; torch lives in the RFD image.
#
# Usage (after rfpeptides smoke):
#   bash scripts/pdl1_smoke_mpnn.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BB_DIR="${ROOT}/results/PDL1/smoke/rfpeptides"
OUT_DIR="${ROOT}/results/PDL1/smoke/mpnn"
LOG_DIR="${ROOT}/results/PDL1/smoke/logs"
MPNN_DIR="${PROTEINMPNN_DIR:-${ROOT}/third_party/ProteinMPNN}"
IMAGE="${RFD_IMAGE:-rosettacommons/rfdiffusion:latest}"
SEQS_PER_BB="${SEQS_PER_BB:-2}"
DESIGN_CHAIN="${DESIGN_CHAIN:-A}"  # cyclic binder is chain A in RFD outputs
FIXED_CHAINS="${FIXED_CHAINS:-B}" # PD-L1

mkdir -p "${OUT_DIR}/seqs" "${OUT_DIR}/designed_pdbs" "${LOG_DIR}"
LOG="${LOG_DIR}/mpnn_smoke_$(date +%Y%m%d_%H%M%S).log"

shopt -s nullglob
BBS=("${BB_DIR}"/pdl1_cyc_*.pdb)
if [[ ${#BBS[@]} -lt 1 ]]; then
  echo "No backbones in ${BB_DIR}. Run bash scripts/pdl1_smoke_rfpeptides.sh first." >&2
  exit 1
fi
if [[ ! -f "${MPNN_DIR}/protein_mpnn_run.py" ]]; then
  echo "Missing ProteinMPNN at ${MPNN_DIR}" >&2
  echo "Run: bash scripts/bootstrap_local_assets.sh" >&2
  exit 1
fi

# Lightweight runner inside RFD image python (has torch)
RUNNER="${OUT_DIR}/_mpnn_runner.py"
cat > "${RUNNER}" <<'PY'
import glob, json, os, subprocess, sys
from pathlib import Path

bb_dir = Path("/backbones")
out_dir = Path("/mpnn_out")
mpnn = Path("/ProteinMPNN")
seqs_per = int(os.environ.get("SEQS_PER_BB", "2"))
design_chain = os.environ.get("DESIGN_CHAIN", "A")
# ProteinMPNN expects a JSONL of parsed PDBs — use helper if present, else call run with pdb path loop.

pdbs = sorted(bb_dir.glob("pdl1_cyc_*.pdb"))
print(f"n_backbones={len(pdbs)} seqs_per_bb={seqs_per}", flush=True)
out_dir.mkdir(parents=True, exist_ok=True)
(out_dir / "seqs").mkdir(exist_ok=True)

# Prefer official helper parse script if available
parse = mpnn / "helper_scripts" / "parse_multiple_chains.py"
assign = mpnn / "helper_scripts" / "assign_fixed_chains.py"
parsed = out_dir / "parsed_pdbs.jsonl"
fixed = out_dir / "assigned_pdbs.jsonl"

pdbs_dir = out_dir / "pdbs"
pdbs_dir.mkdir(exist_ok=True)
for p in pdbs:
    (pdbs_dir / p.name).write_bytes(p.read_bytes())

cmd_parse = [
    sys.executable, str(parse),
    "--input_path", str(pdbs_dir),
    "--output_path", str(parsed),
]
print("RUN", " ".join(cmd_parse), flush=True)
subprocess.check_call(cmd_parse)

cmd_assign = [
    sys.executable, str(assign),
    "--input_path", str(parsed),
    "--output_path", str(fixed),
    "--chain_list", design_chain,
]
print("RUN", " ".join(cmd_assign), flush=True)
subprocess.check_call(cmd_assign)

cmd = [
    sys.executable, str(mpnn / "protein_mpnn_run.py"),
    "--jsonl_path", str(parsed),
    "--chain_id_jsonl", str(fixed),
    "--out_folder", str(out_dir),
    "--num_seq_per_target", str(seqs_per),
    "--sampling_temp", "0.1",
    "--seed", "37",
    "--batch_size", "1",
    "--omit_AAs", "C",
]
print("RUN", " ".join(cmd), flush=True)
subprocess.check_call(cmd)
print("MPNN done", flush=True)
PY

{
  echo "=== S2 ProteinMPNN smoke ==="
  echo "backbones=${#BBS[@]} seqs_per_bb=${SEQS_PER_BB}"
  echo "mpnn=${MPNN_DIR} image=${IMAGE}"
} | tee "${LOG}"

set +e
docker run --rm --gpus all \
  --entrypoint bash \
  -e SEQS_PER_BB="${SEQS_PER_BB}" \
  -e DESIGN_CHAIN="${DESIGN_CHAIN}" \
  -v "${BB_DIR}:/backbones:ro" \
  -v "${MPNN_DIR}:/ProteinMPNN:ro" \
  -v "${OUT_DIR}:/mpnn_out" \
  "${IMAGE}" \
  -lc '/app/RFdiffusion/.venv/bin/python /mpnn_out/_mpnn_runner.py' \
  2>&1 | tee -a "${LOG}"
rc=${PIPESTATUS[0]}
set -e

# Stitch sequences back onto backbone PDBs with uv (biopython)
uv run --project "${ROOT}" python - <<PY
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain
from Bio.Data.IUPACData import protein_letters_1to3
import re

root = Path("${ROOT}")
bb_dir = root / "results/PDL1/smoke/rfpeptides"
seq_dir = root / "results/PDL1/smoke/mpnn/seqs"
out_dir = root / "results/PDL1/smoke/mpnn/designed_pdbs"
out_dir.mkdir(parents=True, exist_ok=True)
aa1to3 = {k.upper(): v.upper() for k, v in protein_letters_1to3.items()}
parser = PDBParser(QUIET=True)

def set_seq_on_chain(chain, seq: str):
    residues = [r for r in chain if r.id[0] == " "]
    if len(residues) != len(seq):
        raise ValueError(f"len mismatch {len(residues)} vs {len(seq)}")
    for r, aa in zip(residues, seq):
        r.resname = aa1to3.get(aa, "UNK")

n = 0
for fa in sorted(seq_dir.glob("*.fa")):
    text = fa.read_text()
    # ProteinMPNN writes recovery + samples; keep sample lines (not original)
    blocks = re.split(r">", text)
    samples = []
    for b in blocks:
        b = b.strip()
        if not b:
            continue
        header, *seq_lines = b.splitlines()
        seq = "".join(seq_lines).replace(" ", "").strip()
        if "sample=" in header or "T=" in header:
            samples.append(seq)
        elif samples == [] and seq and "score=" not in header.lower():
            # sometimes first is native recovery — skip if marked recovered
            if "recovery" in header.lower() or "native" in header.lower():
                continue
    # Fallback: all sequences after first
    if not samples:
        seqs = []
        for b in blocks:
            b = b.strip()
            if not b:
                continue
            header, *seq_lines = b.splitlines()
            seq = "".join(seq_lines).replace(" ", "").strip()
            if seq:
                seqs.append(seq)
        samples = seqs[1:] if len(seqs) > 1 else seqs

    stem = fa.stem  # e.g. pdl1_cyc_0
    bb = bb_dir / f"{stem}.pdb"
    if not bb.exists():
        print("skip missing bb", bb)
        continue
    for i, seq in enumerate(samples, 1):
        s = parser.get_structure("x", str(bb))
        # binder is typically chain A
        if "A" not in s[0]:
            print("no chain A", bb)
            continue
        try:
            set_seq_on_chain(s[0]["A"], seq)
        except ValueError as e:
            print("seq stitch fail", fa.name, e)
            continue
        out = out_dir / f"{stem}_sample{i}.pdb"
        io = PDBIO()
        io.set_structure(s)
        io.save(str(out))
        n += 1
print(f"wrote_designed_pdbs={n}")
PY

n_des="$(find "${OUT_DIR}/designed_pdbs" -name '*.pdb' 2>/dev/null | wc -l)"
echo "exit=${rc} n_designed_pdbs=${n_des}" | tee -a "${LOG}"
if [[ "${n_des}" -lt 1 ]]; then
  echo "MPNN smoke FAILED (no designed PDBs)" >&2
  exit 1
fi
echo "MPNN smoke OK → ${OUT_DIR}/designed_pdbs"
