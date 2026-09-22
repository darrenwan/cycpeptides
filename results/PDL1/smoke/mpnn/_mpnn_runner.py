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
