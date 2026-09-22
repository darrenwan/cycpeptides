# Optional PyRosetta wheel cache

This directory is **not required** for the PD-L1 S2 smoke path (RFpeptides / MPNN / AfCyc).

Used only when building `pdl1-pyrosetta`:

```bash
# Place the official cp310 linux wheel here, then:
bash scripts/build_academic_images.sh pyrosetta
```

Expected filename (override with `PYROSETTA_WHEEL`):

- `pyrosetta-0-cp310-cp310-linux_x86_64.whl`

If the file is missing, the Dockerfile downloads from RosettaCommons quarterly
find-links instead. Override the directory with `PYROSETTA_WHEELS_DIR`.

Do not commit the `.whl` (~1.5 GiB); it is gitignored.
