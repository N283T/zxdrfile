#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "mdtraj",
# ]
# ///
"""Convert XTC files to TRR format for benchmarking.

Reads XTC files from benchmarks/md_data/ and writes corresponding TRR files.
Since XTC has no topology, a minimal topology is created automatically.

Usage:
    ./scripts/convert_xtc_to_trr.py
"""

from pathlib import Path

import mdtraj


def make_topology(n_atoms: int) -> mdtraj.Topology:
    """Create a minimal topology with the given number of atoms."""
    top = mdtraj.Topology()
    chain = top.add_chain()
    for _ in range(n_atoms):
        res = top.add_residue("UNK", chain)
        top.add_atom("X", mdtraj.element.virtual, res)
    return top


def convert_xtc_to_trr(xtc_path: Path) -> None:
    """Convert a single XTC file to TRR."""
    trr_path = xtc_path.with_suffix(".trr")
    if trr_path.exists():
        print(f"  Skipping {xtc_path.name}: {trr_path.name} already exists")
        return

    # Read raw XTC to get natoms
    with mdtraj.formats.XTCTrajectoryFile(str(xtc_path)) as f:
        xyz, time, step, box = f.read()

    n_atoms = xyz.shape[1]
    top = make_topology(n_atoms)

    traj = mdtraj.load(str(xtc_path), top=top)
    traj.save_trr(str(trr_path))

    trr_size_mb = trr_path.stat().st_size / (1024 * 1024)
    print(f"  {xtc_path.name} -> {trr_path.name} ({n_atoms} atoms, {traj.n_frames} frames, {trr_size_mb:.1f} MB)")


def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    md_data_dir = project_root.joinpath("benchmarks", "md_data")

    xtc_files = sorted(md_data_dir.glob("*.xtc"))
    if not xtc_files:
        print(f"No XTC files found in {md_data_dir}")
        return

    print(f"Converting {len(xtc_files)} XTC files to TRR:")
    for xtc_path in xtc_files:
        convert_xtc_to_trr(xtc_path)

    print("Done.")


if __name__ == "__main__":
    main()
