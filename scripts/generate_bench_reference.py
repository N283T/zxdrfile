#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "mdtraj",
# ]
# ///
"""Generate reference JSON for benchmark data files.

Extracts sample frames (first, middle, last) to keep JSON small.
Validates that zxdrfile produces correct results on large files.

Usage:
    ./scripts/generate_bench_reference.py
"""

import json
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


def extract_sample_frames(
    path: Path,
    fmt: str,
) -> dict | None:
    """Extract sample frames from a trajectory file.

    Returns reference dict with first, middle, and last frames.
    """
    if not path.exists():
        print(f"  Skipping {path.name}: not found")
        return None

    if fmt == "xtc":
        with mdtraj.formats.XTCTrajectoryFile(str(path)) as f:
            xyz, time, step, box = f.read()
        lambd = None
    elif fmt == "trr":
        with mdtraj.formats.TRRTrajectoryFile(str(path)) as f:
            xyz, time, step, box, lambd = f.read()
    else:
        raise ValueError(f"Unknown format: {fmt}")

    n_frames = xyz.shape[0]
    n_atoms = xyz.shape[1]

    # Sample frames: first, middle, last
    sample_indices = sorted(set([0, n_frames // 2, n_frames - 1]))

    frames = {}
    for idx in sample_indices:
        frame_data = {
            "time": float(time[idx]),
            "coords": xyz[idx].flatten().tolist(),
        }
        if box is not None:
            frame_data["box"] = box[idx].tolist()
        if lambd is not None:
            frame_data["lambda"] = float(lambd[idx])
        frames[str(idx)] = frame_data

    return {
        "natoms": n_atoms,
        "n_frames": n_frames,
        "format": fmt,
        "source": path.name,
        "sample_frames": frames,
    }


def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    bench_data_dir = project_root.joinpath("benchmarks", "md_data")
    out_dir = project_root.joinpath("benchmarks", "reference")
    out_dir.mkdir(parents=True, exist_ok=True)

    files = [
        ("3tvj_I_R1.xtc", "xtc"),
        ("5wvo_C_R1.xtc", "xtc"),
        ("6sup_A_R1.xtc", "xtc"),
        ("3tvj_I_R1.trr", "trr"),
        ("5wvo_C_R1.trr", "trr"),
        ("6sup_A_R1.trr", "trr"),
    ]

    print(f"Generating benchmark references from {bench_data_dir}:")
    for filename, fmt in files:
        path = bench_data_dir.joinpath(filename)
        ref = extract_sample_frames(path, fmt)
        if ref is None:
            continue

        out_name = path.stem + f"_{fmt}_reference.json"
        out_path = out_dir.joinpath(out_name)
        with open(out_path, "w") as f:
            json.dump(ref, f, indent=2)

        n_samples = len(ref["sample_frames"])
        print(
            f"  {filename} -> {out_name}"
            f" ({ref['natoms']} atoms, {ref['n_frames']} frames,"
            f" {n_samples} samples)"
        )

    print("Done.")


if __name__ == "__main__":
    main()
