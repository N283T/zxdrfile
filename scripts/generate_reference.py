#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "mdtraj",
# ]
# ///
"""Generate reference JSON from XTC/TRR test data using mdtraj.

Usage:
    ./scripts/generate_reference.py

Outputs JSON files to test_data/ for Zig validation tests.
"""

import json
from pathlib import Path

import mdtraj


def trajectory_to_reference(traj: mdtraj.Trajectory) -> dict:
    """Extract reference data from an mdtraj trajectory."""
    frames = []
    for i in range(traj.n_frames):
        frame_data = {
            "step": int(traj.time[i] / traj.timestep) if traj.timestep > 0 else i,
            "time": float(traj.time[i]),
            "box": traj.unitcell_vectors[i].tolist() if traj.unitcell_vectors is not None else None,
            "coords": traj.xyz[i].flatten().tolist(),
        }
        frames.append(frame_data)

    return {
        "natoms": traj.n_atoms,
        "n_frames": traj.n_frames,
        "frames": frames,
    }


def generate_xtc_reference(test_data_dir: Path) -> None:
    """Generate reference for 1l2y.xtc."""
    xtc_path = test_data_dir.joinpath("1l2y.xtc")
    if not xtc_path.exists():
        print(f"Skipping {xtc_path}: not found")
        return

    # mdtraj needs a topology; for raw XTC we create a minimal one
    top = mdtraj.Topology()
    # First read to get natoms
    with mdtraj.formats.XTCTrajectoryFile(str(xtc_path)) as f:
        xyz, time, step, box = f.read()

    n_atoms = xyz.shape[1]
    chain = top.add_chain()
    for _ in range(n_atoms):
        res = top.add_residue("UNK", chain)
        top.add_atom("X", mdtraj.element.virtual, res)

    traj = mdtraj.load(str(xtc_path), top=top)
    ref = trajectory_to_reference(traj)

    out_path = test_data_dir.joinpath("1l2y_xtc_reference.json")
    with open(out_path, "w") as f:
        json.dump(ref, f, indent=2)
    print(f"Written {out_path} ({ref['n_frames']} frames, {ref['natoms']} atoms)")


def generate_trr_reference(test_data_dir: Path) -> None:
    """Generate reference for frame0.trr."""
    trr_path = test_data_dir.joinpath("frame0.trr")
    if not trr_path.exists():
        print(f"Skipping {trr_path}: not found")
        return

    # Read raw TRR to get all data including velocities
    with mdtraj.formats.TRRTrajectoryFile(str(trr_path)) as f:
        xyz, time, step, box, lambd = f.read()

    frames = []
    for i in range(xyz.shape[0]):
        frame_data = {
            "step": int(step[i]),
            "time": float(time[i]),
            "lambda": float(lambd[i]),
            "box": box[i].tolist() if box is not None else None,
            "coords": xyz[i].flatten().tolist(),
        }
        frames.append(frame_data)

    ref = {
        "natoms": int(xyz.shape[1]),
        "n_frames": int(xyz.shape[0]),
        "frames": frames,
    }

    out_path = test_data_dir.joinpath("frame0_trr_reference.json")
    with open(out_path, "w") as f:
        json.dump(ref, f, indent=2)
    print(f"Written {out_path} ({ref['n_frames']} frames, {ref['natoms']} atoms)")


def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    test_data_dir = project_root.joinpath("test_data")

    generate_xtc_reference(test_data_dir)
    generate_trr_reference(test_data_dir)


if __name__ == "__main__":
    main()
