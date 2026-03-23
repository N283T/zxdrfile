#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "mdtraj",
#     "numpy",
# ]
# ///
"""Verify zxdrfile writer output using mdtraj as independent reader.

Writes test files with zxdrfile (via a Zig helper), then reads them
with mdtraj to verify correctness against the original files.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

import mdtraj as md
import numpy as np


def main():
    project_root = Path(__file__).parent.parent
    test_data = project_root / "test_data"

    print("=" * 60)
    print("Cross-format validation with mdtraj")
    print("=" * 60)

    # Test 1: Read original files with mdtraj as reference
    print("\n--- Loading reference files with mdtraj ---")
    trr_orig = md.load(str(test_data / "frame0.trr"), top=_make_topology(22))
    xtc_orig = md.load(str(test_data / "1l2y.xtc"), top=_make_topology(304))
    print(f"frame0.trr: {trr_orig.n_frames} frames, {trr_orig.n_atoms} atoms")
    print(f"1l2y.xtc:   {xtc_orig.n_frames} frames, {xtc_orig.n_atoms} atoms")

    # Test 2: Write TRR → XTC with zxdrfile, read with mdtraj
    print("\n--- Test: TRR → XTC (zxdrfile write, mdtraj read) ---")
    with tempfile.NamedTemporaryFile(suffix=".xtc", delete=False) as f:
        xtc_path = f.name

    _run_zig_converter(project_root, "trr_to_xtc",
                       str(test_data / "frame0.trr"), xtc_path)

    xtc_from_trr = md.load(xtc_path, top=_make_topology(22))
    print(f"  Written XTC: {xtc_from_trr.n_frames} frames, {xtc_from_trr.n_atoms} atoms")

    assert xtc_from_trr.n_frames == trr_orig.n_frames, \
        f"Frame count mismatch: {xtc_from_trr.n_frames} vs {trr_orig.n_frames}"

    # XTC is lossy — compare within precision tolerance
    max_diff = np.max(np.abs(xtc_from_trr.xyz - trr_orig.xyz))
    print(f"  Max coord diff (TRR vs XTC): {max_diff:.6f} nm")
    assert max_diff < 0.002, f"Coordinate difference too large: {max_diff}"
    print("  PASS")

    Path(xtc_path).unlink()

    # Test 3: Write XTC → TRR with zxdrfile, read with mdtraj
    print("\n--- Test: XTC → TRR (zxdrfile write, mdtraj read) ---")
    with tempfile.NamedTemporaryFile(suffix=".trr", delete=False) as f:
        trr_path = f.name

    _run_zig_converter(project_root, "xtc_to_trr",
                       str(test_data / "1l2y.xtc"), trr_path)

    trr_from_xtc = md.load(trr_path, top=_make_topology(304))
    print(f"  Written TRR: {trr_from_xtc.n_frames} frames, {trr_from_xtc.n_atoms} atoms")

    assert trr_from_xtc.n_frames == xtc_orig.n_frames, \
        f"Frame count mismatch: {trr_from_xtc.n_frames} vs {xtc_orig.n_frames}"

    # TRR is lossless after XTC decode — should match exactly
    max_diff = np.max(np.abs(trr_from_xtc.xyz - xtc_orig.xyz))
    print(f"  Max coord diff (XTC vs TRR): {max_diff:.6f} nm")
    assert max_diff < 0.0001, f"Coordinate difference too large: {max_diff}"
    print("  PASS")

    Path(trr_path).unlink()

    # Test 4: Full round-trip TRR → XTC → TRR, verify with mdtraj
    print("\n--- Test: TRR → XTC → TRR round-trip ---")
    with tempfile.NamedTemporaryFile(suffix=".xtc", delete=False) as f:
        xtc_mid = f.name
    with tempfile.NamedTemporaryFile(suffix=".trr", delete=False) as f:
        trr_final = f.name

    _run_zig_converter(project_root, "trr_to_xtc",
                       str(test_data / "frame0.trr"), xtc_mid)
    _run_zig_converter(project_root, "xtc_to_trr", xtc_mid, trr_final)

    trr_roundtrip = md.load(trr_final, top=_make_topology(22))
    print(f"  Round-trip TRR: {trr_roundtrip.n_frames} frames")

    max_diff = np.max(np.abs(trr_roundtrip.xyz - trr_orig.xyz))
    print(f"  Max coord diff (original vs round-trip): {max_diff:.6f} nm")
    assert max_diff < 0.002, f"Round-trip coordinate difference too large: {max_diff}"
    print("  PASS")

    Path(xtc_mid).unlink()
    Path(trr_final).unlink()

    print("\n" + "=" * 60)
    print("All cross-format validation tests PASSED")
    print("=" * 60)


def _make_topology(n_atoms: int) -> md.Topology:
    """Create a minimal topology with the given number of atoms."""
    top = md.Topology()
    chain = top.add_chain()
    res = top.add_residue("UNK", chain)
    for _ in range(n_atoms):
        top.add_atom("X", md.element.carbon, res)
    return top


def _run_zig_converter(project_root: Path, mode: str,
                       input_path: str, output_path: str):
    """Run the Zig converter helper."""
    result = subprocess.run(
        [str(project_root / "zig-out" / "bin" / "converter"),
         mode, input_path, output_path],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"Converter failed: {result.stderr}")
        sys.exit(1)


if __name__ == "__main__":
    main()
