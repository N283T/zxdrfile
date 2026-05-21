# zxdrfile

> **Archived:** This repository is now in maintenance/archive mode.
> XTC/TRR trajectory I/O is managed as part of [ztraj](https://github.com/N283T/ztraj),
> which is the recommended place for new development, bug fixes, and usage.
> This repository remains available as a historical standalone Zig xdrfile implementation.

[![Zig](https://img.shields.io/badge/Zig-0.16.0+-f7a41d?logo=zig&logoColor=white)](https://ziglang.org/)
[![CI](https://github.com/N283T/zxdrfile/actions/workflows/ci.yml/badge.svg)](https://github.com/N283T/zxdrfile/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-BSD_2--Clause-blue.svg)](LICENSE)

A high-performance Zig library for reading and writing GROMACS XDR trajectory files (XTC and TRR formats).

[**API Documentation**](https://n283t.github.io/zxdrfile/)

## Repository Status

This standalone package has been archived because trajectory I/O is now managed
in the unified [ztraj](https://github.com/N283T/ztraj) repository. Use ztraj for
ongoing development and new integrations.

The code here is kept for reference and for existing users pinned to zxdrfile,
but no new feature work is planned in this repository.

## Features

- **XTC reader/writer** -- Compressed coordinate trajectories with lossy 3D compression
- **TRR reader/writer** -- Full-precision trajectories with coordinates, velocities, and forces
- **High performance** -- Buffered I/O, bulk byte-swapping, batched header writes
- **Zero dependencies** -- Pure Zig implementation, no C bindings required
- **Append mode** -- Append frames to existing trajectory files with natoms validation

## Installation

> **For new projects:** prefer [ztraj](https://github.com/N283T/ztraj).
> The instructions below are retained for existing users who intentionally depend
> on the archived standalone `zxdrfile` package.

Add as a dependency in your `build.zig.zon`:

```zon
.dependencies = .{
    .zxdrfile = .{
        .url = "https://github.com/N283T/zxdrfile/archive/<commit>.tar.gz",
        .hash = "...",
    },
},
```

Then in your `build.zig`:

```zig
const zxdrfile = b.dependency("zxdrfile", .{ .target = target, .optimize = optimize });
mod.addImport("zxdrfile", zxdrfile.module("zxdrfile"));
```

## Quick Start

All reader/writer entry points take `io: std.Io` as their first argument
(Zig 0.16's unified I/O interface). Obtain `io` from `std.process.Init.io`
in your `main`:

```zig
pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const allocator = init.gpa;

    // Command-line args (sentinel-terminated strings):
    const args = try init.minimal.args.toSlice(init.arena.allocator());
    _ = args;

    // ... use io and allocator with the readers/writers below ...
}
```

In tests, use `std.testing.io`.

### Reading

```zig
const xdrfile = @import("zxdrfile");

// XTC
var reader = try xdrfile.XtcReader.open(io, allocator, "trajectory.xtc");
defer reader.close();

while (true) {
    var frame = reader.readFrame() catch |err| {
        if (err == xdrfile.XtcError.EndOfFile) break;
        return err;
    };
    defer frame.deinit(allocator);
    // frame.step, frame.time, frame.box, frame.coords, frame.precision
}

// TRR
var reader = try xdrfile.TrrReader.open(io, allocator, "trajectory.trr");
defer reader.close();

while (true) {
    var frame = reader.readFrame() catch |err| {
        if (err == xdrfile.TrrError.EndOfFile) break;
        return err;
    };
    defer frame.deinit(allocator);
    // frame.step, frame.time, frame.lambda, frame.box
    // frame.coords, frame.velocities, frame.forces (optional)
}
```

### Writing

```zig
const xdrfile = @import("zxdrfile");

// XTC (lossy compression — precision controls accuracy)
var writer = try xdrfile.XtcWriter.open(io, allocator, "output.xtc", natoms, .write);
defer writer.close() catch {};

try writer.writeFrame(.{
    .step = 0,
    .time = 0.0,
    .box = box,
    .coords = coords,
    .precision = 1000.0, // ~0.001 nm accuracy
});

// TRR (lossless)
var writer = try xdrfile.TrrWriter.open(io, allocator, "output.trr", natoms, .write);
defer writer.close() catch {};

try writer.writeFrame(.{
    .step = 0,
    .time = 0.0,
    .lambda = 0.0,
    .box = box,
    .has_x = true,
    .has_v = false,
    .has_f = false,
    .coords = coords,
    .velocities = null,
    .forces = null,
});

// Append to existing file
var writer = try xdrfile.TrrWriter.open(io, allocator, "existing.trr", natoms, .append);
```

## Building

```bash
zig build test         # Run unit tests
zig build validate     # Run validation tests (against mdtraj reference)
zig build cross-format # Run cross-format conversion tests
zig build bench        # Run benchmarks (ReleaseFast)
```

## Requirements

- Zig 0.16.0 or later

> **Note:** Zig has not yet reached 1.0 and its standard library API changes
> frequently between versions. This library may not compile on versions other
> than the one specified above. Check the CI status badge for current compatibility.

## Performance

Benchmarked on Apple M4 Pro, reading all frames from trajectory files (ReleaseFast).
C reference uses the original xdrfile library from mdtraj, compiled with `-O2`.

### XTC (compressed)

| File | Atoms | Frames | Size | zxdrfile | C (mdtraj) | Speedup |
|------|------:|-------:|-----:|---------:|-----------:|--------:|
| 3tvj_I | 531 | 1,001 | 2.4 MB | 113 MB/s | 184 MB/s | 0.6x |
| 5wvo_C | 3,858 | 1,001 | 17 MB | 261 MB/s | 246 MB/s | 1.1x |
| 6sup_A | 33,377 | 1,001 | 148 MB | 321 MB/s | 284 MB/s | 1.1x |

XTC performance is comparable to C. The compression algorithm dominates runtime,
so I/O optimizations have limited impact.

### TRR (uncompressed)

| File | Atoms | Frames | Size | zxdrfile | C (mdtraj) | Speedup |
|------|------:|-------:|-----:|---------:|-----------:|--------:|
| 3tvj_I | 531 | 1,001 | 6.2 MB | 1,245 MB/s | 226 MB/s | **5.5x** |
| 5wvo_C | 3,858 | 1,001 | 44 MB | 5,256 MB/s | 199 MB/s | **26x** |
| 6sup_A | 33,377 | 1,001 | 383 MB | 8,298 MB/s | 231 MB/s | **36x** |

TRR is dramatically faster because the C library decodes each 4-byte value individually
via XDR function calls, while zxdrfile reads entire vectors in bulk and byte-swaps in place.

## Differences from the Original C Library

This is not a line-by-line translation. Key differences:

- **Full read/write support** -- Both XTC and TRR formats support reading and writing
- **Zig-native error handling** -- Uses Zig's error union types instead of C-style return codes
- **Allocator-aware** -- All memory allocation goes through a caller-provided `std.mem.Allocator`
- **Buffered I/O** -- 64KB buffer for both reads and writes, reducing syscall overhead
- **Bulk byte-swapping** -- Reads/writes entire vectors at once with in-place byte-swap
- **Batched header writes** -- TRR header packed into a single buffer write
- **Overflow-safe arithmetic** -- Uses `std.math.mul` for bounds checking on atom count calculations
- **Compression bounds checks** -- Validates `smallidx` against `FIRSTIDX`/`LASTIDX` to prevent out-of-bounds access

## Acknowledgments

This library is a Zig port of the xdrfile C library from the
[mdtraj](https://github.com/mdtraj/mdtraj) project, specifically the BSD-licensed
source files in
[`mdtraj/formats/xtc/src/`](https://github.com/mdtraj/mdtraj/tree/master/mdtraj/formats/xtc/src):

- `xdrfile.c` / `xdrfile_xtc.c` / `xdrfile_trr.c`

The original xdrfile library was created by **Erik Lindahl** and **David van der Spoel**
as part of the [GROMACS](https://www.gromacs.org/) project, with modifications by
**Robert T. McGibbon** for mdtraj.

Thanks to all the original authors for making this code available under the BSD 2-Clause license.

This library is based on the following projects:

- McGibbon, R. T. et al. "MDTraj: A Modern Open Library for the Analysis of Molecular Dynamics Trajectories." *Biophys. J.* **109**, 1528–1532 (2015). [doi:10.1016/j.bpj.2015.08.015](https://doi.org/10.1016/j.bpj.2015.08.015)
- Abraham, M. J. et al. "GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers." *SoftwareX* **1–2**, 19–25 (2015). [doi:10.1016/j.softx.2015.06.001](https://doi.org/10.1016/j.softx.2015.06.001)

## License

BSD 2-Clause License. See [LICENSE](LICENSE) for details.
