# zxdrfile

[![Zig](https://img.shields.io/badge/Zig-0.15.2+-f7a41d?logo=zig&logoColor=white)](https://ziglang.org/)
[![CI](https://github.com/N283T/zxdrfile/actions/workflows/ci.yml/badge.svg)](https://github.com/N283T/zxdrfile/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-BSD_2--Clause-blue.svg)](LICENSE)

A Zig library for reading GROMACS XDR trajectory files (XTC and TRR formats).

## Features

- **XTC reader** -- Compressed trajectory format with 3D coordinate decompression
- **TRR reader** -- Uncompressed trajectory format with coordinates, velocities, and forces
- **High performance** -- Buffered I/O with bulk reads and in-place byte-swapping
- **Zero dependencies** -- Pure Zig, no C bindings required

## Usage

Add as a Zig package dependency in your `build.zig.zon`:

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
mod.addImport("xdrfile", zxdrfile.module("xdrfile"));
```

### Reading XTC

```zig
const xdrfile = @import("xdrfile");

var reader = try xdrfile.XtcReader.open(allocator, "trajectory.xtc");
defer reader.close();

while (true) {
    var frame = reader.readFrame() catch |err| {
        if (err == xdrfile.XtcError.EndOfFile) break;
        return err;
    };
    defer frame.deinit(allocator);

    // frame.step, frame.time, frame.box, frame.coords, frame.precision
}
```

### Reading TRR

```zig
const xdrfile = @import("xdrfile");

var reader = try xdrfile.TrrReader.open(allocator, "trajectory.trr");
defer reader.close();

while (true) {
    var frame = reader.readFrame() catch |err| {
        if (err == xdrfile.TrrError.EndOfFile) break;
        return err;
    };
    defer frame.deinit(allocator);

    // frame.step, frame.time, frame.lambda, frame.box
    // frame.coords, frame.velocities, frame.forces (optional, nullable)
}
```

## Building

```bash
zig build test      # Run unit tests + validation tests
zig build validate  # Run validation tests only (against mdtraj reference)
zig build bench     # Run benchmarks (ReleaseFast)
```

## Requirements

- Zig 0.15.2 or later

> **Note:** Zig has not yet reached 1.0 and its standard library API changes
> frequently between versions. This library may not compile on versions other
> than the one specified above. Check the CI status badge for current compatibility.

## Differences from the Original C Library

This is not a line-by-line translation. Key differences:

- **Read-only API** -- Only reading is supported (no writing)
- **Zig-native error handling** -- Uses Zig's error union types instead of C-style return codes
- **Allocator-aware** -- All memory allocation goes through a caller-provided `std.mem.Allocator`
- **Buffered I/O** -- 64KB read buffer via `std.fs.File.Reader`, reducing syscall overhead
- **Bulk reads with in-place byte-swap** -- TRR vectors are read in a single call and byte-swapped in place, instead of one-element-at-a-time XDR decoding
- **Overflow-safe arithmetic** -- Uses `std.math.mul` for bounds checking on atom count calculations
- **Bounds checks on smallidx** -- Validates compression index against `FIRSTIDX`/`LASTIDX` to prevent out-of-bounds access

## Performance

Benchmarked on Apple M4 Pro, reading all frames from trajectory files (ReleaseFast).
C reference uses the original xdrfile library from mdtraj, compiled with `-O2`.

### XTC (compressed)

| File | Atoms | Frames | Size | zxdrfile | C (mdtraj) | Speedup |
|------|------:|-------:|-----:|---------:|-----------:|--------:|
| 3tvj_I | 531 | 1,001 | 2.4 MB | 113 MB/s | 184 MB/s | 0.6x |
| 5wvo_C | 3,858 | 1,001 | 17 MB | 261 MB/s | 246 MB/s | 1.1x |
| 6sup_A | 33,377 | 1,001 | 148 MB | 321 MB/s | 284 MB/s | 1.1x |

XTC performance is comparable to C. The decompression algorithm dominates runtime,
so I/O optimizations have limited impact.

### TRR (uncompressed)

| File | Atoms | Frames | Size | zxdrfile | C (mdtraj) | Speedup |
|------|------:|-------:|-----:|---------:|-----------:|--------:|
| 3tvj_I | 531 | 1,001 | 6.2 MB | 1,245 MB/s | 226 MB/s | **5.5x** |
| 5wvo_C | 3,858 | 1,001 | 44 MB | 5,256 MB/s | 199 MB/s | **26x** |
| 6sup_A | 33,377 | 1,001 | 383 MB | 8,298 MB/s | 231 MB/s | **36x** |

TRR is dramatically faster because the C library decodes each 4-byte value individually
via XDR function calls, while zxdrfile reads entire vectors in bulk and byte-swaps in place.

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

## License

BSD 2-Clause License. See [LICENSE](LICENSE) for details.
