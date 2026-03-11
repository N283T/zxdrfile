# zxdrfile

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
