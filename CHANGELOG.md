# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/).

## [0.3.0] - 2026-04-26

### Changed (BREAKING)

- **Zig 0.16 required**: bumped `minimum_zig_version` from `0.15.2` to `0.16.0`
- **Public API now takes `std.Io`**: as part of the Zig 0.16 unified I/O interface (`std.Io`), all reader/writer entry points have a new first parameter:
  - `XtcReader.open(io, allocator, path)`
  - `XtcWriter.open(io, allocator, path, natoms, mode)`
  - `TrrReader.open(io, allocator, path)`
  - `TrrWriter.open(io, allocator, path, natoms, mode)`
  Obtain `io` from `std.process.Init.io` in your `main` (recommended) or by initializing `std.Io.Threaded`.
- `converter` and `benchmark` migrated to the new `pub fn main(init: std.process.Init)` entry-point convention.
- Internal types: `std.fs.File` → `std.Io.File`, `std.io.Reader` → `std.Io.Reader`.

### Migration

```zig
// 0.2.x
var reader = try XtcReader.open(allocator, "traj.xtc");

// 0.3.0
var reader = try XtcReader.open(io, allocator, "traj.xtc");
```

In tests, use `std.testing.io` as the `io` argument.

## [0.2.0] - 2026-03-23

### Added

- **TRR writer** with frame-by-frame writing and append mode (#17)
- **XTC writer** with full coordinate compression algorithm (#18)
- Cross-format conversion tests and mdtraj validation script (#19)
- Converter CLI tool for TRR↔XTC conversion
- Append mode for both writers with natoms validation against existing files
- Frame data validation in writers (array length, null checks, natoms bounds)

### Changed

- Refactored XTC encoding utilities (`sizeofint`, `sizeofints`, `decodebits`, `decodeints`) from `XtcReader` methods to module-level functions for reader/writer reuse
- TRR writer header writes batched into a single 84-byte buffer write
- TRR `writeFloatsBulk` chunk size increased from 1024 to 4096 floats
- CI now runs validation tests as a separate sequential step to avoid parallel file conflicts
- README rewritten with complete reader/writer API documentation

### Fixed

- `close()` in both writers now releases all resources before returning flush errors (#17, #18)
- `encodebits` bit masks added to match C reference implementation (#18)
- XTC compression bounds check: guard `smallidx` against `LASTIDX` overflow (#18)
- Invalid precision (≤ 0) now returns error instead of silent fallback to 1000 (#18)
- Parallel test file conflicts resolved with unique tmp path prefixes per module (#21)

### Performance

- TRR writer: ~900 MB/s throughput (read+write combined on M4 Pro)
- XTC writer: ~280 MB/s throughput (compression-bound)

## [0.1.1] - 2026-03-11

### Changed

- Unify module name to `zxdrfile` (was `xdrfile` in build.zig, now matches package name)
- Update README usage examples to reflect new module name

## [0.1.0] - 2026-03-11

### Added

- XTC trajectory file reader with 3D coordinate decompression
- TRR trajectory file reader with coordinates, velocities, and forces support
- Validation tests comparing output against mdtraj reference data
- Benchmark infrastructure with coordinate validation
- C reference benchmark using mdtraj xdrfile for comparison
- README with usage examples, benchmarks, and acknowledgments

### Performance

- Buffered I/O (64KB) for both XTC and TRR readers
- TRR bulk read with in-place byte-swap (5x-36x faster than C xdrfile)

### Fixed

- Bounds check for `smallidx` after adjustment in XTC decompression
- Overflow-safe arithmetic for atom count calculations
- Proper `errdefer` cleanup in read paths
- Heap-allocated read buffer to avoid self-referencing pointer bug
