# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/).

## [Archived] - 2026-05-21

### Changed

- Mark this standalone repository as archived. Ongoing trajectory I/O development
  is now managed in [ztraj](https://github.com/N283T/ztraj).

## [0.4.0] - 2026-04-26

### Changed

- **More accurate error reporting**: `XtcError` and `TrrError` gained
  `AccessDenied`, `IsDir`, `NoSpaceLeft`, and `IoError` variants.
  `XtcError` also gained `InvalidHeader` for truncated headers on append.
  Previously these were all collapsed into `FileNotFound`/`ReadError`/`WriteError`,
  which made debugging open/create failures impossible. (#25, #26, #27)
- `XtcReader.open` and `TrrReader.open` now reject directories at open time
  (via `allow_directory = false`) instead of failing later with a generic
  `ReadError`. (#25)
- `TrrWriter.open` append-mode now propagates `InvalidMagic` and `InvalidHeader`
  from the existing-file probe instead of collapsing them to `ReadError`. (#26)
- `XtcWriter.open` append-mode now reports `InvalidHeader` for truncated files
  instead of generic `ReadError`. (#26)
- `XtcReader.open` now classifies a non-positive natoms in the header as
  `InvalidAtomCount` instead of `ReadError`.
- `benchmark` now logs non-FileNotFound failures (corrupt JSON, permission
  errors) instead of silently skipping them. (#28)

### Added

- macOS to CI matrix and run cross-format tests on every push (#31, #29)
- `zig build converter-smoke` end-to-end smoke test of the converter binary
  exercising the full `std.process.Init` plumbing (#33, #30)
- IsDir tests for both readers and writers

### Compatibility

The new error variants are technically additive, but exhaustive switches over
`XtcError` / `TrrError` will need new arms — hence minor bump.

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
- Internal types and APIs:
  - `std.fs.File` / `std.fs.cwd()` → `std.Io.File` / `std.Io.Dir.cwd()`
  - `std.io.Reader` / `std.io.Writer` → `std.Io.Reader` / `std.Io.Writer`
  - `file.getEndPos()` → `file.length(io)`
  - `file.seekTo(0)` (on the File) → `reader.seekTo(0)` (on the buffered Reader)
  - `std.time.Timer` → `std.Io.Timestamp` with `Clock.awake`
  - `std.heap.GeneralPurposeAllocator` → `std.heap.DebugAllocator`
  - `std.process.argsAlloc/argsFree` → `std.process.Init.minimal.args.toSlice(arena)`

### Migration

```zig
// 0.2.x
var reader = try XtcReader.open(allocator, "traj.xtc");
var writer = try XtcWriter.open(allocator, "out.xtc", natoms, .write);

// 0.3.0
var reader = try XtcReader.open(io, allocator, "traj.xtc");
var writer = try XtcWriter.open(io, allocator, "out.xtc", natoms, .write);
```

In tests, use `std.testing.io` as the `io` argument. In binaries, take
`init: std.process.Init` as the `main` parameter and use `init.io`.

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
