# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/).

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
