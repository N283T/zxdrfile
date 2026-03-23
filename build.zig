const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const mod = b.addModule("zxdrfile", .{
        .root_source_file = b.path("src/xdrfile.zig"),
        .target = target,
        .optimize = optimize,
    });

    // Tests
    const tests = b.addTest(.{ .root_module = mod });
    const run_tests = b.addRunArtifact(tests);
    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_tests.step);

    // Validation tests (compare against mdtraj reference)
    const validation_mod = b.addModule("validation_test", .{
        .root_source_file = b.path("src/validation_test.zig"),
        .target = target,
        .optimize = optimize,
    });
    const validation_tests = b.addTest(.{ .root_module = validation_mod });
    const run_validation = b.addRunArtifact(validation_tests);
    const validate_step = b.step("validate", "Run validation tests against mdtraj reference");
    validate_step.dependOn(&run_validation.step);
    // Not included in `zig build test` to avoid tmp file conflicts
    // from parallel test binaries. Run separately with `zig build validate`.

    // Cross-format conversion tests
    const cross_format_mod = b.addModule("cross_format_test", .{
        .root_source_file = b.path("src/cross_format_test.zig"),
        .target = target,
        .optimize = optimize,
    });
    const cross_format_tests = b.addTest(.{ .root_module = cross_format_mod });
    const run_cross_format = b.addRunArtifact(cross_format_tests);
    const cross_format_step = b.step("cross-format", "Run cross-format conversion tests");
    cross_format_step.dependOn(&run_cross_format.step);
    // Not included in `zig build test` to avoid tmp file conflicts
    // from parallel test binaries. Run separately with `zig build cross-format`.

    // Documentation
    const lib = b.addLibrary(.{
        .name = "zxdrfile",
        .root_module = mod,
    });
    const install_docs = b.addInstallDirectory(.{
        .source_dir = lib.getEmittedDocs(),
        .install_dir = .prefix,
        .install_subdir = "docs",
    });
    const docs_step = b.step("docs", "Generate documentation");
    docs_step.dependOn(&install_docs.step);

    // Converter tool (for cross-validation with mdtraj)
    const converter_mod = b.addModule("converter", .{
        .root_source_file = b.path("src/converter.zig"),
        .target = target,
        .optimize = optimize,
    });
    const converter_exe = b.addExecutable(.{
        .name = "converter",
        .root_module = converter_mod,
    });
    b.installArtifact(converter_exe);

    // Benchmark (always ReleaseFast for meaningful results)
    const bench_mod = b.addModule("benchmark", .{
        .root_source_file = b.path("src/benchmark.zig"),
        .target = target,
        .optimize = .ReleaseFast,
    });
    const bench_exe = b.addExecutable(.{
        .name = "benchmark",
        .root_module = bench_mod,
    });
    b.installArtifact(bench_exe);
    const run_bench = b.addRunArtifact(bench_exe);
    const bench_step = b.step("bench", "Run benchmarks (ReleaseFast)");
    bench_step.dependOn(&run_bench.step);
}
