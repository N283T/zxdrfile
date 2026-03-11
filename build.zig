const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const mod = b.addModule("xdrfile", .{
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
    test_step.dependOn(&run_validation.step);
}
