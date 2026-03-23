// Cross-format conversion tests.
//
// Validates that data survives format conversion round-trips:
// - TRR → XTC → read back (lossy: XTC precision applies)
// - XTC → TRR → read back (lossless after initial XTC decode)
// - TRR → XTC → TRR → compare with original
// - XTC → TRR → XTC → compare with original
//
const std = @import("std");
const xdrfile = @import("xdrfile.zig");

const XtcReader = xdrfile.XtcReader;
const XtcWriter = xdrfile.XtcWriter;
const XtcFrame = xdrfile.XtcFrame;
const XtcError = xdrfile.XtcError;

const TrrReader = xdrfile.TrrReader;
const TrrWriter = xdrfile.TrrWriter;
const TrrFrame = xdrfile.TrrFrame;
const TrrError = xdrfile.TrrError;

test "TRR → XTC → read back" {
    const allocator = std.testing.allocator;
    const xtc_path = "test_data/tmp_trr_to_xtc.xtc";

    // Read TRR and write as XTC
    const precision: f32 = 1000.0;
    var natoms: i32 = 0;
    var frame_count: usize = 0;
    {
        var reader = try TrrReader.open(allocator, "test_data/frame0.trr");
        defer reader.close();
        natoms = reader.getNumAtoms();

        var writer = try XtcWriter.open(allocator, xtc_path, natoms, .write);
        defer writer.close() catch {};

        while (true) {
            var frame = reader.readFrame() catch |err| {
                if (err == TrrError.EndOfFile) break;
                return err;
            };
            defer frame.deinit(allocator);

            const xtc_frame = XtcFrame{
                .step = frame.step,
                .time = frame.time,
                .box = frame.box,
                .coords = frame.coords.?,
                .precision = precision,
            };
            try writer.writeFrame(xtc_frame);
            frame_count += 1;
        }
    }
    defer std.fs.cwd().deleteFile(xtc_path) catch {};

    try std.testing.expectEqual(@as(usize, 501), frame_count);

    // Read back XTC and verify
    var xtc_reader = try XtcReader.open(allocator, xtc_path);
    defer xtc_reader.close();

    try std.testing.expectEqual(natoms, xtc_reader.getNumAtoms());

    // Also re-read TRR for comparison
    var trr_reader = try TrrReader.open(allocator, "test_data/frame0.trr");
    defer trr_reader.close();

    const tolerance: f32 = 1.0 / precision + 0.001;
    var read_count: usize = 0;
    while (true) {
        var xtc_frame = xtc_reader.readFrame() catch |err| {
            if (err == XtcError.EndOfFile) break;
            return err;
        };
        defer xtc_frame.deinit(allocator);

        var trr_frame = trr_reader.readFrame() catch |err| {
            return err;
        };
        defer trr_frame.deinit(allocator);

        // Step and time should match exactly
        try std.testing.expectEqual(trr_frame.step, xtc_frame.step);
        try std.testing.expectApproxEqAbs(trr_frame.time, xtc_frame.time, 0.001);

        // Coords within XTC precision tolerance
        const trr_coords = trr_frame.coords.?;
        for (trr_coords, xtc_frame.coords) |t, x| {
            try std.testing.expectApproxEqAbs(t, x, tolerance);
        }

        read_count += 1;
    }
    try std.testing.expectEqual(@as(usize, 501), read_count);
}

test "XTC → TRR → read back" {
    const allocator = std.testing.allocator;
    const trr_path = "test_data/tmp_xtc_to_trr.trr";

    // Read XTC and write as TRR
    var natoms: i32 = 0;
    var frame_count: usize = 0;
    {
        var reader = try XtcReader.open(allocator, "test_data/1l2y.xtc");
        defer reader.close();
        natoms = reader.getNumAtoms();

        var writer = try TrrWriter.open(allocator, trr_path, natoms, .write);
        defer writer.close() catch {};

        while (true) {
            var frame = reader.readFrame() catch |err| {
                if (err == XtcError.EndOfFile) break;
                return err;
            };
            defer frame.deinit(allocator);

            const trr_frame = TrrFrame{
                .step = frame.step,
                .time = frame.time,
                .lambda = 0.0,
                .box = frame.box,
                .has_x = true,
                .has_v = false,
                .has_f = false,
                .coords = frame.coords,
                .velocities = null,
                .forces = null,
            };
            try writer.writeFrame(trr_frame);
            frame_count += 1;
        }
    }
    defer std.fs.cwd().deleteFile(trr_path) catch {};

    try std.testing.expectEqual(@as(usize, 38), frame_count);

    // Read back TRR and compare with original XTC
    var xtc_reader = try XtcReader.open(allocator, "test_data/1l2y.xtc");
    defer xtc_reader.close();

    var trr_reader = try TrrReader.open(allocator, trr_path);
    defer trr_reader.close();

    try std.testing.expectEqual(natoms, trr_reader.getNumAtoms());

    // TRR is lossless, so XTC→TRR should preserve exact float values
    const tolerance: f32 = 0.0001;
    var read_count: usize = 0;
    while (true) {
        var xtc_frame = xtc_reader.readFrame() catch |err| {
            if (err == XtcError.EndOfFile) break;
            return err;
        };
        defer xtc_frame.deinit(allocator);

        var trr_frame = trr_reader.readFrame() catch |err| {
            return err;
        };
        defer trr_frame.deinit(allocator);

        try std.testing.expectEqual(xtc_frame.step, trr_frame.step);
        try std.testing.expectApproxEqAbs(xtc_frame.time, trr_frame.time, tolerance);

        const trr_coords = trr_frame.coords.?;
        for (xtc_frame.coords, trr_coords) |x, t| {
            try std.testing.expectApproxEqAbs(x, t, tolerance);
        }

        read_count += 1;
    }
    try std.testing.expectEqual(@as(usize, 38), read_count);
}

test "TRR → XTC → TRR round-trip" {
    const allocator = std.testing.allocator;
    const xtc_path = "test_data/tmp_roundtrip_xtc.xtc";
    const trr_path = "test_data/tmp_roundtrip_trr.trr";
    const precision: f32 = 1000.0;

    var natoms: i32 = 0;

    // Step 1: TRR → XTC
    {
        var reader = try TrrReader.open(allocator, "test_data/frame0.trr");
        defer reader.close();
        natoms = reader.getNumAtoms();

        var writer = try XtcWriter.open(allocator, xtc_path, natoms, .write);
        defer writer.close() catch {};

        while (true) {
            var frame = reader.readFrame() catch |err| {
                if (err == TrrError.EndOfFile) break;
                return err;
            };
            defer frame.deinit(allocator);
            try writer.writeFrame(XtcFrame{
                .step = frame.step,
                .time = frame.time,
                .box = frame.box,
                .coords = frame.coords.?,
                .precision = precision,
            });
        }
    }
    defer std.fs.cwd().deleteFile(xtc_path) catch {};

    // Step 2: XTC → TRR
    {
        var reader = try XtcReader.open(allocator, xtc_path);
        defer reader.close();

        var writer = try TrrWriter.open(allocator, trr_path, natoms, .write);
        defer writer.close() catch {};

        while (true) {
            var frame = reader.readFrame() catch |err| {
                if (err == XtcError.EndOfFile) break;
                return err;
            };
            defer frame.deinit(allocator);
            try writer.writeFrame(TrrFrame{
                .step = frame.step,
                .time = frame.time,
                .lambda = 0.0,
                .box = frame.box,
                .has_x = true,
                .has_v = false,
                .has_f = false,
                .coords = frame.coords,
                .velocities = null,
                .forces = null,
            });
        }
    }
    defer std.fs.cwd().deleteFile(trr_path) catch {};

    // Step 3: Compare original TRR with TRR-via-XTC
    var orig_reader = try TrrReader.open(allocator, "test_data/frame0.trr");
    defer orig_reader.close();

    var conv_reader = try TrrReader.open(allocator, trr_path);
    defer conv_reader.close();

    const tolerance: f32 = 1.0 / precision + 0.001;
    var count: usize = 0;
    while (true) {
        var orig = orig_reader.readFrame() catch |err| {
            if (err == TrrError.EndOfFile) break;
            return err;
        };
        defer orig.deinit(allocator);

        var conv = conv_reader.readFrame() catch |err| {
            return err;
        };
        defer conv.deinit(allocator);

        try std.testing.expectEqual(orig.step, conv.step);
        try std.testing.expectApproxEqAbs(orig.time, conv.time, 0.001);

        const oc = orig.coords.?;
        const cc = conv.coords.?;
        for (oc, cc) |o, c| {
            try std.testing.expectApproxEqAbs(o, c, tolerance);
        }
        count += 1;
    }
    try std.testing.expectEqual(@as(usize, 501), count);
}
