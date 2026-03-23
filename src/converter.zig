// Simple format converter for cross-validation testing.
// Usage: converter <trr_to_xtc|xtc_to_trr> <input> <output>
//
const std = @import("std");
const xdrfile = @import("xdrfile.zig");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    if (args.len != 4) {
        std.debug.print("Usage: converter <trr_to_xtc|xtc_to_trr> <input> <output>\n", .{});
        std.process.exit(1);
    }

    const mode = args[1];
    const input_path = args[2];
    const output_path = args[3];

    if (std.mem.eql(u8, mode, "trr_to_xtc")) {
        try convertTrrToXtc(allocator, input_path, output_path);
    } else if (std.mem.eql(u8, mode, "xtc_to_trr")) {
        try convertXtcToTrr(allocator, input_path, output_path);
    } else {
        std.debug.print("Unknown mode: {s}\n", .{mode});
        std.process.exit(1);
    }
}

fn convertTrrToXtc(allocator: std.mem.Allocator, input: []const u8, output: []const u8) !void {
    var reader = try xdrfile.TrrReader.open(allocator, input);
    defer reader.close();

    var writer = try xdrfile.XtcWriter.open(allocator, output, reader.getNumAtoms(), .write);
    defer writer.close() catch {};

    var count: usize = 0;
    while (true) {
        var frame = reader.readFrame() catch |err| {
            if (err == xdrfile.TrrError.EndOfFile) break;
            return err;
        };
        defer frame.deinit(allocator);

        try writer.writeFrame(xdrfile.XtcFrame{
            .step = frame.step,
            .time = frame.time,
            .box = frame.box,
            .coords = frame.coords.?,
            .precision = 1000.0,
        });
        count += 1;
    }
    std.debug.print("Converted {d} frames: TRR → XTC\n", .{count});
}

fn convertXtcToTrr(allocator: std.mem.Allocator, input: []const u8, output: []const u8) !void {
    var reader = try xdrfile.XtcReader.open(allocator, input);
    defer reader.close();

    var writer = try xdrfile.TrrWriter.open(allocator, output, reader.getNumAtoms(), .write);
    defer writer.close() catch {};

    var count: usize = 0;
    while (true) {
        var frame = reader.readFrame() catch |err| {
            if (err == xdrfile.XtcError.EndOfFile) break;
            return err;
        };
        defer frame.deinit(allocator);

        try writer.writeFrame(xdrfile.TrrFrame{
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
        count += 1;
    }
    std.debug.print("Converted {d} frames: XTC → TRR\n", .{count});
}
