// Benchmark for XTC and TRR readers.
//
// Reads all frames from trajectory files and reports throughput.
// Build with: zig build bench
//
const std = @import("std");
const xdrfile = @import("xdrfile.zig");
const XtcReader = xdrfile.XtcReader;
const XtcError = xdrfile.XtcError;
const TrrReader = xdrfile.TrrReader;
const TrrError = xdrfile.TrrError;

const BenchResult = struct {
    name: []const u8,
    file_size_mb: f64,
    natoms: i32,
    n_frames: usize,
    elapsed_ms: f64,
    throughput_mbps: f64,
    frames_per_sec: f64,
};

fn benchmarkXtc(allocator: std.mem.Allocator, path: []const u8, name: []const u8) !?BenchResult {
    const file_size = blk: {
        const stat = std.fs.cwd().statFile(path) catch return null;
        break :blk stat.size;
    };

    var reader = XtcReader.open(allocator, path) catch return null;
    defer reader.close();

    const natoms = reader.getNumAtoms();
    var n_frames: usize = 0;

    var timer = try std.time.Timer.start();

    while (true) {
        var frame = reader.readFrame() catch |err| {
            if (err == XtcError.EndOfFile) break;
            return err;
        };
        defer frame.deinit(allocator);
        n_frames += 1;
    }

    const elapsed_ns = timer.read();
    const elapsed_ms = @as(f64, @floatFromInt(elapsed_ns)) / 1_000_000.0;
    const file_size_mb = @as(f64, @floatFromInt(file_size)) / (1024.0 * 1024.0);
    const throughput = file_size_mb / (elapsed_ms / 1000.0);
    const fps = @as(f64, @floatFromInt(n_frames)) / (elapsed_ms / 1000.0);

    return BenchResult{
        .name = name,
        .file_size_mb = file_size_mb,
        .natoms = natoms,
        .n_frames = n_frames,
        .elapsed_ms = elapsed_ms,
        .throughput_mbps = throughput,
        .frames_per_sec = fps,
    };
}

fn benchmarkTrr(allocator: std.mem.Allocator, path: []const u8, name: []const u8) !?BenchResult {
    const file_size = blk: {
        const stat = std.fs.cwd().statFile(path) catch return null;
        break :blk stat.size;
    };

    var reader = TrrReader.open(allocator, path) catch return null;
    defer reader.close();

    const natoms = reader.getNumAtoms();
    var n_frames: usize = 0;

    var timer = try std.time.Timer.start();

    while (true) {
        var frame = reader.readFrame() catch |err| {
            if (err == TrrError.EndOfFile) break;
            return err;
        };
        defer frame.deinit(allocator);
        n_frames += 1;
    }

    const elapsed_ns = timer.read();
    const elapsed_ms = @as(f64, @floatFromInt(elapsed_ns)) / 1_000_000.0;
    const file_size_mb = @as(f64, @floatFromInt(file_size)) / (1024.0 * 1024.0);
    const throughput = file_size_mb / (elapsed_ms / 1000.0);
    const fps = @as(f64, @floatFromInt(n_frames)) / (elapsed_ms / 1000.0);

    return BenchResult{
        .name = name,
        .file_size_mb = file_size_mb,
        .natoms = natoms,
        .n_frames = n_frames,
        .elapsed_ms = elapsed_ms,
        .throughput_mbps = throughput,
        .frames_per_sec = fps,
    };
}

fn printResult(r: BenchResult) void {
    std.debug.print("  {s:<30} {d:6} atoms  {d:5} frames  {d:7.1} MB  {d:8.1} ms  {d:7.1} MB/s  {d:8.0} fps\n", .{
        r.name,
        r.natoms,
        r.n_frames,
        r.file_size_mb,
        r.elapsed_ms,
        r.throughput_mbps,
        r.frames_per_sec,
    });
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    std.debug.print("\n=== zxdrfile benchmark ===\n\n", .{});

    // XTC benchmarks
    std.debug.print("XTC:\n", .{});
    const xtc_files = [_]struct { path: []const u8, name: []const u8 }{
        .{ .path = "test_data/1l2y.xtc", .name = "1l2y (small)" },
        .{ .path = "benchmarks/md_data/3tvj_I_R1.xtc", .name = "3tvj_I (531 atoms)" },
        .{ .path = "benchmarks/md_data/5wvo_C_R1.xtc", .name = "5wvo_C (3858 atoms)" },
        .{ .path = "benchmarks/md_data/6sup_A_R1.xtc", .name = "6sup_A (33377 atoms)" },
    };
    for (xtc_files) |f| {
        if (try benchmarkXtc(allocator, f.path, f.name)) |r| {
            printResult(r);
        }
    }

    std.debug.print("\nTRR:\n", .{});
    const trr_files = [_]struct { path: []const u8, name: []const u8 }{
        .{ .path = "test_data/frame0.trr", .name = "frame0 (small)" },
        .{ .path = "benchmarks/md_data/3tvj_I_R1.trr", .name = "3tvj_I (531 atoms)" },
        .{ .path = "benchmarks/md_data/5wvo_C_R1.trr", .name = "5wvo_C (3858 atoms)" },
        .{ .path = "benchmarks/md_data/6sup_A_R1.trr", .name = "6sup_A (33377 atoms)" },
    };
    for (trr_files) |f| {
        if (try benchmarkTrr(allocator, f.path, f.name)) |r| {
            printResult(r);
        }
    }

    std.debug.print("\n", .{});
}
