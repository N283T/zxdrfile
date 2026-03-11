// Benchmark for XTC and TRR readers.
//
// Reads all frames from trajectory files, reports throughput,
// and optionally validates coordinates against mdtraj reference.
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
    validated: bool,
    validation_ok: bool,
};

// Reference data parsed from JSON
const SampleFrame = struct {
    time: f64,
    coords: []f64,
    box: ?[3][3]f64 = null,
};

const Reference = struct {
    natoms: i64,
    n_frames: i64,
    sample_frames: std.json.ArrayHashMap(SampleFrame),
};

fn loadReference(allocator: std.mem.Allocator, path: []const u8) !?std.json.Parsed(Reference) {
    const file = std.fs.cwd().openFile(path, .{}) catch return null;
    defer file.close();

    const data = file.readToEndAlloc(allocator, 64 * 1024 * 1024) catch return null;
    defer allocator.free(data);

    return std.json.parseFromSlice(Reference, allocator, data, .{
        .allocate = .alloc_always,
        .ignore_unknown_fields = true,
    }) catch null;
}

fn validateFrame(coords: []f32, ref_coords: []f64, tolerance: f32) bool {
    if (coords.len != ref_coords.len) return false;
    for (0..coords.len) |i| {
        const expected: f32 = @floatCast(ref_coords[i]);
        const diff = @abs(coords[i] - expected);
        if (diff > tolerance) return false;
    }
    return true;
}

fn benchmarkXtc(allocator: std.mem.Allocator, path: []const u8, name: []const u8, ref_path: ?[]const u8) !?BenchResult {
    const file_size = blk: {
        const stat = std.fs.cwd().statFile(path) catch return null;
        break :blk stat.size;
    };

    // Load reference if available
    var parsed_ref: ?std.json.Parsed(Reference) = null;
    defer if (parsed_ref) |*p| p.deinit();
    if (ref_path) |rp| {
        parsed_ref = try loadReference(allocator, rp);
    }

    var reader = XtcReader.open(allocator, path) catch return null;
    defer reader.close();

    const natoms = reader.getNumAtoms();
    var n_frames: usize = 0;
    const validated = parsed_ref != null;
    var validation_ok = true;

    var timer = try std.time.Timer.start();

    while (true) {
        var frame = reader.readFrame() catch |err| {
            if (err == XtcError.EndOfFile) break;
            return err;
        };
        defer frame.deinit(allocator);

        // Validate sample frames
        if (parsed_ref) |pr| {
            var buf: [16]u8 = undefined;
            const key = std.fmt.bufPrint(&buf, "{d}", .{n_frames}) catch unreachable;
            if (pr.value.sample_frames.map.get(key)) |ref_frame| {
                if (!validateFrame(frame.coords, ref_frame.coords, 0.001)) {
                    validation_ok = false;
                }
            }
        }

        n_frames += 1;
    }

    const elapsed_ns = timer.read();
    const elapsed_ms = @as(f64, @floatFromInt(elapsed_ns)) / 1_000_000.0;
    const file_size_mb = @as(f64, @floatFromInt(file_size)) / (1024.0 * 1024.0);
    const throughput = file_size_mb / (elapsed_ms / 1000.0);
    const fps = @as(f64, @floatFromInt(n_frames)) / (elapsed_ms / 1000.0);

    // Validate frame count
    if (parsed_ref) |pr| {
        if (n_frames != @as(usize, @intCast(pr.value.n_frames))) {
            validation_ok = false;
        }
    }

    return BenchResult{
        .name = name,
        .file_size_mb = file_size_mb,
        .natoms = natoms,
        .n_frames = n_frames,
        .elapsed_ms = elapsed_ms,
        .throughput_mbps = throughput,
        .frames_per_sec = fps,
        .validated = validated,
        .validation_ok = validation_ok,
    };
}

fn benchmarkTrr(allocator: std.mem.Allocator, path: []const u8, name: []const u8, ref_path: ?[]const u8) !?BenchResult {
    const file_size = blk: {
        const stat = std.fs.cwd().statFile(path) catch return null;
        break :blk stat.size;
    };

    var parsed_ref: ?std.json.Parsed(Reference) = null;
    defer if (parsed_ref) |*p| p.deinit();
    if (ref_path) |rp| {
        parsed_ref = try loadReference(allocator, rp);
    }

    var reader = TrrReader.open(allocator, path) catch return null;
    defer reader.close();

    const natoms = reader.getNumAtoms();
    var n_frames: usize = 0;
    const validated = parsed_ref != null;
    var validation_ok = true;

    var timer = try std.time.Timer.start();

    while (true) {
        var frame = reader.readFrame() catch |err| {
            if (err == TrrError.EndOfFile) break;
            return err;
        };
        defer frame.deinit(allocator);

        // Validate sample frames
        if (parsed_ref) |pr| {
            var buf: [16]u8 = undefined;
            const key = std.fmt.bufPrint(&buf, "{d}", .{n_frames}) catch unreachable;
            if (pr.value.sample_frames.map.get(key)) |ref_frame| {
                if (frame.coords) |coords| {
                    if (!validateFrame(coords, ref_frame.coords, 0.0001)) {
                        validation_ok = false;
                    }
                }
            }
        }

        n_frames += 1;
    }

    const elapsed_ns = timer.read();
    const elapsed_ms = @as(f64, @floatFromInt(elapsed_ns)) / 1_000_000.0;
    const file_size_mb = @as(f64, @floatFromInt(file_size)) / (1024.0 * 1024.0);
    const throughput = file_size_mb / (elapsed_ms / 1000.0);
    const fps = @as(f64, @floatFromInt(n_frames)) / (elapsed_ms / 1000.0);

    if (parsed_ref) |pr| {
        if (n_frames != @as(usize, @intCast(pr.value.n_frames))) {
            validation_ok = false;
        }
    }

    return BenchResult{
        .name = name,
        .file_size_mb = file_size_mb,
        .natoms = natoms,
        .n_frames = n_frames,
        .elapsed_ms = elapsed_ms,
        .throughput_mbps = throughput,
        .frames_per_sec = fps,
        .validated = validated,
        .validation_ok = validation_ok,
    };
}

fn printResult(r: BenchResult) void {
    const status: []const u8 = if (!r.validated) "  " else if (r.validation_ok) "OK" else "NG";
    std.debug.print("  {s:<30} {d:6} atoms  {d:5} frames  {d:7.1} MB  {d:8.1} ms  {d:7.1} MB/s  {d:8.0} fps  [{s}]\n", .{
        r.name,
        r.natoms,
        r.n_frames,
        r.file_size_mb,
        r.elapsed_ms,
        r.throughput_mbps,
        r.frames_per_sec,
        status,
    });
}

const BenchEntry = struct {
    path: []const u8,
    name: []const u8,
    ref_path: ?[]const u8 = null,
};

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    std.debug.print("\n=== zxdrfile benchmark ===\n\n", .{});

    std.debug.print("XTC:\n", .{});
    const xtc_files = [_]BenchEntry{
        .{ .path = "test_data/1l2y.xtc", .name = "1l2y (small)" },
        .{ .path = "benchmarks/md_data/3tvj_I_R1.xtc", .name = "3tvj_I (531 atoms)", .ref_path = "benchmarks/reference/3tvj_I_R1_xtc_reference.json" },
        .{ .path = "benchmarks/md_data/5wvo_C_R1.xtc", .name = "5wvo_C (3858 atoms)", .ref_path = "benchmarks/reference/5wvo_C_R1_xtc_reference.json" },
        .{ .path = "benchmarks/md_data/6sup_A_R1.xtc", .name = "6sup_A (33377 atoms)", .ref_path = "benchmarks/reference/6sup_A_R1_xtc_reference.json" },
    };
    for (xtc_files) |f| {
        if (try benchmarkXtc(allocator, f.path, f.name, f.ref_path)) |r| {
            printResult(r);
        }
    }

    std.debug.print("\nTRR:\n", .{});
    const trr_files = [_]BenchEntry{
        .{ .path = "test_data/frame0.trr", .name = "frame0 (small)" },
        .{ .path = "benchmarks/md_data/3tvj_I_R1.trr", .name = "3tvj_I (531 atoms)", .ref_path = "benchmarks/reference/3tvj_I_R1_trr_reference.json" },
        .{ .path = "benchmarks/md_data/5wvo_C_R1.trr", .name = "5wvo_C (3858 atoms)", .ref_path = "benchmarks/reference/5wvo_C_R1_trr_reference.json" },
        .{ .path = "benchmarks/md_data/6sup_A_R1.trr", .name = "6sup_A (33377 atoms)", .ref_path = "benchmarks/reference/6sup_A_R1_trr_reference.json" },
    };
    for (trr_files) |f| {
        if (try benchmarkTrr(allocator, f.path, f.name, f.ref_path)) |r| {
            printResult(r);
        }
    }

    std.debug.print("\n", .{});
}
