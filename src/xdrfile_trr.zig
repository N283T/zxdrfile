// TRR trajectory file reader
// Zig port of xdrfile_trr from mdtraj:
// https://github.com/mdtraj/mdtraj/tree/master/mdtraj/formats/xtc
//
// TRR is an uncompressed trajectory format storing:
// - Coordinates (x), velocities (v), and/or forces (f)
// - Box matrix, simulation step, time, and lambda
// - Supports both single and double precision
//
// Copyright (c) 2009-2014, Erik Lindahl & David van der Spoel
// Copyright (c) 2014, Robert T. McGibbon (mdtraj modifications)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
const std = @import("std");
const Allocator = std.mem.Allocator;
const native_endian = @import("builtin").cpu.arch.endian();

pub const TrrError = error{
    FileNotFound,
    InvalidMagic,
    InvalidHeader,
    EndOfFile,
    ReadError,
    OutOfMemory,
};

const TRR_MAGIC: i32 = 1993;
const DIM: usize = 3;
const VERSION_STRING = "GMX_trn_file";
const READ_BUF_SIZE = 65536;

/// TRR frame header (matches t_trnheader from GROMACS)
const TrrHeader = struct {
    is_double: bool = false,
    ir_size: i32 = 0,
    e_size: i32 = 0,
    box_size: i32 = 0,
    vir_size: i32 = 0,
    pres_size: i32 = 0,
    top_size: i32 = 0,
    sym_size: i32 = 0,
    x_size: i32 = 0,
    v_size: i32 = 0,
    f_size: i32 = 0,
    natoms: i32 = 0,
    step: i32 = 0,
    nre: i32 = 0,
    time: f32 = 0,
    lambda: f32 = 0,
};

/// TRR frame data
pub const TrrFrame = struct {
    step: i32,
    time: f32,
    lambda: f32,
    box: [3][3]f32,
    has_x: bool,
    has_v: bool,
    has_f: bool,
    coords: ?[]f32, // flat array (natoms * 3), null if not present
    velocities: ?[]f32,
    forces: ?[]f32,

    pub fn deinit(self: *TrrFrame, allocator: Allocator) void {
        if (self.coords) |c| allocator.free(c);
        if (self.velocities) |v| allocator.free(v);
        if (self.forces) |f| allocator.free(f);
    }
};

/// TRR file reader
pub const TrrReader = struct {
    file: std.fs.File,
    reader: std.fs.File.Reader,
    read_buf: *[READ_BUF_SIZE]u8,
    allocator: Allocator,
    natoms: i32,

    const Self = @This();

    pub fn open(allocator: Allocator, path: []const u8) !Self {
        const file = std.fs.cwd().openFile(path, .{}) catch {
            return TrrError.FileNotFound;
        };
        errdefer file.close();

        const read_buf = allocator.create([READ_BUF_SIZE]u8) catch return TrrError.OutOfMemory;
        errdefer allocator.destroy(read_buf);

        var self = Self{
            .file = file,
            .reader = undefined,
            .read_buf = read_buf,
            .allocator = allocator,
            .natoms = 0,
        };
        self.reader = file.reader(read_buf);

        // Read first header to get natoms
        const header = try self.readHeader();
        self.natoms = header.natoms;

        // Reset to beginning
        file.seekTo(0) catch return TrrError.ReadError;
        self.reader = file.reader(read_buf);

        return self;
    }

    pub fn close(self: *Self) void {
        self.allocator.destroy(self.read_buf);
        self.file.close();
    }

    pub fn getNumAtoms(self: *const Self) i32 {
        return self.natoms;
    }

    /// Read next frame from the TRR file.
    pub fn readFrame(self: *Self) !TrrFrame {
        const header = try self.readHeader();

        if (header.natoms != self.natoms) {
            return TrrError.ReadError;
        }

        const natoms_u: usize = @intCast(header.natoms);
        const size3 = natoms_u * DIM;

        // Read box, virial, pressure (skip virial and pressure)
        var box: [3][3]f32 = std.mem.zeroes([3][3]f32);
        if (header.box_size != 0) {
            if (header.is_double) {
                var dbox: [DIM * DIM]f64 = undefined;
                try self.readDoublesBulk(&dbox);
                for (0..DIM) |i| {
                    for (0..DIM) |j| {
                        box[i][j] = @floatCast(dbox[i * DIM + j]);
                    }
                }
            } else {
                var fbox: [DIM * DIM]f32 = undefined;
                try self.readFloatsBulk(&fbox);
                for (0..DIM) |i| {
                    for (0..DIM) |j| {
                        box[i][j] = fbox[i * DIM + j];
                    }
                }
            }
        }

        // Skip virial
        if (header.vir_size != 0) {
            try self.skipBytes(@intCast(header.vir_size));
        }

        // Skip pressure
        if (header.pres_size != 0) {
            try self.skipBytes(@intCast(header.pres_size));
        }

        // Read coordinates
        var coords: ?[]f32 = null;
        if (header.x_size != 0) {
            coords = try self.readVectors(header.is_double, size3);
        }
        errdefer if (coords) |c| self.allocator.free(c);

        // Read velocities
        var velocities: ?[]f32 = null;
        if (header.v_size != 0) {
            velocities = try self.readVectors(header.is_double, size3);
        }
        errdefer if (velocities) |v| self.allocator.free(v);

        // Read forces
        var forces: ?[]f32 = null;
        if (header.f_size != 0) {
            forces = try self.readVectors(header.is_double, size3);
        }

        return TrrFrame{
            .step = header.step,
            .time = header.time,
            .lambda = header.lambda,
            .box = box,
            .has_x = header.x_size != 0,
            .has_v = header.v_size != 0,
            .has_f = header.f_size != 0,
            .coords = coords,
            .velocities = velocities,
            .forces = forces,
        };
    }

    // ============================================
    // Internal I/O
    // ============================================

    fn readHeader(self: *Self) !TrrHeader {
        // Magic
        const magic = self.readInt() catch return TrrError.EndOfFile;
        if (magic != TRR_MAGIC) return TrrError.InvalidMagic;

        // Version string: slen (int) then XDR string (int length + padded data)
        const slen = try self.readInt();
        const expected_len: i32 = @intCast(VERSION_STRING.len + 1);
        if (slen != expected_len) return TrrError.InvalidHeader;

        // XDR string encoding: 4 bytes for string length + data padded to 4-byte boundary
        const str_len: usize = @intCast(try self.readInt());
        const padded_len = ((str_len + 3) / 4) * 4;
        try self.skipBytes(padded_len);

        // Read header fields
        var header = TrrHeader{};
        header.ir_size = try self.readInt();
        header.e_size = try self.readInt();
        header.box_size = try self.readInt();
        header.vir_size = try self.readInt();
        header.pres_size = try self.readInt();
        header.top_size = try self.readInt();
        header.sym_size = try self.readInt();
        header.x_size = try self.readInt();
        header.v_size = try self.readInt();
        header.f_size = try self.readInt();
        header.natoms = try self.readInt();

        // Determine float size (single or double precision)
        const nflsize = blk: {
            if (header.box_size != 0) {
                break :blk @divTrunc(header.box_size, @as(i32, DIM * DIM));
            } else if (header.x_size != 0) {
                break :blk @divTrunc(header.x_size, header.natoms * @as(i32, DIM));
            } else if (header.v_size != 0) {
                break :blk @divTrunc(header.v_size, header.natoms * @as(i32, DIM));
            } else if (header.f_size != 0) {
                break :blk @divTrunc(header.f_size, header.natoms * @as(i32, DIM));
            } else {
                break :blk @as(i32, 4); // default to float
            }
        };
        header.is_double = (nflsize == 8);

        // Step and nre
        header.step = try self.readInt();
        header.nre = try self.readInt();

        // Time and lambda
        if (header.is_double) {
            const td = try self.readDouble();
            header.time = @floatCast(td);
            const ld = try self.readDouble();
            header.lambda = @floatCast(ld);
        } else {
            header.time = try self.readFloat();
            header.lambda = try self.readFloat();
        }

        return header;
    }

    inline fn io(self: *Self) *std.io.Reader {
        return &self.reader.interface;
    }

    fn readInt(self: *Self) !i32 {
        const buf = self.io().takeArray(4) catch |err| return mapIoError(err);
        return @bitCast(std.mem.readInt(u32, buf, .big));
    }

    fn readFloat(self: *Self) !f32 {
        const buf = self.io().takeArray(4) catch |err| return mapIoError(err);
        return @bitCast(std.mem.readInt(u32, buf, .big));
    }

    fn readDouble(self: *Self) !f64 {
        const buf = self.io().takeArray(8) catch |err| return mapIoError(err);
        return @bitCast(std.mem.readInt(u64, buf, .big));
    }

    fn readExact(self: *Self, dest: []u8) !void {
        self.io().readSliceAll(dest) catch |err| return mapIoError(err);
    }

    /// Bulk read f32 array: read raw bytes then byte-swap in place.
    fn readFloatsBulk(self: *Self, dest: []f32) !void {
        const bytes: []u8 = @as([*]u8, @ptrCast(dest.ptr))[0 .. dest.len * 4];
        try self.readExact(bytes);
        if (native_endian != .big) {
            for (dest) |*d| {
                d.* = @bitCast(std.mem.bigToNative(u32, @bitCast(d.*)));
            }
        }
    }

    /// Bulk read f64 array: read raw bytes then byte-swap in place.
    fn readDoublesBulk(self: *Self, dest: []f64) !void {
        const bytes: []u8 = @as([*]u8, @ptrCast(dest.ptr))[0 .. dest.len * 8];
        try self.readExact(bytes);
        if (native_endian != .big) {
            for (dest) |*d| {
                d.* = @bitCast(std.mem.bigToNative(u64, @bitCast(d.*)));
            }
        }
    }

    fn skipBytes(self: *Self, count: usize) !void {
        self.io().discardAll(count) catch |err| return mapIoError(err);
    }

    fn mapIoError(err: std.io.Reader.Error) TrrError {
        return switch (err) {
            error.EndOfStream => TrrError.EndOfFile,
            error.ReadFailed => TrrError.ReadError,
        };
    }

    /// Read a vector array (coords/velocities/forces), handling float/double conversion.
    fn readVectors(self: *Self, is_double: bool, size3: usize) ![]f32 {
        const result = self.allocator.alloc(f32, size3) catch return TrrError.OutOfMemory;
        errdefer self.allocator.free(result);

        if (is_double) {
            // Read doubles in chunks to avoid huge temp allocation
            const chunk_size: usize = 1024;
            var tmp: [chunk_size]f64 = undefined;
            var i: usize = 0;
            while (i < size3) {
                const remaining = size3 - i;
                const n = if (remaining < chunk_size) remaining else chunk_size;
                try self.readDoublesBulk(tmp[0..n]);
                for (0..n) |j| {
                    result[i + j] = @floatCast(tmp[j]);
                }
                i += n;
            }
        } else {
            try self.readFloatsBulk(result);
        }
        return result;
    }
};

// ============================================
// Tests
// ============================================

test "TrrReader open non-existent file" {
    const allocator = std.testing.allocator;
    const result = TrrReader.open(allocator, "non_existent.trr");
    try std.testing.expectError(TrrError.FileNotFound, result);
}

test "read frame0.trr first frame" {
    const allocator = std.testing.allocator;

    var reader = try TrrReader.open(allocator, "test_data/frame0.trr");
    defer reader.close();

    // frame0.trr has 22 atoms
    try std.testing.expectEqual(@as(i32, 22), reader.getNumAtoms());

    var frame = try reader.readFrame();
    defer frame.deinit(allocator);

    try std.testing.expectEqual(@as(i32, 1), frame.step);
    try std.testing.expect(frame.has_x);

    const tolerance: f32 = 0.001;

    // time[0] = 500.00003
    try std.testing.expectApproxEqAbs(@as(f32, 500.0), frame.time, 0.001);
    // lambda[0] = 1.0
    try std.testing.expectApproxEqAbs(@as(f32, 1.0), frame.lambda, tolerance);

    // box = identity matrix
    try std.testing.expectApproxEqAbs(@as(f32, 1.0), frame.box[0][0], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 1.0), frame.box[1][1], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 1.0), frame.box[2][2], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 0.0), frame.box[0][1], tolerance);

    const coords = frame.coords.?;

    // atom[0]: [0.429, 1.31, 0.859]
    try std.testing.expectApproxEqAbs(@as(f32, 0.429), coords[0], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 1.31), coords[1], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 0.859), coords[2], tolerance);

    // atom[21] (last): [1.13, 0.909, 0.83]
    try std.testing.expectApproxEqAbs(@as(f32, 1.13), coords[21 * 3], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 0.909), coords[21 * 3 + 1], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 0.83), coords[21 * 3 + 2], tolerance);
}

test "read frame0.trr all frames" {
    const allocator = std.testing.allocator;

    var reader = try TrrReader.open(allocator, "test_data/frame0.trr");
    defer reader.close();

    var frame_count: usize = 0;
    while (true) {
        var frame = reader.readFrame() catch |err| {
            if (err == TrrError.EndOfFile) break;
            return err;
        };
        defer frame.deinit(allocator);
        frame_count += 1;

        try std.testing.expectEqual(@as(usize, 22 * 3), frame.coords.?.len);
    }

    // frame0.trr has 501 frames
    try std.testing.expectEqual(@as(usize, 501), frame_count);
}

test "read frame0.trr last frame" {
    const allocator = std.testing.allocator;

    var reader = try TrrReader.open(allocator, "test_data/frame0.trr");
    defer reader.close();

    const tolerance: f32 = 0.01;
    var frame_count: usize = 0;

    // Read all frames, keeping only the last one
    var last_coords: ?[]f32 = null;
    defer if (last_coords) |c| allocator.free(c);

    while (true) {
        var frame = reader.readFrame() catch |err| {
            if (err == TrrError.EndOfFile) break;
            return err;
        };
        frame_count += 1;

        // Keep last frame's coords, free previous
        if (last_coords) |c| allocator.free(c);
        last_coords = frame.coords;
        frame.coords = null; // prevent deinit from freeing
        if (frame.velocities) |v| allocator.free(v);
        if (frame.forces) |f| allocator.free(f);
    }

    try std.testing.expectEqual(@as(usize, 501), frame_count);

    const coords = last_coords.?;

    // Frame 501 (index 500): atom[0]: [0.77, 1.01, 0.47]
    try std.testing.expectApproxEqAbs(@as(f32, 0.77), coords[0], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 1.01), coords[1], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 0.47), coords[2], tolerance);

    // atom[21]: [0.81, 1.4, 1.11]
    try std.testing.expectApproxEqAbs(@as(f32, 0.81), coords[21 * 3], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 1.4), coords[21 * 3 + 1], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 1.11), coords[21 * 3 + 2], tolerance);
}
