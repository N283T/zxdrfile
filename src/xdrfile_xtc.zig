// XTC trajectory file reader
// Zig port of xdrfile library
//
// Based on the xdrfile implementation from mdtraj:
// https://github.com/mdtraj/mdtraj/tree/master/mdtraj/formats/xtc
//
// XTC is a compressed trajectory format using:
// - XDR (External Data Representation) for portable binary I/O
// - Custom 3D coordinate compression with delta encoding
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

pub const XtcError = error{
    FileNotFound,
    InvalidMagic,
    EndOfFile,
    ReadError,
    DecompressionError,
    BufferTooSmall,
    OutOfMemory,
};

/// XTC file magic number
const XTC_MAGIC: i32 = 1995;

/// Magic integers for coordinate compression (from GROMACS xdrfile)
const magicints = [_]u32{
    0,        0,        0,        0,       0,       0,       0,       0,       0,       8,
    10,       12,       16,       20,      25,      32,      40,      50,      64,      80,
    101,      128,      161,      203,     256,     322,     406,     512,     645,     812,
    1024,     1290,     1625,     2048,    2580,    3250,    4096,    5060,    6501,    8192,
    10321,    13003,    16384,    20642,   26007,   32768,   41285,   52015,   65536,   82570,
    104031,   131072,   165140,   208063,  262144,  330280,  416127,  524287,  660561,  832255,
    1048576,  1321122,  1664510,  2097152, 2642245, 3329021, 4194304, 5284491, 6658042, 8388607,
    10568983, 13316085, 16777216,
};

const FIRSTIDX: usize = 9;
const LASTIDX: usize = magicints.len;

/// XTC frame data
pub const XtcFrame = struct {
    step: i32,
    time: f32,
    box: [3][3]f32,
    coords: []f32, // flat array of x,y,z coordinates (length = natoms * 3)
    precision: f32,

    pub fn deinit(self: *XtcFrame, allocator: Allocator) void {
        allocator.free(self.coords);
    }
};

/// XTC file reader
pub const XtcReader = struct {
    file: std.fs.File,
    allocator: Allocator,
    natoms: i32,
    // buf1: integer coordinate buffer for decompression (thiscoord storage)
    // buf2: compressed bitstream with 3-element header for bit decoding state
    buf1: []i32,
    buf2: []i32,

    const Self = @This();

    pub fn open(allocator: Allocator, path: []const u8) !Self {
        const file = std.fs.cwd().openFile(path, .{}) catch {
            return XtcError.FileNotFound;
        };

        var reader = Self{
            .file = file,
            .allocator = allocator,
            .natoms = 0,
            .buf1 = &[_]i32{},
            .buf2 = &[_]i32{},
        };

        // Read first frame header to get natoms
        const magic = reader.readInt() catch return XtcError.ReadError;
        if (magic != XTC_MAGIC) {
            file.close();
            return XtcError.InvalidMagic;
        }

        reader.natoms = reader.readInt() catch return XtcError.ReadError;

        // Reset to beginning
        file.seekTo(0) catch return XtcError.ReadError;

        // Allocate decompression buffers
        const size3: usize = @intCast(reader.natoms * 3);
        reader.buf1 = allocator.alloc(i32, size3) catch return XtcError.OutOfMemory;
        // buf2: size3 * 1.2 for worst-case compression + 3 for bit decoder header
        const buf2_size: usize = size3 + size3 / 5;
        reader.buf2 = allocator.alloc(i32, buf2_size + 3) catch {
            allocator.free(reader.buf1);
            return XtcError.OutOfMemory;
        };

        return reader;
    }

    pub fn close(self: *Self) void {
        self.allocator.free(self.buf1);
        self.allocator.free(self.buf2);
        self.file.close();
    }

    pub fn getNumAtoms(self: *const Self) i32 {
        return self.natoms;
    }

    /// Read next frame from the XTC file.
    pub fn readFrame(self: *Self) !XtcFrame {
        // Read header: magic, natoms, step, time
        const magic = self.readInt() catch return XtcError.EndOfFile;
        if (magic != XTC_MAGIC) {
            return XtcError.InvalidMagic;
        }

        const natoms = try self.readInt();
        if (natoms != self.natoms) {
            return XtcError.ReadError;
        }

        const step = try self.readInt();
        const time = try self.readFloat();

        // Read box (3x3 matrix)
        var box: [3][3]f32 = undefined;
        for (0..3) |i| {
            for (0..3) |j| {
                box[i][j] = try self.readFloat();
            }
        }

        // Allocate output coordinates
        const size3: usize = @intCast(natoms * 3);
        const coords = self.allocator.alloc(f32, size3) catch return XtcError.OutOfMemory;
        errdefer self.allocator.free(coords);

        // Decompress coordinates
        var precision: f32 = 0;
        const result = self.decompressCoords(coords, &precision);
        if (result < 0) {
            return XtcError.DecompressionError;
        }

        return XtcFrame{
            .step = step,
            .time = time,
            .box = box,
            .coords = coords,
            .precision = precision,
        };
    }

    // ============================================
    // XDR I/O (big-endian, 4 bytes per element)
    // ============================================

    fn readInt(self: *Self) !i32 {
        var buf: [4]u8 = undefined;
        const n = self.file.read(&buf) catch return XtcError.ReadError;
        if (n < 4) return XtcError.EndOfFile;
        return @bitCast(std.mem.readInt(u32, &buf, .big));
    }

    fn readFloat(self: *Self) !f32 {
        var buf: [4]u8 = undefined;
        const n = self.file.read(&buf) catch return XtcError.ReadError;
        if (n < 4) return XtcError.EndOfFile;
        return @bitCast(std.mem.readInt(u32, &buf, .big));
    }

    fn readOpaque(self: *Self, dest: []u8) !void {
        const n = self.file.readAll(dest) catch return XtcError.ReadError;
        if (n < dest.len) return XtcError.EndOfFile;
        // XDR opaque data is padded to 4-byte boundary
        const padding = (4 - (dest.len % 4)) % 4;
        if (padding > 0) {
            var pad_buf: [3]u8 = undefined;
            _ = self.file.read(pad_buf[0..padding]) catch return XtcError.ReadError;
        }
    }

    // ============================================
    // Bit-level encoding/decoding utilities
    // ============================================

    /// Return the number of bits needed to represent `size` values (0..size-1).
    fn sizeofint(size: u32) u32 {
        var num: u32 = 1;
        var num_of_bits: u32 = 0;
        while (size >= num and num_of_bits < 32) {
            num_of_bits += 1;
            num <<= 1;
        }
        return num_of_bits;
    }

    /// Return the number of bits needed to represent the product of `sizes`.
    fn sizeofints(num_of_ints: usize, sizes: []const u32) u32 {
        var bytes: [32]u32 = undefined;
        var num_of_bytes: usize = 1;
        bytes[0] = 1;

        for (0..num_of_ints) |i| {
            var tmp: u32 = 0;
            for (0..num_of_bytes) |bytecnt| {
                tmp = bytes[bytecnt] * sizes[i] + tmp;
                bytes[bytecnt] = tmp & 0xff;
                tmp >>= 8;
            }
            while (tmp != 0) {
                bytes[num_of_bytes] = tmp & 0xff;
                num_of_bytes += 1;
                tmp >>= 8;
            }
        }

        var num: u32 = 1;
        var num_of_bits: u32 = 0;
        num_of_bytes -= 1;
        while (bytes[num_of_bytes] >= num) {
            num_of_bits += 1;
            num *= 2;
        }
        return num_of_bits + @as(u32, @intCast(num_of_bytes)) * 8;
    }

    /// Decode `num_of_bits` bits from the compressed buffer.
    /// buf[0] = byte counter, buf[1] = lastbits, buf[2] = lastbyte
    fn decodebits(buf: []i32, num_of_bits_arg: u32) i32 {
        var num_of_bits = num_of_bits_arg;
        const cbuf: [*]u8 = @ptrCast(@alignCast(buf.ptr + 3));

        var cnt: usize = @intCast(@as(u32, @bitCast(buf[0])));
        var lastbits: u32 = @bitCast(buf[1]);
        var lastbyte: u32 = @bitCast(buf[2]);

        const mask: u32 = if (num_of_bits_arg >= 32)
            0xFFFFFFFF
        else
            (@as(u32, 1) << @as(u5, @intCast(num_of_bits_arg))) - 1;

        var num: u32 = 0;
        while (num_of_bits >= 8) {
            lastbyte = (lastbyte << 8) | cbuf[cnt];
            cnt += 1;
            const shift_r: u5 = @intCast(lastbits & 31);
            const shift_l: u5 = @intCast((num_of_bits - 8) & 31);
            num |= (lastbyte >> shift_r) << shift_l;
            num_of_bits -= 8;
        }
        if (num_of_bits > 0) {
            if (lastbits < num_of_bits) {
                lastbits += 8;
                lastbyte = (lastbyte << 8) | cbuf[cnt];
                cnt += 1;
            }
            lastbits -= num_of_bits;
            const shift: u5 = @intCast(lastbits & 31);
            const mask_bits: u5 = @intCast(num_of_bits & 31);
            num |= (lastbyte >> shift) & ((@as(u32, 1) << mask_bits) -% 1);
        }

        num &= mask;
        buf[0] = @bitCast(@as(u32, @intCast(cnt)));
        buf[1] = @bitCast(lastbits);
        buf[2] = @bitCast(lastbyte);
        return @bitCast(num);
    }

    /// Decode multiple packed integers from the compressed buffer.
    fn decodeints(buf: []i32, num_of_ints: usize, num_of_bits_arg: u32, sizes: []const u32, nums: []i32) void {
        var bytes: [32]i32 = undefined;
        bytes[1] = 0;
        bytes[2] = 0;
        bytes[3] = 0;

        var num_of_bytes: usize = 0;
        var num_of_bits = num_of_bits_arg;
        while (num_of_bits > 8) {
            bytes[num_of_bytes] = decodebits(buf, 8);
            num_of_bytes += 1;
            num_of_bits -= 8;
        }
        if (num_of_bits > 0) {
            bytes[num_of_bytes] = decodebits(buf, num_of_bits);
            num_of_bytes += 1;
        }

        var i = num_of_ints - 1;
        while (i > 0) : (i -= 1) {
            var num: i32 = 0;
            var j = num_of_bytes;
            while (j > 0) {
                j -= 1;
                num = (num << 8) | bytes[j];
                const size_i: i32 = @intCast(sizes[i]);
                const p = @divTrunc(num, size_i);
                bytes[j] = p;
                num = num - p * size_i;
            }
            nums[i] = num;
        }
        nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
    }

    // ============================================
    // 3D coordinate decompression (XTC-specific)
    // ============================================

    /// Decompress XTC-compressed 3D coordinates.
    /// Based on xdrfile_decompress_coord_float() from mdtraj/xdrfile.
    fn decompressCoords(self: *Self, coords: []f32, precision: *f32) i32 {
        const buf1 = self.buf1;
        const buf2 = self.buf2;

        // Read number of atoms
        const lsize: i32 = self.readInt() catch return -1;
        if (lsize <= 0) return -1;
        const size3: usize = @intCast(lsize * 3);

        // Don't bother with compression for 9 atoms or less
        if (lsize <= 9) {
            for (0..size3) |i| {
                coords[i] = self.readFloat() catch return -1;
            }
            return @divTrunc(lsize, 3);
        }

        // Read precision
        precision.* = self.readFloat() catch return -1;

        buf2[0] = 0;
        buf2[1] = 0;
        buf2[2] = 0;

        // Read min/max integer bounds
        var minint: [3]i32 = undefined;
        var maxint: [3]i32 = undefined;
        for (0..3) |i| {
            minint[i] = self.readInt() catch return -1;
        }
        for (0..3) |i| {
            maxint[i] = self.readInt() catch return -1;
        }

        var sizeint: [3]u32 = undefined;
        var bitsizeint: [3]u32 = .{ 0, 0, 0 };
        sizeint[0] = @intCast(maxint[0] - minint[0] + 1);
        sizeint[1] = @intCast(maxint[1] - minint[1] + 1);
        sizeint[2] = @intCast(maxint[2] - minint[2] + 1);

        // Check if one of the sizes is too big to be multiplied
        var bitsize: u32 = 0;
        if ((sizeint[0] | sizeint[1] | sizeint[2]) > 0xffffff) {
            bitsizeint[0] = sizeofint(sizeint[0]);
            bitsizeint[1] = sizeofint(sizeint[1]);
            bitsizeint[2] = sizeofint(sizeint[2]);
            bitsize = 0; // flag the use of large sizes
        } else {
            bitsize = sizeofints(3, &sizeint);
        }

        // Read smallidx and compute bounds (mdtraj style)
        var smallidx: i32 = self.readInt() catch return -1;
        if (smallidx < FIRSTIDX or smallidx >= LASTIDX) return -1;

        var tmp: i32 = smallidx + 8;
        const maxidx: usize = if (LASTIDX < @as(usize, @intCast(tmp))) LASTIDX else @intCast(tmp);
        const minidx: usize = maxidx - 8; // often equals smallidx
        _ = minidx;

        tmp = smallidx - 1;
        if (tmp < FIRSTIDX) tmp = @intCast(FIRSTIDX);
        var smaller: i32 = @intCast(magicints[@intCast(tmp)] / 2);
        var smallnum: i32 = @intCast(magicints[@intCast(smallidx)] / 2);
        var sizesmall: [3]u32 = .{
            magicints[@intCast(smallidx)],
            magicints[@intCast(smallidx)],
            magicints[@intCast(smallidx)],
        };

        // Read compressed data length and payload
        const data_len = self.readInt() catch return -1;
        const data_len_u: usize = @intCast(data_len);

        const cbuf: [*]u8 = @ptrCast(@alignCast(buf2.ptr + 3));
        self.readOpaque(cbuf[0..data_len_u]) catch return -1;
        buf2[0] = 0;
        buf2[1] = 0;
        buf2[2] = 0;

        const inv_precision = 1.0 / precision.*;
        var run: i32 = 0;
        var i: usize = 0;
        var lfp: usize = 0;
        var prevcoord: [3]i32 = undefined;

        while (i < @as(usize, @intCast(lsize))) {
            // Use buf1 as thiscoord storage (like C version: thiscoord = (int*)(lip) + i*3)
            const thiscoord = buf1[i * 3 ..][0..3];

            if (bitsize == 0) {
                thiscoord[0] = decodebits(buf2, bitsizeint[0]);
                thiscoord[1] = decodebits(buf2, bitsizeint[1]);
                thiscoord[2] = decodebits(buf2, bitsizeint[2]);
            } else {
                decodeints(buf2, 3, bitsize, &sizeint, thiscoord);
            }

            i += 1;
            thiscoord[0] += minint[0];
            thiscoord[1] += minint[1];
            thiscoord[2] += minint[2];

            prevcoord[0] = thiscoord[0];
            prevcoord[1] = thiscoord[1];
            prevcoord[2] = thiscoord[2];

            const flag = decodebits(buf2, 1);
            var is_smaller: i32 = 0;
            if (flag == 1) {
                run = decodebits(buf2, 5);
                is_smaller = @mod(run, 3);
                run -= is_smaller;
                is_smaller -= 1;
            }

            // Buffer overrun check (mdtraj style)
            if (lfp + @as(usize, @intCast(@max(0, run))) > size3) {
                return -1;
            }

            if (run > 0) {
                // Advance thiscoord pointer into buf1 for run-length decoded atoms
                const run_start = i * 3;
                var k: i32 = 0;
                while (k < run) : (k += 3) {
                    const k_u: usize = @intCast(@divTrunc(k, 3));
                    const rc = buf1[run_start + k_u * 3 ..][0..3];
                    decodeints(buf2, 3, @intCast(smallidx), &sizesmall, rc);
                    i += 1;
                    rc[0] += prevcoord[0] - smallnum;
                    rc[1] += prevcoord[1] - smallnum;
                    rc[2] += prevcoord[2] - smallnum;
                    if (k == 0) {
                        // Interchange first with second atom for better
                        // compression of water molecules
                        const t0 = rc[0];
                        rc[0] = prevcoord[0];
                        prevcoord[0] = t0;
                        const t1 = rc[1];
                        rc[1] = prevcoord[1];
                        prevcoord[1] = t1;
                        const t2 = rc[2];
                        rc[2] = prevcoord[2];
                        prevcoord[2] = t2;

                        coords[lfp] = @as(f32, @floatFromInt(prevcoord[0])) * inv_precision;
                        coords[lfp + 1] = @as(f32, @floatFromInt(prevcoord[1])) * inv_precision;
                        coords[lfp + 2] = @as(f32, @floatFromInt(prevcoord[2])) * inv_precision;
                        lfp += 3;
                    } else {
                        prevcoord[0] = rc[0];
                        prevcoord[1] = rc[1];
                        prevcoord[2] = rc[2];
                    }

                    coords[lfp] = @as(f32, @floatFromInt(rc[0])) * inv_precision;
                    coords[lfp + 1] = @as(f32, @floatFromInt(rc[1])) * inv_precision;
                    coords[lfp + 2] = @as(f32, @floatFromInt(rc[2])) * inv_precision;
                    lfp += 3;
                }
            } else {
                coords[lfp] = @as(f32, @floatFromInt(thiscoord[0])) * inv_precision;
                coords[lfp + 1] = @as(f32, @floatFromInt(thiscoord[1])) * inv_precision;
                coords[lfp + 2] = @as(f32, @floatFromInt(thiscoord[2])) * inv_precision;
                lfp += 3;
            }

            // Adjust smallidx based on compression efficiency
            smallidx += is_smaller;
            if (is_smaller < 0) {
                smallnum = smaller;
                if (smallidx > FIRSTIDX) {
                    smaller = @intCast(magicints[@intCast(smallidx - 1)] / 2);
                } else {
                    smaller = 0;
                }
            } else if (is_smaller > 0) {
                smaller = smallnum;
                smallnum = @intCast(magicints[@intCast(smallidx)] / 2);
            }
            sizesmall[0] = magicints[@intCast(smallidx)];
            sizesmall[1] = magicints[@intCast(smallidx)];
            sizesmall[2] = magicints[@intCast(smallidx)];

            // Zero-size check (mdtraj bugfix: prevents divide-by-zero)
            if (sizesmall[0] == 0 or sizesmall[1] == 0 or sizesmall[2] == 0) {
                return -1;
            }
        }

        return lsize;
    }
};

// ============================================
// Tests
// ============================================

test "sizeofint" {
    try std.testing.expectEqual(@as(u32, 0), XtcReader.sizeofint(0));
    try std.testing.expectEqual(@as(u32, 1), XtcReader.sizeofint(1));
    try std.testing.expectEqual(@as(u32, 8), XtcReader.sizeofint(255));
    try std.testing.expectEqual(@as(u32, 9), XtcReader.sizeofint(256));
}

test "sizeofints" {
    const sizes1 = [_]u32{ 21801, 21008, 15514 };
    try std.testing.expectEqual(@as(u32, 43), XtcReader.sizeofints(3, &sizes1));

    const sizes2 = [_]u32{ 2048, 2048, 2048 };
    try std.testing.expectEqual(@as(u32, 34), XtcReader.sizeofints(3, &sizes2));
}

test "magicints table" {
    try std.testing.expectEqual(@as(u32, 8), magicints[9]);
    try std.testing.expectEqual(@as(u32, 1024), magicints[30]);
    try std.testing.expectEqual(@as(u32, 16777216), magicints[72]);
}

test "read 1l2y.xtc first frame" {
    const allocator = std.testing.allocator;

    var reader = try XtcReader.open(allocator, "test_data/1l2y.xtc");
    defer reader.close();

    const natoms = reader.getNumAtoms();
    try std.testing.expectEqual(@as(i32, 304), natoms);

    var frame = try reader.readFrame();
    defer frame.deinit(allocator);

    try std.testing.expectEqual(@as(i32, 1), frame.step);

    const tolerance: f32 = 0.0001;

    // atom[0]
    try std.testing.expectApproxEqAbs(@as(f32, -0.8901), frame.coords[0], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 0.4127), frame.coords[1], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, -0.0555), frame.coords[2], tolerance);

    // atom[1]
    try std.testing.expectApproxEqAbs(@as(f32, -0.8608), frame.coords[3], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 0.3135), frame.coords[4], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, -0.1618), frame.coords[5], tolerance);

    // atom[152]
    try std.testing.expectApproxEqAbs(@as(f32, -0.3502), frame.coords[152 * 3], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, -0.415), frame.coords[152 * 3 + 1], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, -0.4813), frame.coords[152 * 3 + 2], tolerance);

    // atom[302]
    try std.testing.expectApproxEqAbs(@as(f32, 0.1636), frame.coords[302 * 3], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 1.1959), frame.coords[302 * 3 + 1], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 0.1824), frame.coords[302 * 3 + 2], tolerance);

    // atom[303] (last atom)
    try std.testing.expectApproxEqAbs(@as(f32, 0.2831), frame.coords[303 * 3], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 1.004), frame.coords[303 * 3 + 1], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 0.2676), frame.coords[303 * 3 + 2], tolerance);
}

test "read 1l2y.xtc all frames" {
    const allocator = std.testing.allocator;

    var reader = try XtcReader.open(allocator, "test_data/1l2y.xtc");
    defer reader.close();

    var frame_count: usize = 0;
    while (true) {
        var frame = reader.readFrame() catch |err| {
            if (err == XtcError.EndOfFile) break;
            return err;
        };
        defer frame.deinit(allocator);
        frame_count += 1;

        try std.testing.expectEqual(@as(usize, 304 * 3), frame.coords.len);
    }

    try std.testing.expectEqual(@as(usize, 38), frame_count);
}

test "read large xtc (6qfk 90MB)" {
    const allocator = std.testing.allocator;

    var reader = XtcReader.open(allocator, "benchmarks/md_data/6qfk_A_analysis/6qfk_A_R1.xtc") catch |err| {
        if (err == error.FileNotFound) return;
        return err;
    };
    defer reader.close();

    const natoms: usize = @intCast(reader.getNumAtoms());
    try std.testing.expectEqual(@as(usize, 20391), natoms);

    const tolerance: f32 = 0.0001;
    var frame_count: usize = 0;

    while (true) {
        var frame = reader.readFrame() catch |err| {
            if (err == XtcError.EndOfFile) break;
            return err;
        };
        defer frame.deinit(allocator);
        frame_count += 1;

        if (frame_count == 1) {
            try std.testing.expectApproxEqAbs(@as(f32, 9.561), frame.coords[0], tolerance);
            try std.testing.expectApproxEqAbs(@as(f32, 0.786), frame.coords[1], tolerance);
            try std.testing.expectApproxEqAbs(@as(f32, 1.222), frame.coords[2], tolerance);
            try std.testing.expectApproxEqAbs(@as(f32, 9.318), frame.coords[20390 * 3], tolerance);
            try std.testing.expectApproxEqAbs(@as(f32, 5.342), frame.coords[20390 * 3 + 1], tolerance);
            try std.testing.expectApproxEqAbs(@as(f32, 7.609), frame.coords[20390 * 3 + 2], tolerance);
        } else if (frame_count == 100) {
            try std.testing.expectApproxEqAbs(@as(f32, 10.7114), frame.coords[0], tolerance);
            try std.testing.expectApproxEqAbs(@as(f32, 2.4758), frame.coords[1], tolerance);
            try std.testing.expectApproxEqAbs(@as(f32, 1.8902), frame.coords[2], tolerance);
        } else if (frame_count == 1000) {
            try std.testing.expectApproxEqAbs(@as(f32, 9.830999), frame.coords[0], tolerance);
            try std.testing.expectApproxEqAbs(@as(f32, 1.4728), frame.coords[1], tolerance);
            try std.testing.expectApproxEqAbs(@as(f32, 3.9008), frame.coords[2], tolerance);
        }
    }

    try std.testing.expectEqual(@as(usize, 1001), frame_count);
}

test "read very large xtc (5ltj 511MB)" {
    const allocator = std.testing.allocator;

    var reader = XtcReader.open(allocator, "benchmarks/md_data/5ltj_A_protein/5ltj_A_prod_R1_fit.xtc") catch |err| {
        if (err == error.FileNotFound) return;
        return err;
    };
    defer reader.close();

    const natoms: usize = @intCast(reader.getNumAtoms());
    try std.testing.expectEqual(@as(usize, 11487), natoms);

    const tolerance: f32 = 0.0001;
    var frame_count: usize = 0;

    while (true) {
        var frame = reader.readFrame() catch |err| {
            if (err == XtcError.EndOfFile) break;
            return err;
        };
        defer frame.deinit(allocator);
        frame_count += 1;

        if (frame_count == 1) {
            try std.testing.expectApproxEqAbs(@as(f32, 1.626), frame.coords[0], tolerance);
            try std.testing.expectApproxEqAbs(@as(f32, 5.291), frame.coords[1], tolerance);
            try std.testing.expectApproxEqAbs(@as(f32, 5.627), frame.coords[2], tolerance);
            try std.testing.expectApproxEqAbs(@as(f32, 8.840), frame.coords[11486 * 3], tolerance);
            try std.testing.expectApproxEqAbs(@as(f32, 2.650), frame.coords[11486 * 3 + 1], tolerance);
            try std.testing.expectApproxEqAbs(@as(f32, -0.004), frame.coords[11486 * 3 + 2], tolerance);
        } else if (frame_count == 1000) {
            try std.testing.expectApproxEqAbs(@as(f32, 0.902), frame.coords[0], tolerance);
            try std.testing.expectApproxEqAbs(@as(f32, 5.5208), frame.coords[1], tolerance);
            try std.testing.expectApproxEqAbs(@as(f32, 6.3102), frame.coords[2], tolerance);
        }
    }

    try std.testing.expectEqual(@as(usize, 10001), frame_count);
}
