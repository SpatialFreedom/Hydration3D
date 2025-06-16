"""Hydration3D - Simple 3D Coordinate Compression and Decompression

Author: John Hilton, Spatial Freedom Pty Ltd
Date: 10 June 2025
Sydney, Australia

HISTORY
=======
10 June 2025 - John Hilton - Initial version.

PURPOSE
=======
This program compresses and decompresses 3D coordinates using a simple
compression algorithm that reduces the size of 3D coordinate data by
almost a third. It transforms 3D coordinates from a float32 format to a
packed uint64 format, allowing for efficient storage and transmission of
3D data. The six float32 values representing the original bounding box
are stored alongside the packed coordinates, enabling the reconstruction
of the original coordinates to a 21-bit accuracy.

BACKGROUND
==========
Most 3D apps use 32-bit floating point coordinates to represent 3D
points. These coordinates are often stored as a sequence of three
float32 values, one for each of the x, y and z axes. The effective
resolution of these values is 24 bits, as the IEEE 754 float32 format
uses one sign bit, eight exponent bits and 23 mantissa bits. The
resolution derives from the single sign bit and the 23 mantissa bits -
24 bits in total.

Transforming the coordinates to a 21-bit unsigned integer range allows
each 3D coordinate to be packed into a 64-bit value. A 21-bit resolution
is more than sufficient for most 3D applications apart from Computer
Aided Design apps.

Simple 3D Coordinate Compression is conceptually described by four steps:
    1. Take any set of single or double precision 3D coordinates.
    2. Find the x, y and z extents.
    3. Calculate the transformation matrix, using the extents, to translate
       and scale the whole set of coordinate into the range [1.0 .. 2.0)
       where "[1.0" is inclusive of 1.0 and "2.0)" is exclusive of 2.0.
       Store the three translation and one (or three) scale values to be used
       when reversing this transformation.
    4. All values are positive and all exponents are exactly the same so pack
       the mantissas together and throw away the sign bit and the exponents
       - voila!

Results
    32 bits reduces to 23 - a 28% reduction
    64 bits reduces to 52 - a 19% reduction

    Packing three 21-bit values into a single uint64 reduces the size by 33%
    for a loss of 3 bits of resolution but a conveniently packed 64-bit
    value.

Reddit URL - June 6, 2025:
    https://www.reddit.com/r/GraphicsProgramming/comments/1l4kyqm/
    simple_3d_coordinate_compression_duh_what_do_you/

ALGORITHM
=========
Dehydration (Compression):
    1. Read a file containing an array of float32 x,y,z coordinates.
    2. Find the minimum and maximum values of each coordinate axis.
    3. Scale and translate the coordinates so that they fit into a
       [(1.75,1.75,1.75) .. (2.0,2.0,2.0)) cube - the top 11 bits of all values
       are 0b00111111111.
    5. Extract each 3D coordinate's lower 21 mantissa bits and pack them into
       a uint64.
    6. Write to a file the original bounding box as six float32 values followed
       by the packed uint64 values.

Rehydration (Decompression):
    1. Read a file containing the original bounding box as six float32 values
       followed by the packed uint64 values.
    2. Unpack the three 21-bit integers from each 3D coordinate's uint64 into
       three int32s.
    3. Scale and translate these integral coordinates back to their original
       float32 bounding box values.
    4. Write the resulting float32 x,y,z coordinates to a file.

NOTES
=====
Some IEEE 754 float32 values are represented as follows:
    0x40000000 = 2.0
    0x3fffffff = 1.9999998807907104
    0x3fe00000 = 1.75
"""

import numpy as np


# Provide an empty utility class.
class Blank:
    pass


# Provide a utility global that can save objects for later investigation.
global g
g = Blank()

# Provide the float32 values
ALMOST_2 = np.array(0x40000000 - 1, dtype=np.uint32).view(np.float32).item()
ALMOST_PT25 = np.array(0x3E800000 - 1, dtype=np.uint32).view(np.float32).item()
MAXUINT21 = (1 << 21) - 1  # 0x001fffff


def read_coords32_file(in32_file_path):
    """
    Read a binary file containing an array of float32 x,y,z coordinates
    into a numpy array.
    """
    data = np.fromfile(in32_file_path, dtype=np.float32)
    num_coords = data.size // 3
    return data.reshape((num_coords, 3))


def write_coords32_file(out32_file_path, coords32):
    """
    Write an array of float32 x,y,z coordinates to a binary file.
    """
    coords32.tofile(out32_file_path)


def read_coords64_file(in64_file_path):
    """
    Read a binary file containing the two float32 3D coordinates specifying
    the original bounding box and an array of uint64 packed 21-bit integral
    3D coordinates.
    """
    with open(in64_file_path, "rb") as f:
        bounds = Blank()
        bounds.min_vals = np.fromfile(f, dtype=np.float32, count=3)
        bounds.max_vals = np.fromfile(f, dtype=np.float32, count=3)
        cube64 = np.fromfile(f, dtype=np.uint64)
    return bounds, cube64


def write_coords64_file(out64_file_path, bounds, cube64):
    """
    Write min_vals, max_vals and cube64, the packed uint64 3D coordinate array
    to a binary file.
    """
    with open(out64_file_path, "wb") as f:
        bounds.min_vals.tofile(f)
        bounds.max_vals.tofile(f)
        cube64.tofile(f)


def dehydrate(in32_file_path, out64_file_path):
    """
    Compress a file containing a list of float32 3D coordinates into a file
    containing the original bounding box and the 3D coordinates packed as
    21-bit unsigned integers into 64-bit values.
    """
    coords32 = read_coords32_file(in32_file_path)
    g.coords32a = coords32
    bounds, coords64 = dry(coords32)
    g.boundsa = bounds
    g.coords64a = coords64
    write_coords64_file(out64_file_path, bounds, coords64)


def rehydrate(in64_file_path, out32_file_path):
    """
    Decompress a file containing a the original bounding box corners, as
    float32s, and the list of uint64 values containing three 21-bit unsigned
    integer coordinate values.
    """
    bounds, coords64 = read_coords64_file(in64_file_path)
    g.boundsb = bounds
    g.coords64b = coords64
    coords32 = soak(bounds, coords64)
    g.coords32b = coords32
    write_coords32_file(out32_file_path, coords32)


def dry(coords32):
    """
    Nonuniformly scale and translate a list of float32 x, y, z coordinates into
    21-bit unsigned integers then pack three values to a uint64.
    Return the original bounds and the packed array.
    """
    # Find min and max of each coordinate axis
    bounds = Blank()
    bounds.min_vals = np.min(coords32, axis=0)
    bounds.max_vals = np.max(coords32, axis=0)

    # Nonuniformly scale then translate coords32 into a cube with corners
    #   (1.75,1.75,1.75) (ALMOST_2,ALMOST_2,ALMOST_2).
    # The sign will be 0, the mantissa will be 0x7f and the top two mantissa
    # bits will be 0b11. i.e. binary values of 0x3fe00000 to 0x3fffffff.
    scale_f2i = ALMOST_PT25 / (bounds.max_vals - bounds.min_vals)
    g.scale_f2i = scale_f2i
    translate_f2i = -bounds.min_vals * scale_f2i + 1.75
    g.translate_f2i = translate_f2i
    coords_i = (coords32 * scale_f2i + translate_f2i).view(np.uint32)
    g.coords_ia = coords_i

    # Overflow can produce 2.0 values - clip these then mask off lower 21 bits
    coords_i = np.clip(coords_i, 0x3FE00000, 0x3FFFFFFF) & 0x1FFFFF
    g.coords_ib = coords_i

    # Compose the packed 64-bit integers: (x << 42) | (y << 21) | z
    coords64 = (
        (coords_i[:, 0].astype(np.uint64) << 42)
        | (coords_i[:, 1].astype(np.uint64) << 21)
        | coords_i[:, 2].astype(np.uint64)
    )

    return bounds, coords64  # shape (n,)


def soak(bounds, coords64):
    """
    Unpack then nonuniformly scale and translate the three 21-bit coordinates
    from each uint64 value in coords64 back to their original bounding box.
    """
    scale_i2f = (bounds.max_vals - bounds.min_vals) * np.float32(1.0 / MAXUINT21)
    g.scale_i2f = scale_i2f
    xcoords = np.array((coords64 >> 42) & MAXUINT21, dtype=np.float32)
    ycoords = np.array((coords64 >> 21) & MAXUINT21, dtype=np.float32)
    zcoords = np.array(coords64 & MAXUINT21, dtype=np.float32)
    g.xcoords = xcoords
    g.ycoords = ycoords
    g.zcoords = zcoords
    coords_j = np.stack([xcoords, ycoords, zcoords], axis=1)
    return coords_j * scale_i2f + bounds.min_vals


def generate_random_3d_coordinates(num_points=1000, between=(-500, 500)):
    """
    Generate a NumPy array of shape (num_points, 3) with random 3D coordinates.
    """
    a = np.random.random((2,)) * (between[1] - between[0]) + between[0]
    return np.random.uniform(np.min(a), np.max(a), (num_points, 3)).astype(np.float32)


if __name__ == "__main__":

    # Specify the basename for the input and output files.
    basename = "Random1K"

    # Generate random 3D coordinates and write them to a file.
    g.rndm = generate_random_3d_coordinates()
    write_coords32_file(basename + ".f32", g.rndm)

    # Compress the 3D coordinates to a packed uint64 format.
    dehydrate(basename + ".f32", basename + ".u64")

    # Decompress the packed uint64 format back to float32 coordinates.
    rehydrate(basename + ".u64", basename + "_rehydrated.f32")
