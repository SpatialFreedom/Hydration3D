# Hydration3D

## Purpose
This program compresses and decompresses 3D coordinates using a simple
compression algorithm that reduces the size of 3D coordinate data by
almost a third. It transforms 3D coordinates from a float32 format to a
packed uint64 format, allowing for efficient storage and transmission of
3D data. The six float32 values representing the original bounding box
are stored alongside the packed coordinates, enabling the reconstruction
of the original coordinates to a 21-bit accuracy.

## Background
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

## Results
* 32 bits reduces to 23 - a 28% reduction
* 64 bits reduces to 52 - a 19% reduction
* Packing three 21-bit values into a single uint64 reduces the size by 33%
    for a loss of 3 bits of resolution but a conveniently packed 64-bit
    value.

## Reddit URL - June 6, 2025
[Simple 3D Coordinate Compression - Duh! What Do You Think?](https://www.reddit.com/r/GraphicsProgramming/comments/1l4kyqm/simple_3d_coordinate_compression_duh_what_do_you)

Algorithm
=========
Dehydration (Compression):
1. Read a file containing an array of float32 x,y,z coordinates.
2. Find the minimum and maximum values of each coordinate axis.
3. Scale and translate the coordinates so that they fit into a
       [(1.75,1.75,1.75) .. (2.0,2.0,2.0)) cube - the top 11 bits of all values
       are 0b00111111111.
4. Extract each 3D coordinate's lower 21 mantissa bits and pack them into
       a uint64.
5. Write to a file the original bounding box as six float32 values followed
       by the packed uint64 values.

Rehydration (Decompression):
1. Read a file containing the original bounding box as six float32 values
       followed by the packed uint64 values.
2. Unpack the three 21-bit integers from each 3D coordinate's uint64 into
       three int32s.
3. Scale and translate these integral coordinates back to their original
       float32 bounding box values.
4. Write the resulting float32 x,y,z coordinates to a file.

## Notes
Some IEEE 754 float32 values are represented as follows:\
    0x40000000 = 2.0\
    0x3fffffff = 1.9999998807907104\
    0x3fe00000 = 1.75