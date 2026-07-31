#!/usr/bin/env python3
import sys

rows = 64
cols = 65
nz = rows + rows # identity + last col

print("%%MatrixMarket matrix coordinate integer general")
print(f"{rows} {cols} {nz}")

# Identity matrix (pivots)
for i in range(1, rows + 1):
    print(f"{i} {i} 1")

# Last column (non-pivot, all 1s)
for i in range(1, rows + 1):
    print(f"{i} {cols} 1")
