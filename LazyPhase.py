#!/usr/bin/env python3


"""
Calculate radius of gyration from a PDB file.
All atoms are treated as having equal mass.
"""

"""
Count contacts between atoms of different types in a PDB file.
Type determination: the B‑factor (temperature factor) column must contain
the integer type (1, 2, or 3). All other atoms are ignored.
Radii: r1 = 2, r2 = 2, r3 = 6 (in arbitrary units, typically Å).
Contact: dist < 2 * (2^(1/6)*r_i + 2^(1/6)*r_j).
"""


import sys
import os
import math


#os.system('python generate_polymer.py && ./README_setup.sh && ./README_run.sh')


#def parse_pdb(filename: str) -> list[tuple[float, float, float]]:
#    """
#    Extract (x, y, z) coordinates from ATOM and HETATM lines.
#    Standard PDB column slicing is used (columns 31-38, 39-46, 47-54).
#    """
#    coordinates = []
#    with open(filename, 'r') as f:
#        for line in f:
#            if line.startswith(('ATOM', 'HETATM')):
#                try:
#                    x = float(line[30:38])
#                    y = float(line[38:46])
#                    z = float(line[46:54])
#                    coordinates.append((x, y, z))
#                except ValueError:
                    # Skip malformed lines
#                    continue
#    return coordinates


# Radii for each type (type 1, 2, 3)
RADII = {1: 2.0, 2: 2.0, 3: 6.0}
SIXTH_ROOT_OF_TWO = 2.0 ** (1.0 / 6.0)


def parse_pdb(filename):
    """
    Extract (x, y, z) and integer type (1,2,3) from B‑factor of ATOM/HETATM lines.
    Standard PDB format:
        x: 31-38, y: 39-46, z: 47-54, bfactor: 61-66.
    Returns list of dicts with 'coords' (tuple) and 'type' (int).
    """
    atoms = []
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            try:
                tokens = line.split()
                x, y, z = float(tokens[5]), float(tokens[6]), float(tokens[7])
                bfactor = float(line[60:66])
                atype = int(tokens[2])
                if atype not in RADII:   # only types 1,2,3
                    continue
                atoms.append({'coords': (x, y, z), 'type': atype})
            except (ValueError, IndexError):
                continue
    return atoms

def contact_threshold(type_i, type_j):
    """Contact distance threshold for two atom types."""
    r_i = RADII[type_i]
    r_j = RADII[type_j]
    return 2.0 * (SIXTH_ROOT_OF_TWO * r_i + SIXTH_ROOT_OF_TWO * r_j)

def count_contacts(atoms):
    """
    Count pairs (i<j) with different types and distance < threshold.
    """

    count = dict()

    n = len(atoms)

    for i in range(n):
        ai = atoms[i]
        ti = ai['type']
        xi, yi, zi = ai['coords']
        for j in range(i + 1, n):
            aj = atoms[j]
            tj = aj['type']
            xj, yj, zj = aj['coords']
            dx = xi - xj
            dy = yi - yj
            dz = zi - zj
            dist_sq = dx*dx + dy*dy + dz*dz
            thresh = contact_threshold(ti, tj)
            if dist_sq < thresh * thresh:
                pair = (min((ti, tj)), max((ti, tj)))
                if pair in count:
                    count[pair] += 1
                else:
                    count[pair] = 1
    return sum(count.values())


def radius_of_gyration(coords: list[tuple[float, float, float]]) -> float:
    """
    Compute the radius of gyration Rg = sqrt( (1/N) * Σ|r_i - r_com|² ).
    All atoms are equally weighted (mass = 1).
    """
    n = len(coords)
    if n == 0:
        return 0.0

    # Centre of mass (geometric centre)
    cx = sum(p[0] for p in coords) / n
    cy = sum(p[1] for p in coords) / n
    cz = sum(p[2] for p in coords) / n

    # Sum of squared distances from the centre
    sum_sq = sum((p[0] - cx)**2 + (p[1] - cy)**2 + (p[2] - cz)**2 for p in coords)

    rg = (sum_sq / n) ** 0.5
    return rg


#if __name__ == '__main__':
#    if len(sys.argv) != 2:
#        print("Usage: python radius_of_gyration.py <pdb_file>")
#        sys.exit(1)

#    pdb_path = sys.argv[1]
#    atoms = parse_pdb(pdb_path)
#    if not atoms:
#        print("No atomic coordinates found in the file.")
#        sys.exit(1)

#    rg = radius_of_gyration(atoms)
#    print(f"Radius of gyration: {rg:.3f} Å")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python count_contacts.py <pdb_file>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    atoms = parse_pdb(pdb_file)

    if not atoms:
        print("No atoms of type 1, 2, or 3 found.")
        sys.exit(1)

    contacts = count_contacts(atoms)
    print(f"Number of contacts between different types: {contacts}")
