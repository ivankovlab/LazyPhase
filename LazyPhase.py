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
import argparse


# Radii for each type (type 1, 2, 3)
RADII = {1: 2.0, 2: 2.0, 3: 6.0}
SIXTH_ROOT_OF_TWO = 2.0 ** (1.0 / 6.0)


def get_seq_het(seq: str, w: int) -> dict:
    """
    Calculate asymmetry score.
    seq : str
        Sequence with one-letter codes of residues.
    w : int
        Length of the window.
    """

    s, S = dict(), dict()

    for a in alphabet:
        s[a], S[a] = 0, []

    L = len(seq)

    # Count residues in each sliding window.
    for i in range(0, L - w + 1):
        if i == 0:  # Initiate the sliding window.
            for j in range(0, w):
                s[seq[j]] += 1
        else:  # Update the sliding window.
            s[seq[i+w-1]] += 1
            s[seq[i-1]] -= 1

        # Add the current sliding window statistics to the pool of the
        # sequence.
        for a in alphabet:
            S[a].append(s[a])

    for a in alphabet:
        ss = sum(x**2 for x in S[a])
        if ss > 0:
            S[a] = np.var(S[a]) / ss
        else:
            S[a] = 0

    return S


def parse_pdb(filename: str):
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


def count_contacts(atoms: list):
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

    # Center of masses as a geometric center.
    COM_x = sum(p[0] for p in coords) / n
    COM_y = sum(p[1] for p in coords) / n
    COM_z = sum(p[2] for p in coords) / n

    # Sum of squared distances from the center.
    sum_sq = sum((p[0] - COM_x)**2 + (p[1] - COM_y)**2 + (p[2] - COM_z)**2 \
                 for p in coords)

    R_gyr = math.sqrt(sum_sq / n)

    return R_gyr


def clustering(coords: list[tuple[float, float, float]], R: float) -> float:
    """
    Compute clustering of the 3D points in spheres of radius R.
    """

    num_points = len(coords)
    count = 0

    for i in range(num_points):
        for j in range(i+1, num_points):
            a, b = coords[i], coords[j]
            r = math.sqrt((a[0] - b[0])**2 + \
                          (a[1] - b[1])**2 + \
                          (a[2] - b[2])**2)
            if r < R:
                count += 1

    return 2 * count / N / (N - 1)


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
    parser = argparse.ArgumentParser(
                prog='LazyPhase',
                description='''Prepare, run, and analyse high-throughput
                               simulations of liquid condensates formed by
                               two-bead-per-residue polypeptide-like model
                               polymers.''',
                epilog='(C) Egor Vasilenko, 2026'
             )

    subparsers = parser.add_subparsers()

    parser_prepare = subparsers.add_parser("prepare",
                                           help="Prepate the simulations.")

    parser_prepare.add_argument("--dir", type=str,
                                help="Directory path for the simulations.")

    parser_prepare.add_argument("--repeats", type=int, default=3,
                                help="""Number of simulation repeats for each
                                        sequence.""")

    parser_prepare.add_argument("--time", type=float, default=1.0,
                                help="""Condensate simulation time in
                                        microseconds.""")

    parser_generate_sequences = subparsers.add_parser("generate_sequences",
                                    help="""Generate sequences for the
                                            simulations.""")

    parser_generate_sequences.add_argument(
        "--dir", type=str,
        help="Directory path for the simulations."
    )

    parser_generate_sequences.add_argument(
        "--num", type=str,
        help="Number of new sequences to generate."
    )

    parser_generate_sequences.add_argument(
        "--composition", type=str,
        help="""Composition of the sequences in format, where A20.B2.C6
                corresponds to 20 A, 2 B and 6 C beads, for example."""
    )

    parser_run = subparsers.add_parser("run",
                                        help="Run the simulations.")

    parser_run.add_argument("--dir", type=str,
                            help="Directory path for the simulations.")

    parser_analyse = subparsers.add_parser("analyse",
                                           help="Analyse the simulations.")

    parser_analyse.add_argument("--dir", type=str,
                                help="Directory path of the simulations.")

    parser_analyse.add_argument("--cluster", type=int,
                                help="Clustering.")

    parser_analyse.add_argument("--gyration", type=int,
                                help="Radius of gyration.")

    parser_analyse.add_argument("--contacts", type=int,
                                help="Contact statistics.")

    parser.print_help()

    exit(0)
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
