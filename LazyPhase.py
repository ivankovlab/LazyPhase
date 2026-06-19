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
import math
import argparse
from random import shuffle
import numpy as np
from pathlib import Path


# Radii for each type (type 1, 2, 3)
RADII = {1: 2.0, 2: 2.0, 3: 6.0}
SIXTH_ROOT_OF_TWO = 2.0 ** (1.0 / 6.0)


def generate_substituted_file(path_in, path_out, pattern, substit):
    file_in = open(path_in, 'r')
    file_out = open(path_out, 'w')
    for line in file_in:
        file_out.write(line.replace(pattern, substit))
    file_in.close()
    file_out.close()


def get_seq_het(seq: str, w: int, alphabet: set) -> dict:
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

        # Add the current sliding window statistics for each residue from the
        # alphabet to the pool of the sequence.
        for a in alphabet:
            S[a].append(s[a])

    # Calculate resulting asymmetry score for each residue from the alphabet.
    for a in alphabet:
        ss = sum(s**2 for s in S[a])
        if ss > 0:
            S[a] = np.var(S[a]) / ss
        else:
            S[a] = 0

    return S


def generate_random_string(composition: list) -> str:
    """
    Generate a string with exactly N_S 'S' characters and S_L 'L' characters,
    arranged in a random order.
    """

    # Build the blocky string.
    symbols = ''.join([monomer[0] * monomer[1] for monomer in composition])
    symbols = list(symbols)  # Convert the string to list to allow shuffling.
    shuffle(symbols)

    return ''.join(symbols)


#def generate_seqs():
#    file_seqs = open('seqs.txt', 'w')

#    seqs = set([
#        'SL' * 32,
#        'SSLL' * 16,
#        ('S' * 4 + 'L' * 4) * 8,
#        ('S' * 8 + 'L' * 8) * 4,
#        'S' * 16 + 'L' * 16 + 'S' * 16 + 'L' * 16,
#        'S' * 32 + 'L' * 32,
#        ('S' * 3 + 'L' * 1) * 16,
#        'S' * 48 + 'L' * 16,
#        'L' * 8 + 'S' * 48 + 'L' * 8,
#        ('S' * 24 + 'L' * 8) * 2,
#        ('L' * 6 + 'S' * 20 + 'L' * 2 + 'S' * 4) * 2,
#        ('S' * 12 + 'L' * 4) * 4,
#        'L' * 8 + 'S' * 24 + ('S' * 3 + 'L' * 1) * 8
#    ])

#    for i in range(10):
#        N_S = randint(16, 32)
#        seqs.add(generate_random_string(N_S, 64 - N_S))

#    for seq in seqs:
#        file_seqs.write(seq + '\n')

#    file_seqs.close()


def parse_pdb(filename: str):
    """
    Extract (x, y, z) and integer type (1,2,3) of ATOM/HETATM lines.
    Standard PDB format.
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


def generate_sequences(args):
    dir = args.dir
    if dir[-1] != '/':
        dir += '/'

    try:
        seqs_file = open(dir + 'seqs.txt', 'r')
        seqs = seqs_file.readlines()
        seqs_file.close()
    except:
        print('Can not read seqs.txt file in directory ' + dir + '.\n' + \
              'Probably it does not exist or is unreadable.\n' + \
              'The sequences file will be created now.')

        seqs = []

    seqs = [seq[:-1] for seq in seqs]
    seqs = set(seqs)

    num_seqs_before = len(seqs)
    print(num_seqs_before, 'different sequences before generation.')

    composition = args.composition.split('.')
    composition = [(monomer[0], int(monomer[1:])) for monomer in composition]

    num_to_generate = args.num

    print('Generate sequences until', num_to_generate,
          'new sequences will be obtained...')

    while len(seqs) < num_seqs_before + num_to_generate:
        seqs.add(generate_random_string(composition))

    print(len(seqs), 'different sequences after generation.')

    seqs_file = open(dir + 'seqs.txt', 'w')
    for seq in seqs:
        seqs_file.write(seq + '\n')
    seqs_file.close()

    print('The sequences are written in ' + dir + 'seqs.txt.')

    #hets = [get_seq_het(seq, 8, {'L', 'S'}) for seq in seqs]

    #print(hets)


def prepare(args):
    print('Preparing...')

    dir = args.dir
    if dir[-1] != '/':
        dir += '/'

    substitutions = [
        (, 'DT'),
        (, 'LOG_STEPS'),
        (, 'TEMPERATURE'),
        (, 'DAMPING_TIME')
    ]

    # Prepare simulations pipeline.
    run_file = open(args.dir + 'template.run.in.nvt', 'w')
    run_file.close()

    return


def run(args):
    print('Started the high-throughput condensate simulations...')
    return


def analyse(args):
    print('Preparing report on simulations...')

    dir_path = Path(args.dir)

    if dir_path.is_dir():
        print('Report directory already exists.')
    else:
        path.mkdir(parents=False, exist_ok=False)
        print('Report directory will is created.')

    report = dict()
    for seq in seqs:
        report[seq] = dict()
        report[seq]['contacts'] = count_contacts()
        report[seq]['R_gyr'] = radius_of_gyration()
        report[seq]['clustering'] = clustering()


if __name__ == '__main__':  # If run as CLI tool.
    ###########################
    # Create the main parser. #
    ###########################

    # Welcome message of the CLI.
    parser = argparse.ArgumentParser(
                prog='LazyPhase',
                description='''Prepare, run, and analyse high-throughput
                               simulations of liquid condensates formed by
                               two-bead-per-residue polypeptide-like model
                               polymers.''',
                epilog='(C) Egor Vasilenko, 2026'
             )

    # Add commands to the CLI tool.
    subparsers = parser.add_subparsers()

    #####################################
    # Generate sequences for screening. #
    #####################################

    parser_generate_sequences = subparsers.add_parser(
        "generate_sequences",
        help="Generate sequences for the simulations."
    )

    parser_generate_sequences.add_argument(
        "--dir", type=str,
        help="Directory path for the simulations."
    )

    parser_generate_sequences.add_argument(
        "--num", type=int,
        help="Number of new sequences to generate."
    )

    parser_generate_sequences.add_argument(
        "--composition", type=str,
        help="""Composition of the sequences in format, where A20.B2.C6
                corresponds to 20 A, 2 B and 6 C beads, for example."""
    )

    parser_generate_sequences.set_defaults(func=generate_sequences)

    #########################################################################
    # Prepare the simulations for the sequences generated on previous step. #
    #########################################################################

    parser_prepare = subparsers.add_parser(
        "prepare",
        help="Prepare the simulations."
    )

    parser_prepare.add_argument(
        "--dir", type=str,
        help="Directory path for the simulations."
    )

    parser_prepare.add_argument(
        "--repeats", type=int, default=3,
        help="Number of simulation repeats for each sequence. Default=3."
    )

    parser_prepare.add_argument(
        "--solution_time", type=float, default=0.0001,
        help="Solution simulation time in microseconds. Default=0.0001."
    )

    parser_prepare.add_argument(
        "--trap_time", type=float, default=0.0002,
        help="Trap simulation time in microseconds. Default=0.0002."
    )

    parser_prepare.add_argument(
        "--relaxation_time", type=float, default=1.0,
        help="Condensate simulation time in microseconds. Default=1."
    )

    parser_prepare.add_argument(
        "--dt", type=float, default=2.0,
        help="Simulation time step in femtoseconds. Default=2."
    )

    parser_prepare.add_argument(
        "--temperature", type=float, default=310.15,
        help="Temperature for NVT conditions. Default=37°C."
    )

    parser_prepare.add_argument(
        "--row", type=int, default=3,
        help="""Number of polymers per cell row. The total number of
                polymers will be cubic power of that. Default=3."""
    )

    parser_prepare.add_argument(
        "--log_steps", type=int, default=5000,
        help="""Number of steps to make next record into the log and trajectory.
                Default=5000."""
    )

    parser_prepare.add_argument(
        "--damping_time", type=int, default=1000,
        help="""Damping time of Langevin thermostat in picoseconds.
                Default=1000."""
    )

    parser_prepare.add_argument(
        "--seed", type=int, default=None,
        help="""Seed for the thermostat."""
    )

    parser_prepare.add_argument(
        "--trap_force", type=float, default=0.1,
        help="""Force to attract to the condensate sphere."""
    )

    parser_prepare.add_argument(
        "--residue_radius_for_trap", type=float, default=0.1,
        help="""Single residue effective radius from which the trap condensate
                radius will be calculated."""
    )

    parser_prepare.set_defaults(func=prepare)

    ########################################################################
    # Run the simulations with settings defined on the previous two steps. #
    ########################################################################

    parser_run = subparsers.add_parser(
        "run",
        help="Run the simulations."
    )

    parser_run.add_argument(
        "--dir", type=str,
        help="Directory path for the simulations."
    )

    parser_run.set_defaults(func=run)

    ######################################
    # Analyse the performed simulations. #
    ######################################

    parser_analyse = subparsers.add_parser(
        "analyse",
        help="Analyse the simulations."
    )

    parser_analyse.add_argument(
        "--dir", type=str,
        help="Directory path of the simulations."
    )

    parser_analyse.add_argument(
        "--cluster", type=int,
        help="Clustering."
    )

    parser_analyse.add_argument(
        "--gyration", type=int,
        help="Radius of gyration."
    )

    parser_analyse.add_argument(
        "--contacts", type=int,
        help="Contact statistics."
    )

    parser.print_help()

    # Parse the input and follow the commands.
    args = parser.parse_args()
    args.func(args)

    #exit(0)

    #if len(sys.argv) != 2:
    #    print("Usage: python count_contacts.py <pdb_file>")
    #    sys.exit(1)

    #pdb_file = sys.argv[1]
    #atoms = parse_pdb(pdb_file)

    #if not atoms:
    #    print("No atoms of type 1, 2, or 3 found.")
    #    sys.exit(1)

    #contacts = count_contacts(atoms)
    #print(f"Number of contacts between different types: {contacts}")
