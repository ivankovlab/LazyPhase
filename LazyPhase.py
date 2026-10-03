#!/usr/bin/env python3


"""
Calculate radius of gyration from a PDB file.
All atoms are treated as having equal mass.

Count contacts between atoms of different types in a PDB file.
Type determination: the B‑factor (temperature factor) column must contain
the integer type (1, 2, or 3). All other atoms are ignored.
Radii: r1 = 2, r2 = 2, r3 = 6 (in arbitrary units, typically Å).
Contact: dist_min (zero-crossing) < dist < 2 * (2^(1/6)*r_i + 2^(1/6)*r_j).
"""


import math
import argparse
from random import shuffle
import numpy as np
from scipy.stats import linregress
from pathlib import Path
import os
from matplotlib import pyplot as plot
import matplotlib
import scipy.spatial as spatial
import ast

matplotlib.rcParams['figure.dpi'] = 1000
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rc('font', family='STIXGeneral')
matplotlib.rc('font', weight='ultralight')


# Pre-calculate the sixth root of two for dipole-dipole attractions, also known
# as dispersion forces.
SIXTH_ROOT_OF_TWO = 2.0 ** (1.0 / 6.0)


def get_condensate_radius(N: int, c: int, b: float) -> float:
    """Compute the radius of the trapped condensate.

    :param N: Number of monomers in a single polymer.
    :type N: int
    :param c: Number of polymers along one edge of the cubic lattice;
        the system therefore contains :math:`c^3` polymers in total.
    :type c: int
    :param b: Side length of a cube that approximates the size of
        an individual monomer.
    :type b: float
    :return: Radius of the trapped condensate.
    :rtype: float
    """

    return (3 / 4 / math.pi * N)**(1/3) * c * b


def generate_substituted_file(path_in: str, path_out: str,
                              substitutions: list):
    """Substitute all patterns in the template file and write the resulting
       specifications file.

    :param path_in: Path to the input template file.
    :type path_in: str
    :param path_out: Path for the output specifications file.
    :type path_out: str
    :param substitutions: A list of substitution pairs, each expressed
        as a tuple of the form ``('pattern', 'word')``.
    :type substitutions: list
    :return: ``None``.
    :rtype: NoneType
    """

    file_in = open(path_in, 'r')
    file_out = open(path_out, 'w')

    for line in file_in:
        line_specified = line
        for pair in substitutions:
            pattern, substitution = pair[0], pair[1]
            line_specified = line_specified.replace(pattern, substitution)
        file_out.write(line_specified)

    file_in.close()
    file_out.close()


import run_condensates_throughput, continue_condensates_throughput


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

    composition : list
        List of tuples ('one-code residue name', # occurences in the sequence).
    """

    # Build the blocky string.
    symbols = ''.join([monomer[0] * monomer[1] for monomer in composition])
    symbols = list(symbols)  # Convert the string to list to allow shuffling.
    shuffle(symbols)

    return ''.join(symbols)


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
                if atype not in RADII:  # only types 1,2,3
                    continue
                atoms.append({'coords': (x, y, z), 'type': atype})
            except (ValueError, IndexError):
                continue

    return atoms


def contact_threshold(type_i: int, type_j: int):
    """Contact distance threshold for two atom types.
    type_i : int
        The first bead type.
    type_j : int
        The second bead type.
    """

    r_i = RADII[type_i]
    r_j = RADII[type_j]

    return 2.0 * (SIXTH_ROOT_OF_TWO * r_i + SIXTH_ROOT_OF_TWO * r_j)


def count_contacts(atoms: list):
    """
    Count pairs (i<j) with different types and distance < threshold.
    atoms : list
        List of all atoms described as dictionaries with keys 'type' and
        'coords'.
    """

    count = dict()

    n = len(atoms)

    for i in range(n):  # Read the first atom.
        ai = atoms[i]

        # Read its type and coordinates.
        ti = ai['type']
        xi, yi, zi = ai['coords']

        for j in range(i + 1, n):  # Read the second atom, the pairs count once.
            aj = atoms[j]

            # Read its type and coordinates.
            tj = aj['type']
            xj, yj, zj = aj['coords']

            # Calculate the distance between the beads.
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
    """Compute the radius of gyration Rg = sqrt( (1/N) * Σ|r_i - r_com|² ). All
    atoms are equally weighted (mass = 1).

    coords : list[tuple[float, float, float]]
        Coordinates of the bead centers.
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


def clustering(coords: list[tuple[int, int, float, float, float]], R: float,
               M: int, bead_types: list[int]) -> float:
    """Compute clustering of the 3D points in spheres of radius R.

    coords : list[tuple[int, int, float, float, float]]
        Molecule IDs, bead types, and coordinates of the bead centers.
    R : float
        Radius of the clustering to search for the adjacent beads.
    M : int
        Number of polymers in the system.
    bead_types : list[int]
        Bead types to calculate clustering for.
    """

    num_points = len(coords)

    # Build a KD-tree for the points.
    tree = spatial.cKDTree([(bead[2], bead[3], bead[4]) for bead in coords])

    # Counter of the intermolecular contacts.
    count = 0

    # Number of appropriate beads in a polymer.
    L = 0

    for j in range(num_points):
        if coords[j][1] in bead_types and coords[j][0] == 1:
            L += 1
        else:
            continue

        point = (coords[j][2], coords[j][3], coords[j][4])

        indices = tree.query_ball_point(point, R)

        for i in indices:
            if coords[i][0] != coords[j][0] and coords[i][1] in bead_types:
                count += 1

    # Maximal number of intermolecular contacts.
    max_inter_contacts = M * L * (M - 1) * L / 2

    # Calculate and return the clustering coefficient.
    return count / 2 / max_inter_contacts


def generate_sequences(args):
    """Generate sequences in seqs.txt in the target directory of the
       simulations.

    args
        CLI arguments of the command generate_sequences.
    """

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

    # Generate the set of different sequences.
    if args.blocky is None:  # Generate random sequence of desired composition.
        composition = args.composition.split('.')
        composition = [
            (monomer[0], int(monomer[1:])) for monomer in composition
        ]

        num_to_generate = args.num

        print('Generate sequences until', num_to_generate,
              'new sequences will be obtained...')

        while len(seqs) < num_seqs_before + num_to_generate:
            seqs.add(generate_random_string(composition))
    else:  # Generate all blocky sequences.
        N = args.blocky

        if N <= 0 or (N & (N - 1)) != 0:
            print('The blocky polymer length must be a power of two.')
            return

        block = 1
        while block <= N // 2:
            seqs.add(('S' * block + 'L' * block) * (N // (block * 2)))
            block *= 2

    print(len(seqs), 'different sequences after generation.')

    seqs_file = open(dir + 'seqs.txt', 'w')
    for seq in seqs:
        seqs_file.write(seq + '\n')
    seqs_file.close()

    print('The sequences are written in ' + dir + 'seqs.txt.')

    #hets = [get_seq_het(seq, 8, {'L', 'S'}) for seq in seqs]

    #print(hets)


def prepare(args):
    """Prepare the target run settings in the local template.run.in.nvt file."""

    print('Preparing...')

    dir = args.dir
    if dir[-1] != '/':
        dir += '/'

    # Read sequences of the polypeptides that will be simulated.
    seqs_file = open(dir + 'seqs.txt', 'r')
    seqs = seqs_file.readlines()
    seqs_file.close()

    # Collect all lengths of the sequences. Verify that the length is the same.
    seqs_lens = set([len(seq[:-1]) for seq in seqs])
    print(len(seqs), 'sequences in ' + dir + 'seqs.txt.')

    if len(seqs_lens) > 1:
        print('The sequences must be of the same length. Terminate.')
        return

    # Identify the common length of all sequences.
    N = seqs_lens.pop()
    print('All sequences are', N, 'residues long.')

    # Convert all arguments to the appropriate units.
    dt                      = str(args.dt)
    repeats                 = str(args.repeats)
    solution_steps          = str(round(args.solution_time * 10**9 / args.dt))
    trap_steps              = str(round(args.trap_time * 10**9 / args.dt))
    relaxation_steps        = str(round(args.relaxation_time * 10**9 / args.dt))
    temperature             = str(args.temperature)
    row                     = args.row
    log_steps               = str(args.log_steps)
    damping_time            = str(args.damping_time)
    seed                    = str(args.seed)
    trap_force              = str(args.trap_force)
    residue_radius_for_trap = args.residue_radius_for_trap

    # Length of a delineated polymer along the main backbone axis.
    polymer_length = (N - 1) * 3.2

    # Define the cell boundaries for the simulation series.
    cell_file = open(dir + 'cell.lt', 'w')
    cell_file.write('write_once("Data Boundary") {\n')
    # The subcells for the polymers are sized with paddings of roughly two
    # monomer lengths.
    cell_file.write('    0 ' + str(round((polymer_length + 6.4) * row, 1)) + \
                    ' xlo xhi\n')
    cell_file.write('    0 ' + str(round((polymer_length + 6.4) * row, 1)) + \
                    ' ylo yhi\n')
    cell_file.write('    0 ' + str(round((polymer_length + 6.4) * row, 1)) + \
                    ' zlo zhi\n')
    cell_file.write('}\n\n')
    cell_file.close()

    # Specify all parameters in the run settings file that are specific for the
    # serie.
    substitutions = [
        ('DT', dt),
        ('LOG_STEPS', log_steps),
        ('TEMPERATURE', temperature),
        ('DAMPING_TIME', damping_time),
        ('SOLUTION_TIME', solution_steps),
        ('TRAP_FORCE', trap_force),
        ('TRAP_TIME', trap_steps),
        ('RELAXATION_TIME', relaxation_steps),
        (
        'RADIUS',
        str(round(get_condensate_radius(N, row, residue_radius_for_trap), 1))
        ),
        ('HALF_CELL_SIZE', str(round((polymer_length + 6.4) * row / 2, 1)))
    ]

    # Prepare simulations pipeline.
    local_template_path = dir + 'template.run.in.nvt'

    # Prepare the template of running settings.
    generate_substituted_file('template.run.in.nvt', local_template_path,
                              substitutions)

    # Create a directory for the outputs.
    os.system('mkdir ' + dir + 'output')


def run(args):
    """Run all simulations defined on the sequence generation and preparation
       steps."""

    print('Started the high-throughput condensate simulations...')

    dir = args.dir
    repeats = args.repeats
    lmp = args.lmp

    run_condensates_throughput.run_condensates_throughput(dir, lmp, repeats)

    print('The simulations finished.')


def cmd_continue(args):
    """Continue all simulations."""

    print('Continued the high-throughput condensate simulations...')

    dir = args.dir
    if dir[-1] != '/':
        dir += '/'

    template_file = open(dir + 'template.run.in.nvt', 'r')
    template_lines = template_file.readlines()
    template_file.close()

    continue_file = open(dir + 'continue.run.in.nvt', 'w')

    for line in template_lines:
        if len(line) > 3 and '#' not in line and 'run' not in line \
                         and 'write_data' not in line \
                         and 'fix trap' not in line:
             continue_file.write(
                line.replace(
                    'SEQUENCE', 'SEQUENCE_prolong').replace(
                        'system.data', 'output/SEQUENCE.data') + '\n')

        if 'timestep' in line:
            steps = round(args.time * 10**9 / float(line.split()[-1]))

    continue_file.write('''# Continue relaxation to equilibrate the pre-formed
                             condensate.\n''')
    continue_file.write('run ' + str(steps) + '\n\n')
    continue_file.write('write_data output/SEQUENCE_prolong.data\n')
    continue_file.close()

    repeats = args.repeats
    lmp = args.lmp
    interrupted = args.interrupted

    continue_condensates_throughput.continue_condensates_throughput(
        dir, lmp, repeats, interrupted
    )

    print('The simulations finished.')


def calculate_trajectory_gyrations(traj_path: str):
    """Read coordinates of all beads for all saved trajectory states."""

    traj_file = open(traj_path, 'r')
    traj_lines = traj_file.readlines()
    traj_file.close()

    states = []

    for line in traj_lines:
        if 'ITEM: ATOMS id mol type x y z ix iy iz' in line:
            states.append([])
        if '0 0 0' in line:
            tokens = line.split()
            states[-1].append((float(tokens[3]),
                              float(tokens[4]),
                              float(tokens[5])))

    gyrations = [radius_of_gyration(state) for state in states]

    return gyrations


def calculate_contacts(beads: list[tuple[int, float, float, float]],
                       bead_types: set, bead_borders: dict) -> dict:
    """Count contacts for all type pairs.

    beads : list[tuple[int, float, float, float]]
        Types and coordinates of the beads.
    bead_types : set[int]
        All identified bead types in the system.
    bead_borders : dict
        Lower and upper borders for contact distances.
    """

    contacts = dict()
    for bead_pair in bead_borders:
        contacts[bead_pair] = 0

    tree = spatial.cKDTree([(bead[1], bead[2], bead[3]) for bead in beads])

    for bead in beads:
        indices = tree.query_ball_point((bead[1], bead[2], bead[3]), 12.0)
        adjacent_points = [beads[i] for i in indices]

        for point in adjacent_points:
            type_i, type_j = bead[0], point[0]

            lower_border, upper_border = bead_borders[(min(type_i, type_j),
                                                       max(type_i, type_j))]

            distance = (bead[1] - point[1])**2 + \
                       (bead[2] - point[2])**2 + \
                       (bead[3] - point[3])**2

            if lower_border**2 < distance < upper_border**2:
                contacts[(min(type_i, type_j), max(type_i, type_j))] += 1

    for pair in contacts:
        contacts[pair] //= 2

    return contacts


def calculate_trajectory_contacts(traj_path: str) -> list:
    """Calculate contact statistics for all saved trajactory states."""

    traj_file = open(traj_path, 'r')
    traj_lines = traj_file.readlines()
    traj_file.close()

    states = []

    bead_types = set()

    for line in traj_lines:
        if 'ITEM: ATOMS id mol type x y z ix iy iz' in line:
            states.append([])
        if '0 0 0' in line:
            tokens = line.split()
            states[-1].append((int(tokens[2]),
                               float(tokens[3]),
                               float(tokens[4]),
                               float(tokens[5])))
            bead_types.add(int(tokens[2]))

    print('Identified bead types ' + str(bead_types) + '.')

    bead_borders = {
        (1, 1): (4.0, 4.0),
        (1, 2): (3.0, 3.0),
        (1, 3): (4.0, 5.0),
        (2, 2): (2.0, 3.0),
        (2, 3): (3.0, 4.0),
        (3, 3): (4.0, 12.0)
    }

    print('Contacts for each of', len(states), 'states will be counted.')

    contacts = []

    num_state = 0
    for state in states:
        contacts.append(calculate_contacts(state, bead_types, bead_borders))
        print('Contacts for state', num_state, 'counted.')
        num_state += 1

    return contacts


def get_diffusivity(timestamps: list[float],
                    gyrations_concatenated: list[float]) -> float:
    """Calculate diffusivity from three-dimensional diffusion equation
       R^2 = 6Dt."""

    D = linregress(np.array(timestamps),
                   np.array(gyrations_concatenated)**2 / 6)

    return D.slope, D.intercept


def analyse(args):
    """Calculate statistics for all available simulations.

    args
        CLI arguments.

    Returns the report for all sequences as a dictionary of dictionaries with
    the sequences as main keys and all calculated statistics as subkeys.
    """

    contacts_flag = args.contacts
    gyration_flag = args.gyration
    cluster_flag = args.cluster
    calculate_flag = args.calculate

    block_length_colors = {
        1: '#669bbc',
        2: '#003049',
        4: '#588157',
        8: '#344e41',
        16: '#dd2d4a',
        32: '#880d1e'
    }

    dir = 'examples/stickers_spacers/'

    # Calculate clustering scales, if asked.
    if cluster_flag:
        block_slopes_dict = dict()

        for bead_types in [[1, 2, 3], [1], [2], [3]]:
            bead_types_str = ''.join([str(bead_type) \
                                      for bead_type in bead_types])

            print('@@@@@@@@')
            print('Bead types ' \
                  + ', '.join([str(bead_type) for bead_type in bead_types]) \
                  + '.')
            print('@@@@@@@@')

            block_slopes = []

            for block_length in [1, 2, 4, 8, 16, 32]:
                print('#### Block length ' + str(block_length) + '. ####')

                clusts_pool = []

                for n in range(1, 10+1):
                    print('Simulation number ' + str(n) + '.')

                    if calculate_flag:
                        cluster_file = open(dir + 'statistics/cluster_' \
                                                + bead_types_str + '_' \
                                                + ('S' * block_length \
                                                + 'L' * block_length) \
                                                * (64 // (block_length * 2)) \
                                                + '_' + str(n) + '.txt', 'w')

                        # Read all beads from a .data file to calculate their
                        # clustering on different scales.
                        beads = read_beads_from_LAMMPS_data(
                            dir + 'output/' \
                                + ('S' * block_length + 'L' * block_length) \
                                * (64 // (block_length * 2)) + '_' + str(n) \
                                + '.data'
                        )

                        Rs, clusts = np.logspace(0, 2, 10), []

                        for R in Rs:
                            print('Radius ' + str(R) + ' Å.')
                            c = clustering(beads, R, 27, bead_types)
                            clusts.append(c)
                            cluster_file.write(str(R) + ' ' + str(c) + '\n')

                        cluster_file.close()
                    else:
                        cluster_file = open(
                            dir + 'statistics/cluster_' \
                                + bead_types_str + '_' \
                                + ('S' * block_length + 'L' * block_length) \
                                * (64 // (block_length * 2)) + '_' + str(n) \
                                + '.txt', 'r')
                        clustering_lines = cluster_file.readlines()
                        cluster_file.close()

                        Rs, clusts = [], []
                        for line in clustering_lines:
                            tokens = line.strip().split()
                            Rs.append(float(tokens[0]))
                            clusts.append(float(tokens[1]))

                    clusts_pool.append(clusts)

                # Calculate and plot the average clusterings.

                clusts_averaged, std_lower, std_upper = [], [], []

                for i in range(len(clusts)):
                    sum_at_distance = sum(
                        [clusts_pool[n][i] for n in range(10)]
                    )
                    std = np.std([clusts_pool[n][i] for n in range(10)])
                    clusts_averaged.append(sum_at_distance / 10)
                    std_lower.append(sum_at_distance / 10 - std)
                    std_upper.append(sum_at_distance / 10 + std)

                min_idx = 0
                while clusts_averaged[min_idx] == 0:
                    min_idx += 1

                plot.plot(Rs[min_idx:], clusts_averaged[min_idx:],
                          color=block_length_colors[block_length], alpha=0.5,
                          zorder=-1)

                clust_slope = linregress(
                    np.log(Rs[min_idx+1:]), np.log(clusts_averaged[min_idx+1:])
                ).slope

                block_slopes.append(clust_slope)

                #plot.fill_between(Rs, std_lower, std_upper,
                #                  color=block_length_colors[block_length],
                #                  alpha=0.25)

                plot.scatter(
                    Rs[min_idx:], clusts_averaged[min_idx:],
                    color=block_length_colors[block_length], alpha=1.0,
                    label=str(block_length)
                )

            plot.xscale('log')
            plot.yscale('log')
            plot.xlabel('R (Å)', fontsize=16)
            plot.ylabel('Clustering coefficient', fontsize=16)
            plot.legend(title='Block length', fontsize=16, title_fontsize=16)
            plot.tight_layout()
            plot.savefig('clustering_' + bead_types_str + '.png')
            plot.clf()

            block_slopes_dict[tuple(bead_types)] = block_slopes

        bead_types_names = {
            (1,): 'Backbone',
            (2,): 'Small',
            (3,): 'Large',
            (1,2,3,): 'All'
        }

        bead_types_colors = {
            (1,): '#9e2a2b',
            (2,): '#335c67',
            (3,): '#e09f3e',
            (1,2,3,): '#001219'
        }

        plot.xscale('linear')
        plot.yscale('linear')

        for bead_types in block_slopes_dict:
            block_slopes = block_slopes_dict[bead_types]

            plot.plot(
                [1, 2, 4, 8, 16, 32], block_slopes,
                color=bead_types_colors[bead_types],
                zorder=sum(bead_types)
            )

            plot.scatter(
                [1, 2, 4, 8, 16, 32], block_slopes,
                label=bead_types_names[bead_types],
                color=bead_types_colors[bead_types],
                zorder=sum(bead_types)
            )

        plot.xticks([1, 2, 4, 8, 16, 32], fontsize=12)
        plot.yticks(fontsize=12)
        plot.xlabel('Block length', fontsize=16)
        plot.ylabel('Slope', fontsize=16)
        plot.legend(fontsize=16)
        plot.tight_layout()
        plot.savefig('clustering_slopes.png')
        plot.clf()

    # Calculate radii of gyration, if asked.
    if gyration_flag:
        if calculate_flag:
            print('Calculating the radii of gyration...')

            for block_length in [1, 2, 4, 8, 16, 32]:
                print('---- Block length ' + str(block_length) + '. ----')

                for n in range(1, 10+1):
                    print('n = ' + str(n))

                    seq = ('S' * block_length + 'L' * block_length) \
                          * (64 // (block_length * 2))

                    statistics_report_file = open(
                        dir + 'statistics/gyration_' + seq + '_' + str(n) \
                            + '.txt', 'w'
                    )

                    gyrations = calculate_trajectory_gyrations(
                        dir + 'output/' + seq + '_' + str(n) + '.lammpstrj'
                    )

                    gyrations_prolong = calculate_trajectory_gyrations(
                        dir + 'output/' + seq + '_' + str(n) \
                            + '_prolong.lammpstrj'
                    )

                    gyrations_concatenated = (gyrations + gyrations_prolong[1:])

                    for gyration in gyrations_concatenated:
                        statistics_report_file.write(str(gyration) + '\n')

                    statistics_report_file.close()

        timestamps = np.array([i for i in range(301 + 1001 - 1)],
                              dtype=float) / 10.0

        print('Plotting the radii of gyration...')

        for block_length in [1, 2, 4, 8, 16, 32]:
            print('---- Block length ' + str(block_length) + '. ----')

            D_1_list, D_bias_list = [], []

            for n in range(1, 10+1):
                print('n = ' + str(n))

                seq_1 = ('S' * 1 + 'L' * 1) * (64 // (1 * 2))

                statistics_report_file = open(
                    dir + 'statistics/gyration_' + seq_1 + '_' + str(n) \
                        + '.txt', 'r'
                )
                statistics_lines = statistics_report_file.readlines()
                statistics_report_file.close()

                gyrations_concatenated = [
                    float(line.strip()) for line in statistics_lines
                ]

                D_1, D_bias = get_diffusivity(timestamps[200:],
                                              gyrations_concatenated[200:])

                D_1_list.append(D_1)
                D_bias_list.append(D_bias)

                if n == 1:
                    label = '1'
                else:
                    label = None

                #plot.scatter(timestamps[200:], gyrations[200:], label=label,
                #             color='#81b29a', alpha=0.5)
                plot.plot(timestamps[200:], gyrations_concatenated[200:],
                          color='#81b29a', alpha=0.5, label=label, zorder=0)

            plot.plot(
                timestamps[200:],
                np.sqrt(6 * (np.mean(D_bias_list) + np.mean(D_1_list) \
                          * np.array(timestamps[200:]))),
                color='#31572c', linestyle='-', zorder=1, linewidth=2
            )

            D_list, D_bias_list = [], []

            for n in range(1, 10+1):
                seq = ('S' * block_length + 'L' * block_length) \
                      * (64 // (block_length * 2))

                statistics_report_file = open(
                    dir + 'statistics/gyration_' + seq + '_' + str(n) + '.txt',
                    'r'
                )
                statistics_lines = statistics_report_file.readlines()
                statistics_report_file.close()

                gyrations_concatenated = [
                    float(line.strip()) for line in statistics_lines
                ]

                D, D_bias = get_diffusivity(timestamps[200:],
                                            gyrations_concatenated[200:])

                D_list.append(D)
                D_bias_list.append(D_bias)

                if n == 1:
                    label = str(block_length)
                else:
                    label = None

                #plot.scatter(timestamps[200:], gyrations[200:], label=label,
                #             color='#e07a5f', alpha=0.5)
                plot.plot(timestamps[200:], gyrations_concatenated[200:],
                          color='#e07a5f', alpha=0.5, label=label, zorder=0)

            plot.plot(
                timestamps[200:],
                np.sqrt(6 * (np.mean(D_bias_list) + np.mean(D_list) \
                          * np.array(timestamps[200:]))),
                color='#bc4749', linestyle='-', zorder=1, linewidth=2
            )

            plot.xlabel('Time (ns)', fontsize=16)
            plot.ylabel('Radius of gyration (Å)', fontsize=16)
            plot.legend(title='Block length', fontsize=16, title_fontsize=16)
            plot.tight_layout()
            plot.savefig('gyration_' + str(seq) + '.png')
            plot.clf()

    # Calculate contact statistics, if asked.
    if contacts_flag:
        #if calculate_flag:
            #for n in range(1, 10+1):
            #    trajectory_contacts_nucleus = calculate_trajectory_contacts('examples/stickers_spacers/output/' + ('S' * 32 + 'L' * 32) * 1 + '_' + str(n) + '.lammpstrj')
            #    trajectory_contacts_prolong = calculate_trajectory_contacts('examples/stickers_spacers/output/' + ('S' * 32 + 'L' * 32) * 1 + '_' + str(n) + '_prolong.lammpstrj')[1:]
            #    trajectory_contacts = trajectory_contacts_nucleus + trajectory_contacts_prolong

            #    statistics_report_file = open('examples/stickers_spacers/statistics/contacts_' + ('S' * 32 + 'L' * 32) * 1 + '_' + str(n) + '.txt', 'w')

            #    for contacts in trajectory_contacts:
            #        statistics_report_file.write(str(contacts) + '\n')

            #    statistics_report_file.close()

        print('Plotting contacts...')

        timestamps = np.array([i for i in range(301 + 1001 - 1)],
                              dtype=float) / 10.0

        bead_type_name_dict = {
            1: 'B',
            2: 'S',
            3: 'L'
        }

        for bead_type_1 in range(1, 3+1):
            for bead_type_2 in range(bead_type_1, 3+1):
                bead_type_name_1 = bead_type_name_dict[bead_type_1]
                bead_type_name_2 = bead_type_name_dict[bead_type_2]

                print('### Bead types ' + bead_type_name_1 + ' and ' \
                                        + bead_type_name_2 + '. ###')

                for block_length in [1, 2, 4, 8, 16, 32]:
                    print('---- Block length ' + str(block_length) + '. ----')

                    seq = ('S' * block_length + 'L' * block_length) \
                          * (64 // (2 * block_length))

                    block_statistics_states = []

                    label = str(block_length)

                    for n in range(1, 10+1):
                        print('n = ' + str(n))

                        statistics_report_file = open(
                            dir + 'statistics/contacts_' + seq + '_' + str(n) \
                                + '.txt', 'r'
                        )
                        statistics_lines = statistics_report_file.readlines()
                        statistics_report_file.close()

                        statistics_states = []
                        for line in statistics_lines[301:]:
                            contacts_dict = ast.literal_eval(line)
                            statistics_states.append(
                                contacts_dict[(bead_type_1, bead_type_2)]
                            )

                        #if n == 1:
                        #    label = str(block_length)
                        #else:
                        #    label = None

                        block_statistics_states.append(statistics_states)

                    maxima = [max(col) for col in zip(*block_statistics_states)]
                    minima = [min(col) for col in zip(*block_statistics_states)]
                    means  = [sum(col) / len(col) \
                              for col in zip(*block_statistics_states)]

                    plot.plot(timestamps[301:], maxima,
                              label=label,
                              color=block_length_colors[block_length],
                              alpha=0.5)

                    plot.plot(timestamps[301:], minima,
                              color=block_length_colors[block_length],
                              alpha=0.5)

                    plot.plot(timestamps[301:], means,
                              color=block_length_colors[block_length])

                    #for statistics_states in block_statistics_states:
                    #    plot.plot(timestamps[301:], statistics_states,
                    #              color=block_length_colors[block_length],
                    #              linestyle='--', alpha=0.1)

                    plot.fill_between(timestamps[301:], minima, maxima,
                                      color=block_length_colors[block_length],
                                      alpha=0.1)

                plot.xlabel('Time (ns)', fontsize=16)
                plot.ylabel('Number of contacts', fontsize=16)
                plot.legend(title='Block length',
                            fontsize=16, title_fontsize=16)
                plot.tight_layout()
                plot.savefig('contacts_' + str(bead_type_name_1) + '_' \
                                         + str(bead_type_name_2) + '.png')
                plot.clf()

    #print('Preparing report on simulations...')

    #dir_path = Path(args.dir)

    #if dir_path.is_dir():
    #    print('Report directory already exists.')
    #else:
    #    dir_path.mkdir(parents=False, exist_ok=False)
    #    print('Report directory created.')

    #report = dict()

    #seqs = {'SL' * 32, 'SSLL' * 16}

    #for seq in seqs:
        #report[seq] = dict()
        #report[seq]['contacts'] = count_contacts()
        #report[seq]['R_gyr'] = radius_of_gyration()
        #report[seq]['clustering'] = clustering()

    # To use as an importable module.
    #return report


def read_beads_from_LAMMPS_data(path: str) -> \
    list[tuple[int, float, float, float]]:
    """"Read molecule ID and coordinates of each bead in LAMMPS data file."""

    beads = []

    LAMMPS_data_file = open(path, 'r')
    LAMMPS_data_lines = LAMMPS_data_file.readlines()
    LAMMPS_data_file.close()

    read_atoms = False

    for line in LAMMPS_data_lines:
        if 'Atoms' in line:  # Atom section started.
            read_atoms = True

        if read_atoms:  # If the current section is atoms, read atomic params.
            tokens = line.strip().split()

            if len(tokens) == 10:
                molecule_id = int(tokens[1])
                bead_type = int(tokens[2])
                x, y, z = float(tokens[4]), float(tokens[5]), float(tokens[6])
                beads.append((molecule_id, bead_type, x, y, z))

        if 'Velocities' in line:  # Atom section ended.
            break

    return beads


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

    # Add specific commands to the CLI tool.
    subparsers = parser.add_subparsers()

    ################################################
    # Command to generate sequences for screening. #
    ################################################

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

    parser_generate_sequences.add_argument(
        "--blocky", type=int, default=None,
        help="""Generate all blocky sequences of the desired length. The block
                lengths will be all possible power of two, starting from one."""
    )

    parser_generate_sequences.set_defaults(func=generate_sequences)

    #########################################################################
    # Command to prepare the simulations for the sequences generated on the #
    # previous step.                                                        #
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
        "--repeats", type=int, default=10,
        help="Number of simulation repeats for each sequence. Default=10."
    )

    parser_prepare.add_argument(
        "--solution_time", type=float, default=0.01,
        help="Solution simulation time in microseconds. Default=0.01."
    )

    parser_prepare.add_argument(
        "--trap_time", type=float, default=0.01,
        help="Trap simulation time in microseconds. Default=0.01."
    )

    parser_prepare.add_argument(
        "--relaxation_time", type=float, default=0.01,
        help="Condensate simulation time in microseconds. Default=0.01."
    )

    parser_prepare.add_argument(
        "--damping_time", type=int, default=100,
        help="""Damping time of Langevin thermostat in picoseconds.
                Default=100."""
    )

    parser_prepare.add_argument(
        "--dt", type=float, default=10.0,
        help="Simulation time step in femtoseconds. Default=10."
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
        "--log_steps", type=int, default=10000,
        help="""Number of steps to make next record into the log and trajectory.
                Default=10000."""
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
        "--residue_radius_for_trap", type=float, default=10,
        help="""Single residue effective radius from which the trap condensate
                radius will be calculated [angstroms]."""
    )

    parser_prepare.set_defaults(func=prepare)

    ############################################################################
    # Command to run the simulations with settings defined on the previous two #
    # steps.                                                                   #
    ############################################################################

    parser_run = subparsers.add_parser(
        "run",
        help="Run the simulations."
    )

    parser_run.add_argument(
        "--dir", type=str,
        help="Directory path for the simulations."
    )

    parser_run.add_argument(
        "--repeats", type=int, default=10,
        help="Number of repeats per sequence."
    )

    parser_run.add_argument(
        "--lmp", type=str, default='lmp_mpi',
        choices=['lmp', 'lmp_kokkos',
                 'lmp_mpi', 'lmp_mpi_kokkos',
                 'lmp_serial', 'lmp_serial_kokkos'],
        help="Setup of LAMMPS run."
    )

    parser_run.set_defaults(func=run)

    ##################################################
    # Command to continue the performed simulations. #
    ##################################################

    parser_continue = subparsers.add_parser(
        "continue",
        help="Continue the simulations."
    )

    parser_continue.add_argument(
        "--dir", type=str,
        help="Directory path for the simulations."
    )

    parser_continue.add_argument(
        "--repeats", type=int, default=10,
        help="Number of repeats per sequence."
    )

    parser_continue.add_argument(
        "--time", type=float, default=0.1,
        help="Condensate simulation time in microseconds. Default=0.1."
    )

    parser_continue.add_argument(
        "--lmp", type=str, default='lmp_mpi',
        choices=['lmp', 'lmp_kokkos',
                 'lmp_mpi', 'lmp_mpi_kokkos',
                 'lmp_serial', 'lmp_serial_kokkos'],
        help="Setup of LAMMPS run."
    )

    parser_continue.add_argument(
        "--interrupted", type=bool, default=True,
        help="Whether to calculate only the interrupted trajectories."
    )

    parser_continue.set_defaults(func=cmd_continue)

    #################################################
    # Command to analyse the performed simulations. #
    #################################################

    parser_analyse = subparsers.add_parser(
        "analyse",
        help="Analyse the simulations."
    )

    parser_analyse.add_argument(
        "--dir", type=str,
        help="Directory path of the simulations."
    )

    parser_analyse.add_argument(
        "--cluster", action='store_true',
        help="Clustering."
    )

    parser_analyse.add_argument(
        "--gyration", action='store_true',
        help="Radius of gyration."
    )

    parser_analyse.add_argument(
        "--contacts", action='store_true',
        help="Contact statistics."
    )

    parser_analyse.add_argument(
        "--calculate", action='store_true',
        help="Calculate, not only plot."
    )

    parser_analyse.set_defaults(func=analyse)

    #####################
    # Run the CLI tool. #
    #####################

    # Parse the input and follow the commands.
    args = parser.parse_args()
    args.func(args)
