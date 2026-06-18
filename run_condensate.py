'''
Generate polypeptide-like chain in Moltemplate format and prepare the simulation
cell for the condensate studies.
'''

import numpy as np
from random import uniform
import os


def get_seq_het(seq: str, w: int) -> dict:
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


def generate_substituted_file(path_in, path_out, pattern, substit):
    file_in = open(path_in, 'r')
    file_out = open(path_out, 'w')
    for line in file_in:
        file_out.write(line.replace(pattern, substit))
    file_in.close()
    file_out.close()


def run_condensate(seq: str, num: int):
    """Run condensate simulation for one sequence.
    seq : str
        Sequence of the polymer written as one-letter monomer names.
    """

    # Generate input for the sequence to write all output will sequence name.
    generate_substituted_file('template.run.in.nvt', 'run.in.nvt', 'NAME',
                              seq + '_' + str(num))

    ###############################################
    # Generate a template for the single polymer. #
    ###############################################

    # Build the polymer from monomers.
    output = '  mon1 = new ' + seq[0] + '\n'  # Put the first monomer in origin.
    for n in range(1, len(seq)):
        if n % 2:
            rot = '180.0'
        else:
            rot = '0.0'

        mov = str(round(3.2 * n, 1))  # Put the next two-bead monomer by cosine.

        output += '  mon' + str(n+1) + ' = new ' + str(seq[n]) + \
                  '.rot(' + rot + ', 1,0,0).move(' + mov + ',0,0)' + '\n'

    output += '\n'

    # Write backbone bonds connecting C-alpha beads.
    output += '  write("Data Bond List") {\n'

    for n in range(1, len(seq)):
        output += '    $bond:backbone' + str(n) + '  $atom:mon' + str(n) + \
                  '/ca   $atom:mon' + str(n+1) + '/ca\n'

    output += '  }\n'

    polymer_file = open('moltemplate_files/polymer.lt', 'w')
    polymer_file.write('import "monomer_S.lt"\n')
    polymer_file.write('import "monomer_L.lt"\n')
    polymer_file.write('\n')
    polymer_file.write('Polymer {\n')
    polymer_file.write('\n')
    polymer_file.write('  create_var {$mol}\n')
    polymer_file.write('\n')
    polymer_file.write(output)
    polymer_file.write('\n}')
    polymer_file.close()

    ###########################################################################
    # Create the simulating cell and fill it with randomly oriented polymers. #
    ###########################################################################

    polymer_length = (len(seq) - 1) * 3.2

    f = open('moltemplate_files/system.lt', 'w')

    # Number of polymers in line.
    c = 3

    f.write('import "polymer.lt"\n\n')
    f.write('write_once("Data Boundary") {\n')
    f.write('0 ' + str((polymer_length + 6.4) * c) + ' xlo xhi\n')
    f.write('0 ' + str((polymer_length + 6.4) * c) + ' ylo yhi\n')
    f.write('0 ' + str((polymer_length + 6.4) * c) + ' zlo zhi\n')
    f.write('}\n\n')

    polymer_num = 0
    for i in range(0, c):
        for j in range(0, c):
            for k in range(0, c):
                # Polymer ID to distinguish.
                polymer_num += 1

                # Coordinates of the subcell edge.
                x = (i + 0.5) * (polymer_length + 6.4)
                y = (j + 0.5) * (polymer_length + 6.4)
                z = (k + 0.5) * (polymer_length + 6.4)

                random_angle_1 = uniform(0, 2 * np.pi)
                random_angle_2 = uniform(0, 2 * np.pi)
                random_angle_3 = uniform(0, 360)
                random_x = np.cos(random_angle_1) * np.cos(random_angle_2)
                random_y = np.cos(random_angle_1) * np.sin(random_angle_2)
                random_z = np.sin(random_angle_1)

                # Initialize a new polymer.
                init_string = 'polymer_' + str(polymer_num) + ' = new Polymer'

                # Rotate randomly around its center.
                rot_string = '.move(' + \
                    str(-polymer_length / 2) + ',' + str(0) + ',' + str(0) + \
                             ').rot(' + str(random_angle_3) + ', ' + \
                    str(random_x) + ', ' + \
                    str(random_y) + ', ' + \
                    str(random_z) + ')'

                # Place in the dedicated subcell.
                placement_string = '.move(' + str(x) + ', ' + str(y) + ', ' + \
                                              str(z) + ')'

                # Write all geometrical transformations.
                f.write(init_string + rot_string + placement_string + '\n')

    f.close()

    os.system('./README_setup.sh && ./README_run.sh')

run_condensate(('S' * 1 + 'L' * 1) * 32, 1)
