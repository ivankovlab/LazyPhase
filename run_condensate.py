'''
Generate polypeptide-like chain in Moltemplate format and prepare the simulation
cell for the condensate studies.
'''

import numpy as np
from random import uniform
import os
from LazyPhase import generate_substituted_file


def run_condensate(seq: str, num: int, c=3):
    """Run condensate simulation for one sequence. Each run assumes creation
    of a new running protocol specific for the sequence and containing random
    seed. The polymers are placed randomly. One simulation will be run, the
    names of all output files before extensions will be in format
    SEQUENCE_NUMBER, like SEQUENCE_NUMBER.data, SEQUENCE_NUMBER.lammpstrj,
    SEQUENCE_NUMBER.log.

    seq : str
        Sequence of the polymer written as one-letter monomer names.
    num : int
        Number of repeats for the sequence. Each repeat will be run with the
        same sequence, but newly generated seed.
    c : int
        Number of polymers in line.
    """

    ###########################################################################
    # Generate input for the sequence to write all output will sequence name. #
    ###########################################################################

    # Maximal seed value is taken as the maximal 32-bit integer.
    # The simulations of one serie differ only in sequence and seed.
    generate_substituted_file('template.run.in.nvt', 'run.in.nvt',
                              [('SEQUENCE', seq + '_' + str(num)),
                               ('SEED', random.uniform(1, 2147483647))])

    ###############################################
    # Generate a template for the single polymer. #
    ###############################################

    # Build the polymer from monomers.
    output = '  mon1 = new ' + seq[0] + '\n'  # Put the first monomer in origin.
    for n in range(1, len(seq)):
        # Snake-like geometry.
        if n % 2:
            rot = '180.0'
        else:
            rot = '0.0'

        # Put the next two-bead monomer by cosine to make initial angle BBB
        # approximately 150 degrees.
        mov = str(round(3.2 * n, 1))

        # Write the monomer type and position.
        output += '  mon' + str(n+1) + ' = new ' + str(seq[n]) + \
                  '.rot(' + rot + ', 1,0,0).move(' + mov + ',0,0)' + '\n'

    output += '\n'

    # Write backbone bonds connecting C-alpha beads. They are not included in
    # the monomer descriptions, which include only bonds between backbonde and
    # side chain beads.
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

    # Length of a delineated polymer along the main backbone axis.
    polymer_length = (len(seq) - 1) * 3.2

    f = open('moltemplate_files/system.lt', 'w')

    f.write('import "polymer.lt"\n\n')
    f.write('write_once("Data Boundary") {\n')
    # The subcells for the polymers are sized with paddings of roughly two
    # monomer lengths.
    f.write('0 ' + str(round((polymer_length + 6.4) * c, 1)) + ' xlo xhi\n')
    f.write('0 ' + str(round((polymer_length + 6.4) * c, 1)) + ' ylo yhi\n')
    f.write('0 ' + str(round((polymer_length + 6.4) * c, 1)) + ' zlo zhi\n')
    f.write('}\n\n')

    polymer_num = 0  # Current polymer ID.

    # Iterate over all lattice units.
    for i in range(0, c):
        for j in range(0, c):
            for k in range(0, c):
                # Polymer ID to distinguish.
                polymer_num += 1

                # Coordinates of the subcell center, the center of polymer
                # will be placed there.
                x = (i + 0.5) * (polymer_length + 6.4)
                y = (j + 0.5) * (polymer_length + 6.4)
                z = (k + 0.5) * (polymer_length + 6.4)

                # Orient a polymer randomly using spherical coordinates.
                random_angle_1 = uniform(0, 2 * np.pi)
                random_angle_2 = uniform(0, 2 * np.pi)
                random_angle_3 = uniform(0, 360)
                random_x = np.cos(random_angle_1) * np.cos(random_angle_2)
                random_y = np.cos(random_angle_1) * np.sin(random_angle_2)
                random_z = np.sin(random_angle_1)

                # Initialize a new polymer.
                init_string = 'polymer_' + str(polymer_num) + ' = new Polymer'

                # Rotate the polymer randomly around its center.
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

    #########################################################
    # Run the final system assembly and run the simulation. #
    #########################################################

    os.system('./README_setup.sh && ./README_run.sh')
