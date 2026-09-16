'''
Generate polypeptide-like chain in Moltemplate format and prepare the simulation
cell for the condensate studies.
'''


import numpy as np
from random import randint, uniform
import os
from LazyPhase import generate_substituted_file
import subprocess


def continue_condensate(dir: str, lmp: str, seq: str, num: int):
    """Run condensate simulation for one sequence. Each run assumes creation
    of a new running protocol specific for the sequence and containing random
    seed. The polymers are placed randomly. One simulation will be run, the
    names of all output files before extensions will be in format
    SEQUENCE_NUMBER, like SEQUENCE_NUMBER.data, SEQUENCE_NUMBER.lammpstrj,
    SEQUENCE_NUMBER.log.

    seq : str
        Sequence of the polymer written as one-letter monomer names.
    lmp : str
        Specifications on LAMMPS launch.
    num : int
        Number of repeats for the sequence. Each repeat will be run with the
        same sequence, but newly generated seed.
    """

    ###########################################################################
    # Generate input for the sequence to write all output will sequence name. #
    ###########################################################################

    # Maximal seed value is taken as the maximal 32-bit integer.
    # The simulations of one serie differ only in sequence and seed.
    generate_substituted_file(dir + 'continue.run.in.nvt', dir + 'run.in.nvt', [
        ('SEQUENCE', seq + '_' + str(num)),
        ('SEED', str(randint(1, 65535)))
    ])

    ############################
    # Continue the simulation. #
    ############################

    print('Continue the simulation...')
    subprocess.run(['bash', 'run.sh', dir, lmp])
