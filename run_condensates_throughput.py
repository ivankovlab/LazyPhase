from random import randint, shuffle
import run_condensate


def generate_random_string(N_S: int, S_L: int) -> str:
    """
    Generate a string with exactly N_S 'S' characters and S_L 'L' characters,
    arranged in a random order.
    """
    symbols = ['S'] * N_S + ['L'] * S_L
    shuffle(symbols)
    return ''.join(symbols)


def generate_seqs():
    file_seqs = open('seqs.txt', 'w')

    seqs = set([
        'SL' * 32,
        'SSLL' * 16,
        ('S' * 4 + 'L' * 4) * 8,
        ('S' * 8 + 'L' * 8) * 4,
        'S' * 16 + 'L' * 16 + 'S' * 16 + 'L' * 16,
        'S' * 32 + 'L' * 32,
        ('S' * 3 + 'L' * 1) * 16,
        'S' * 48 + 'L' * 16,
        'L' * 8 + 'S' * 48 + 'L' * 8,
        ('S' * 24 + 'L' * 8) * 2,
        ('L' * 6 + 'S' * 20 + 'L' * 2 + 'S' * 4) * 2,
        ('S' * 12 + 'L' * 4) * 4,
        'L' * 8 + 'S' * 24 + ('S' * 3 + 'L' * 1) * 8
    ])

    for i in range(10):
        N_S = randint(16, 32)
        seqs.add(generate_random_string(N_S, 64 - N_S))

    for seq in seqs:
        file_seqs.write(seq + '\n')

    file_seqs.close()


def run_condensates_throughput(repeats=3):
    """Run high-throughput condensate simulations for a set of sequences.
    repeats : int
        How many times to repeat the simulation for the single sequence.
    """

    file_seqs = open('seqs.txt', 'r')
    seqs = file_seqs.readlines()
    file_seqs.close()

    seqs = [seq[:-1] for seq in seqs]

    for seq in seqs:
        if len(seq) != 64:
            raise ValueError('Sequences must be 64 residues long.')

    alphabet = ['S', 'L']

    print('Condensate simulations for', len(seqs), 'sequences will be run.')

    for seq in seqs:
        for num in range(1, repeats+1):
            run_condensate.run_condensate(seq, num)


run_condensates_throughput()
