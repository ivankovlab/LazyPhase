import os
import continue_condensate


def continue_condensates_throughput(dir: str, lmp: str, repeats: int=10,
                                    interrupted: bool=True):
    """Run high-throughput condensate simulations for a set of sequences.
    dir : str
    lmp : str
    repeats : int
        How many times to repeat the simulation for the single sequence.
    interrupted : bool
        Whether to calculate only interrupted.
    """

    if dir[-1] != '/':
        dir += '/'

    try:
        file_seqs = open(dir + 'seqs.txt', 'r')
        seqs = file_seqs.readlines()
        file_seqs.close()
    except:
        print('Can not open ' + dir + 'seqs.txt' + '.')
        return

    seqs = [seq[:-1] for seq in seqs]

    num_seqs = len(seqs)

    print('Condensate simulations for', num_seqs, 'sequences will be run.')

    cnt_seq = 0

    for seq in seqs:
        for num in range(1, repeats+1):
            if os.path.isfile(dir + 'output/' + seq + '_' + str(num) + \
                              '_prolong.data'):
                print(seq, num, 'exists.')
            else:
                print(seq, num, 'does not exist.')
                print('Continue one simulation...')
                continue_condensate.continue_condensate(dir, lmp, seq, num)

        cnt_seq += 1
        print(cnt_seq, '/', num_seqs, 'sequences proccessed.')
