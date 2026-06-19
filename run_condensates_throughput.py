import run_condensate


def run_condensates_throughput(dir: str, repeats=3):
    """Run high-throughput condensate simulations for a set of sequences.
    repeats : int
        How many times to repeat the simulation for the single sequence.
    """

    if dir[-1] != '/':
        dir += '/'

    try:
        file_seqs = open('seqs.txt', 'r')
        seqs = file_seqs.readlines()
        file_seqs.close()
    except:
        print('Can not open ' + dir + 'seqs.txt' + '.')
        return

    seqs = [seq[:-1] for seq in seqs]

    # if length_filter:
    #for seq in seqs:
    #    if len(seq) != 64:
    #        raise ValueError('Sequences must be 64 residues long.')

    #if alphabet_filter:
    #alphabet = ['S', 'L']

    num_seqs = len(seqs)

    print('Condensate simulations for', num_seqs, 'sequences will be run.')

    cnt_seq = 0

    for seq in seqs:
        for num in range(1, repeats+1):
            run_condensate.run_condensate(seq, num)

        cnt_seq += 1
        print(cnt_seq, '/', num_seqs, 'sequences proccessed.')
