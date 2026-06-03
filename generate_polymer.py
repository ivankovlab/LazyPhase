'''
Generate polypeptide-like chain in Moltemplate format.
'''

import numpy as np

alphabet = ['S', 'L']

seq_1_1 = 'SL' * 8
seq_4_4 = ('S' * 4 + 'L' * 4) * 2
seq_3_1_1_3 = 'SSSL' * 2 + 'SLLL' * 2
seq_8_8 = 'S' * 8 + 'L' * 8

seq = seq_8_8

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


print('1 1', len(seq_1_1))
print('4 4', len(seq_4_4))
print('3 1 1 3', len(seq_3_1_1_3))
print('8 8', len(seq_8_8))

print('w = 4')

print('1 1',     get_seq_het(seq_1_1,     4))
print('4 4',     get_seq_het(seq_4_4,     4))
print('3 1 1 3', get_seq_het(seq_3_1_1_3, 4))
print('8 8',     get_seq_het(seq_8_8,     4))

print('w = 8')

print('1 1',     get_seq_het(seq_1_1,     8))
print('4 4',     get_seq_het(seq_4_4,     8))
print('3 1 1 3', get_seq_het(seq_3_1_1_3, 8))
print('8 8',     get_seq_het(seq_8_8,     8))

# Build the polymer from monomers.
output = '  mon1 = new ' + seq[0] + '\n'
for n in range(1, len(seq)):
    if n % 2:
        rot = '180.0'
    else:
        rot = '0.0'

    mov = str(round(3.2 * n, 1))

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
