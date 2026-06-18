'''
HOW TO PREPARE INPUT FOR pdb_fine.py?

(1) Open VMD.
In VMD terminal:
    (2) > topo readlammpsdata path_data/result.data

    To remove PBC effect and join all clusters together:
    (1) Open VMD.
    In VMD terminal:
        (2) > pbc join connected -all
        (3) > pbc wrap -centersel all -center com -compound residue -all
        Somteimes needed:
        pbc unwrap
        pbc wrap
        pbc set {a b c}
    Repeat until all clusters will be connected.

    (3) > set sel [atomselect top all]
    (4) > $sel writepdb path_pdb/result.pdb
(5) Run this script "python pdb_fine.py" or "python3 pdb_fine.py" with desired
    path_data/result.data and path_pdb/result.pdb, used during the steps above.
'''

import argparse

cli = argparse.ArgumentParser(
                    prog='PDBfine',
                    description='Add bonds to PDBs prepared from LAMMPS \
                                 data using VMD.',
                    epilog='(C) Egor Vasilenko, 2026')

cli.add_argument('-p', '--pdb')
cli.add_argument('-d', '--data')

args = cli.parse_args()

pdb = args.pdb
data = args.data

# Read PDB data.
f = open(pdb, 'r')
lines = f.readlines()
f.close()

# Read LAMMPS data.
f_LAMMPS = open(data, 'r')
lines_LAMMPS = f_LAMMPS.readlines()
f_LAMMPS.close()

bonds = []  # Bonds in LAMMPS data as pairs (i, j).

bond_flag = False  # Whether the state is bond reading.

for line in lines_LAMMPS:  # Read the whole LAMMPS data file.
    tokens = line.split()

    if bond_flag and len(tokens) > 1:
        bonds.append((tokens[2], tokens[3]))  # Add a new bond to the list.

    if len(tokens) > 0 and tokens[0] == 'Bonds':  # Start reading bonds.
        bond_flag = True
    elif len(tokens) > 0 and tokens[0] == 'Angles':  # 'Angles' section follows
                                                     # 'Bonds' section.
        bond_flag = False

# Define all atoms as heteroatoms.
f = open(pdb, 'w')
for line in lines[:-1]:
    tokens = line.split()
    if line[0] == 'A':
        f.write('HETATM' + line[6:])  # Change 'ATOM' to 'HETATM' because we
                                      # define the beads arbitrary.
    else:
        f.write(line)  # Conserve all other settings.

# Append the bond lines to the PDB structure file.
for bond in bonds:
    f.write('CONECT ' + (4 - len(str(bond[0]))) * ' ' + str(bond[0]) + ' ' + \
            (4 - len(str(bond[1]))) * ' ' + str(bond[1]) + '\n')

f.write(lines[-1])
f.close()
