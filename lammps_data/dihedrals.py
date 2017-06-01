

def get_dihedrals(bonds):
    """ Iterate through bonds to get dihedrals.

    Bonds should contain no duplicates.
    """
    dihedrals = []
    for i1, middle_bond in enumerate(bonds):
        atom1, atom2 = middle_bond
        atom1_bonds, atom2_bonds = [], []
        for i2, other_bond in enumerate(bonds):
            if atom1 in other_bond:
                atom1_bonds.append(other_bond)
            elif atom2 in other_bond:
                atom2_bonds.append(other_bond)

        for bond1 in atom1_bonds:
            for bond2 in atom2_bonds:
                atom0 = [b for b in bond1 if b != atom1][0]
                atom3 = [b for b in bond2 if b != atom2][0]
                dihedral = (min(atom0,atom3), atom1, atom2, max(atom0,atom3))
                # Make sure dihedral is not a loop
                if len(set(dihedral)) == 4:
                    dihedrals.append(tuple(dihedral))

    return sorted(dihedrals)
