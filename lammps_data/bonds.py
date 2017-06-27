import math

# From Pyykko doi: 10.1002/chem.200800987
rcov = [
    0.18, 0.32, 0.46,
    1.33, 1.02, 0.85, 0.75, 0.71, 0.63, 0.64, 0.67,
    1.55, 1.39, 1.26, 1.16, 1.11, 1.03, 0.99, 0.96,
    1.96, 1.71, 1.48,
    1.36, 1.34, 1.22, 1.19, 1.16, 1.11, 1.10, 1.12, 1.18,
    1.24, 1.21, 1.21, 1.16, 1.14, 1.17,
    2.10, 1.85,
    1.63, 1.54, 1.47, 1.38, 1.28, 1.25, 1.25, 1.20, 1.28, 1.36,
    1.42, 1.40, 1.40, 1.36, 1.33, 1.31,
    2.32, 1.96,
    1.80, 1.63, 1.76, 1.74, 1.73, 1.72, 1.68,
    1.69, 1.68, 1.67, 1.66, 1.65, 1.64, 1.70,
    1.62, 1.52, 1.46, 1.37, 1.31, 1.29, 1.22, 1.23, 1.24, 1.33,
    1.44, 1.44, 1.51, 1.45, 1.47, 1.42,
    2.23, 2.01,
    1.86, 1.75, 1.69, 1.70, 1.71, 1.72, 1.66,
    1.66, 1.68, 1.68, 1.65, 1.67, 1.73, 1.76,
    1.61, 1.57, 1.49, 1.43, 1.41, 1.34, 1.29, 1.28, 1.21, 1.22,
    1.36, 1.43, 1.62, 1.75, 1.65, 1.57
]


# RADIUS_BUFFER = 0.45


def extrapolate_bonds(atoms, RADIUS_BUFFER=0.45):
    bonds = []

    atoms_z = list(enumerate(atoms))
    atoms_z = sorted(atoms_z, key=lambda x: x[1][2])

    max_rad = max([rcov[a[1][3]] for a in atoms_z])

    for i in range(0, len(atoms_z)):
        max_cutoff = rcov[atoms_z[i][1][3]] + max_rad + RADIUS_BUFFER
        for j in range(i + 1, len(atoms_z)):
            if abs(atoms_z[j][1][2] - atoms_z[i][1][2]) > max_cutoff:
                break

            distance = ((atoms_z[j][1][0] - atoms_z[i][1][0])**2 +
                        (atoms_z[j][1][1] - atoms_z[i][1][1])**2 +
                        (atoms_z[j][1][2] - atoms_z[i][1][2])**2) ** 0.5

            max_bond_distance = (rcov[atoms_z[i][1][3]] + rcov[atoms_z[j][1][3]] + RADIUS_BUFFER)
            if distance >= 0.16 and distance <= max_bond_distance:
                bonds.append(tuple(sorted((atoms_z[i][0], atoms_z[j][0]))))

    return bonds


def extrapolate_periodic_bonds(atoms, cell_vectors, RADIUS_BUFFER=0.45):
    bonds = []

    atoms_z = list(enumerate(atoms))

    max_rad = max([rcov[a[1][3]] for a in atoms_z])

    translations = [[0, 0, 0]]
    for v in cell_vectors:
        translations.append(v)
        translations.append([-v[0], -v[1], -v[2]])

    for i in range(0, len(atoms_z)):
        max_cutoff = rcov[atoms_z[i][1][3]] + max_rad + RADIUS_BUFFER
        atom1 = atoms_z[i][1][:3]
        for j in range(i + 1, len(atoms_z)):
            atom2 = atoms_z[j][1][:3]
            max_bond_distance = (rcov[atoms_z[i][1][3]] + rcov[atoms_z[j][1][3]] + RADIUS_BUFFER)
            for t in translations:
                atom2_t = [atom2[0] + t[0], atom2[1] + t[1], atom2[2] + t[2]]
                distance = ((atom1[0] - atom2_t[0]) ** 2 +
                            (atom1[1] - atom2_t[1]) ** 2 +
                            (atom1[2] - atom2_t[2]) ** 2) ** 0.5
                if distance >= 0.16 and distance <= max_bond_distance:
                    bonds.append(tuple(sorted((atoms_z[i][0], atoms_z[j][0]))))
                    break

    return bonds


def calculate_distance(p1, p2):
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    distance = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)


def calculate_bond_length(bond, atoms, limit=(0.5, 5)):
    """ Calculate bond distance for given bond and atom list """
    a1, a2 = bond
    x1, y1, z1 = atoms[a1][:3]
    x2, y2, z2 = atoms[a2][:3]
    distance = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
    if distance > limit[1]:
        raise BondLengthLimitError("Bond length of %s between %i and %i is higher than limit %s"
                                   % (str(distance), bond[0], bond[1], str(limit[1])))
    elif distance < limit[0]:
        raise BondLengthLimitError("Bond length of %s between %i and %i is lower than limit %s"
                                   % (str(distance), bond[0], bond[1], str(limit[0])))
    else:
        return distance


class BondLengthLimitError(Exception):
    """ Raised when bond length is higher or lower than allowed limit """
    pass
