import numpy as np


def get_angles(bonds):
    """ Iterate through bonds to get angles.

    Bonds should contain no duplicates.
    """
    angles = []
    for i1, b1 in enumerate(bonds):
        for i2, b2 in enumerate(bonds):
            if i2 > i1:
                shared_atom = list(set(b1) & set(b2))
                if len(shared_atom) > 0:
                    atom1 = [b for b in b1 if b != shared_atom[0]][0]
                    atom2 = [b for b in b2 if b != shared_atom[0]][0]
                    other_atoms = sorted([atom1, atom2])
                    angles.append((other_atoms[0], shared_atom[0], other_atoms[1]))
    return sorted(angles)


def calculate_angle(p1, p2, p3):
    """ Calculate angle for three given points in space
      p2 ->  o
            / \
    p1 ->  o   o  <- p3
    """
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)

    v21 = p1 - p2
    v23 = p3 - p2
    angle = np.arccos(np.dot(v21, v23) / (np.linalg.norm(v21) * np.linalg.norm(v23)))
    return np.degrees(angle)
