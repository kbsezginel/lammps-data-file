

class UffTable:
    def __init__(self):
        return self

    def lookup(element, n_bonds, angle, tolerance=10):
        """
        Look up element from UFF table by element name, number of bonds and angle.
        For angle matching absolute tolerance can be adjusted using tolerance keyword argument.
        """
        element_rows = self.table[element]


class NonExistentElementError(Exception):
    pass


class NonExistentCoordinationError(Exception):
    pass


class NonExistentAngleError(Exception):
    pass


class LookupMultiplicityError(Exception):
    pass
"""

1. Center atom
2. Coordination number
3. Lookup all angle within tolerance of angle
- O rows: error
- 1 row: return rows
- >1 rows: error

Remove rows with same coordination numbers and angles
Add rows for square planar geometry for 180 degree angles

Warning for same atom being assigned multiple atom types

3 ways to assign:
- All UFF parameters
- All distance and angle from cif
- Distance and angle from UFF to calculate force constants, equilibrium cif values for other terms

Recognize symmetrical copies and making sure we use same parameters for those
"""
