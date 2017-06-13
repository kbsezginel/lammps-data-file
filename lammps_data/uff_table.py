import csv
import math


class UffTable:
    def __init__(self, uff_file='uff-parameters.csv'):
        self.load(uff_file)

    def load(self, uff_file='uff-parameters.csv'):
        """
        Load UFF parameters file.
        """
        uff = []
        with open(uff_file, 'r') as csvfile:
            csvreader = csv.reader(csvfile, delimiter=',')
            for row in csvreader:
                uff.append(row)

        self.parameters = dict()
        for atom_data in uff[2:]:
            row = dict(type=atom_data[0],
                       element=atom_data[1],
                       coordination=int(atom_data[2]),
                       bond=float(atom_data[3]),
                       angle=float(atom_data[4]),
                       distance=float(atom_data[5]),
                       energy=float(atom_data[6]),
                       scale=float(atom_data[7]),
                       charge=float(atom_data[8]))

            if row['element'] not in self.parameters:
                self.parameters[row['element']] = {}
            if row['coordination'] not in self.parameters[row['element']]:
                self.parameters[row['element']][row['coordination']] = {}
            if not self.parameters[row['element']][row['coordination']]:
                self.parameters[row['element']][row['coordination']] = []

            self.parameters[row['element']][row['coordination']].append(row)

    def lookup(self, element, n_bonds, angle, tolerance=10):
        """
        Look up element from UFF table by element name, number of bonds and angle.
        For angle matching absolute tolerance can be adjusted using tolerance keyword argument.
        """
        if element not in self.parameters:
            raise NonExistentElementError('Element not found in table')
        elif n_bonds not in self.parameters[element]:
            raise NonExistentCoordinationError('Coordination of %s not found for element %s' % (str(n_bonds), element))
        else:
            rows = []
            for r in self.parameters[element][n_bonds]:
                if math.isclose(r['angle'], angle, abs_tol=tolerance):
                    rows.append(r)
            if len(rows) == 0:
                raise NonExistentAngleError('Angle of %s not found for coordination %s of element %s' % (str(angle), str(n_bonds), element))
            elif len(rows) > 1:
                raise LookupMultiplicityError('Multiple (%i) angles found for tolerance of: %s' % (len(rows), str(tolerance)))
            else:
                return rows[0]


class NonExistentElementError(Exception):
    """ Raised when input element cannot be found in UFF parameters table """
    pass


class NonExistentCoordinationError(Exception):
    """ Raised when input coordination for input element cannot be found in UFF table """
    pass


class NonExistentAngleError(Exception):
    """ Raised when input angle for element and coordination is not within the tolerance value """
    pass


class LookupMultiplicityError(Exception):
    """ Raised when multiple rows are to be returned """
    pass
