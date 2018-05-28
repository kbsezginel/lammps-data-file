import math
import os
import copy
import numpy as np

from lammps_data import ase
from lammps_data.uff_table import UffTable
from lammps_data.angles import get_angles, calculate_angle


class MOF:
    """
    MOF class that holds coordinate, atom name, and unit cell information
    """
    def __init__(self, file_path, file_format='cif'):
        """
        Initialize MOF name, unit cell volume and parameters, atom names and coordinates, and
        unique atom names and coordinates.
        file_format -> ('cif' / 'dict' / 'mol2' / ...)
        """
        if file_format == 'dict':
            molecule = file_path
            self.atom_coors = molecule['atom_coors']
            self.atom_names = molecule['atom_names']
            self.name = molecule['name']
            self.path = file_path
        else:
            self.path = file_path
            self.name, file_format = os.path.splitext(os.path.basename(file_path))
            file_format = file_format[1:]
            self.ase_atoms, molecule = ase.read(file_path, input_format=file_format)

            self.atom_coors = molecule['atom_coors']
            self.atom_names = molecule['atom_names']
            self.uc_size = molecule['uc_size']
            self.uc_angle = molecule['uc_angle']
            self.atoms = molecule['atoms']
            self.separate_atoms()
            self.unit_cell_volume()

    def __repr__(self):
        return "<MOF object: %s>" % (self.name)

    def __str__(self):
        return self.name

    def __len__(self):
        return len(self.atom_coors)

    def __eq__(self, other):
        return (isinstance(other, self.__class__) and self.name == other.name)

    def __ne__(self, other):
        return not self.__eq__(other)

    def set_force_field(self, ff_type):
        """
        Initializes force field parameters according to unique atom names.
        """
        ff_param = get_ff_parameters(self.uniq_atom_names, ff_type)
        self.sigma = []
        self.epsilon = []
        for i in range(len(ff_param)):
            self.uniq_atom_names[i] = ff_param[i][0]
            self.sigma.append(ff_param[i][1])
            self.epsilon.append(ff_param[i][2])

    def calculate_vectors(self):
        """ Calculates unit cell vectors and edge points """
        self.uc_vectors = Packing.uc_vectors(self.uc_size, self.uc_angle)
        self.edge_points = Packing.edge_points(self.uc_vectors)

    def extend_unit_cell(self, cut_off=None, pack=None):
        """
        Extends unit cell of a MOF object according to a given cut_off value or packing list.
        If both are given pack argument has priority over cut_off.
        - pack=[x, y, z]
            Extends unit cell [(x) x (y) x (z)]
        - cut_off=*number
            The structure is extended so that all the points in the center unit cell are
            at least *cut_off Angstroms away.
        This enables energy map calculation by providing coordinates for surrounding atoms.
        Also enables checking for collisions in the extended unit cell of mobile layer
        and energy map.
        The *ext_cut_off parameter in the sim_par input file determines the amount of packing.
        - A high cut off value such as 100 Angstrom can be used to ensure there is no collisions
        between the interpenetrating layers.
        """
        self.calculate_vectors()
        if pack is not None:
            self.packing_factor = pack
        elif cut_off is not None:
            self.packing_factor = Packing.factor(self.uc_size, cut_off)
        else:
            print('Please enter packing factor or cut_off value!')
        trans_vec = Packing.translation_vectors(self.packing_factor, self.uc_vectors)
        self.packed_coors = Packing.uc_coors(trans_vec, self.packing_factor, self.uc_vectors, self.atom_coors)

        extended_structure = {'atom_names': [], 'atom_coors': [], 'atoms': [], 'name': self.name}
        for unit_cell in self.packed_coors:
            for coor_index, coor in enumerate(unit_cell):
                atom_name = self.atom_names[coor_index]
                atom_number = self.atoms[coor_index][3]
                extended_structure['atom_names'].append(atom_name)
                extended_structure['atom_coors'].append(coor)
                extended_structure['atoms'].append((coor[0], coor[1], coor[2], atom_number))

        packed_mof = copy.deepcopy(self)
        packed_mof.uc_size = np.array(self.uc_size) * np.array(self.packing_factor)
        packed_mof.uc_vectors = np.array(self.uc_vectors) * np.array(self.packing_factor)
        packed_mof.atoms = extended_structure['atoms']
        packed_mof.atom_names = extended_structure['atom_names']
        packed_mof.atom_coors = extended_structure['atom_coors']
        return packed_mof

    def separate_atoms(self):
        """
        Identifies different types of atoms in given atom coordinates and atom names
        Returns unique atom coordinates and unique atom names
        """
        uniq_atom_names = []
        for atom in self.atom_names:
            if atom not in uniq_atom_names:
                uniq_atom_names.append(atom)

        uniq_atom_coors = [[] for i in range(len(uniq_atom_names))]
        for atom_index in range(len(self.atom_coors)):
            for uniq_atom_index, uniq_atom in enumerate(uniq_atom_names):
                if self.atom_names[atom_index] == uniq_atom:
                    uniq_atom_coors[uniq_atom_index].append(self.atom_coors[atom_index])

        self.uniq_atom_names = uniq_atom_names
        self.uniq_atom_coors = uniq_atom_coors

    def calculate_cut_off(self):
        """
        Calculate cut-off radius as Rc = L/2 from a given MOF object.
        """
        width_a = self.ucv / (self.uc_size[1] * self.uc_size[2] / math.sin(math.radians(self.uc_angle[0])))
        width_b = self.ucv / (self.uc_size[0] * self.uc_size[2] / math.sin(math.radians(self.uc_angle[1])))
        width_c = self.ucv / (self.uc_size[0] * self.uc_size[1] / math.sin(math.radians(self.uc_angle[2])))
        self.cut_off = min(width_a / 2, width_b / 2, width_c / 2)

    def unit_cell_volume(self):
        """
        Calculates unit cell volume of a given MOF object.
        """
        a = self.uc_size[0]
        b = self.uc_size[1]
        c = self.uc_size[2]
        alp = math.radians(self.uc_angle[0])
        bet = math.radians(self.uc_angle[1])
        gam = math.radians(self.uc_angle[2])

        volume = 1 - math.cos(alp)**2 - math.cos(bet)**2 - math.cos(gam)**2
        volume += 2 * math.cos(alp) * math.cos(bet) * math.cos(gam)
        volume = a * b * c * math.sqrt(volume)
        frac_volume = volume / (a * b * c)

        self.ucv = volume
        self.frac_ucv = frac_volume

    def pbc_parameters(self):
        """
        Calculates constants used in periodic boundary conditions.
        """
        uc_cos = [math.cos(math.radians(a)) for a in self.uc_angle]
        uc_sin = [math.sin(math.radians(a)) for a in self.uc_angle]
        a, b, c = self.uc_size
        v = self.frac_ucv

        xf1 = 1 / a
        xf2 = - uc_cos[2] / (a * uc_sin[2])
        xf3 = (uc_cos[0] * uc_cos[2] - uc_cos[1]) / (a * v * uc_sin[2])
        yf1 = 1 / (b * uc_sin[2])
        yf2 = (uc_cos[1] * uc_cos[2] - uc_cos[0]) / (b * v * uc_sin[2])
        zf1 = uc_sin[2] / (c * v)
        self.to_frac = [xf1, xf2, xf3, yf1, yf2, zf1]

        xc1 = a
        xc2 = b * uc_cos[2]
        xc3 = c * uc_cos[1]
        yc1 = b * uc_sin[2]
        yc2 = c * (uc_cos[0] - uc_cos[1] * uc_cos[2]) / uc_sin[2]
        zc1 = c * v / uc_sin[2]
        self.to_car = [xc1, xc2, xc3, yc1, yc2, zc1]

    def clone(self):
        """
        Clones MOF object into a new MOF object
        """
        if type(self.path) is str:
            clone_mof = MOF(self.path)
        elif type(self.path) is dict:
            clone_mof = MOF(self.path, file_format='dict')
        if hasattr(self, 'packing_factor'):
            clone_mof.extend_unit_cell(pack=self.packing_factor)
        return clone_mof

    def join(self, new_mof, colorify=False, atom_color=['C', 'O']):
        """
        Combines atom names and coordinates of two given structure dictionaries.
        The structure dictionaries should be in following format:
         >>> structure = {'atom_names': [*atom names], 'atom_coors':[*atom coors]}
        """
        base_atom_name = atom_color[0]
        new_atom_name = atom_color[1]
        joined_atom_names = []
        joined_atom_coors = []
        joined_ase_atoms = self.ase_atoms.copy()

        if colorify:
            joined_ase_atoms.set_chemical_symbols([base_atom_name] * len(self.atom_names))
            for coor in new_mof.atom_coors:
                joined_atom_names.append(new_atom_name)
                joined_atom_coors.append(coor)
                joined_ase_atoms.append(ase.ase_atom(symbol=new_atom_name, position=coor))

            for coor in self.atom_coors:
                joined_atom_names.append(base_atom_name)
                joined_atom_coors.append(coor)

        else:
            for atom, coor in zip(new_mof.atom_names, new_mof.atom_coors):
                joined_atom_names.append(atom)
                joined_atom_coors.append(coor)
                joined_ase_atoms.append(ase.ase_atom(symbol=atom, position=coor))

            for atom, coor in zip(self.atom_names, self.atom_coors):
                joined_atom_names.append(atom)
                joined_atom_coors.append(coor)

        joined_structure = {'atom_names': joined_atom_names, 'atom_coors': joined_atom_coors}
        joined_structure['name'] = self.name + '_' + new_mof.name
        joined_mof = MOF(joined_structure, file_format='dict')
        joined_mof.ase_atoms = joined_ase_atoms

        return joined_mof

    def export(self, export_dir, file_format='xyz'):
        """
        Export MOF atom coordinates and names in .xyz format.

        Example usage:
         >>> mof.export(export_dir, file_format='xyz')
        """
        if file_format == 'xyz':
            xyz_path = os.path.join(export_dir, self.name + '.xyz')
            if os.path.exists(xyz_path):
                os.remove(xyz_path)

            with open(xyz_path, 'w') as xyz_file:
                xyz_file.write(str(len(self.atom_coors)) + '\n')
                xyz_file.write(self.name + '\n')
                for atom, coor in zip(self.atom_names, self.atom_coors):
                    xyz_file.write(atom + ' ' + str(coor[0]) + ' ' + str(coor[1]) + ' ' + str(coor[2]) + '\n')

        else:
            file_path = os.path.join(export_dir, self.name + '.' + file_format)
            ase.write(file_path, self.ase_atoms, file_format=file_format)

    def topology(self):
        """
        Get topology of the MOF object.
        --------------------------------------------------------------------------------------------
        Instead of determining bonds and agles separately it might be better idea to get them all
        together here since we need them before atom typing. Also we need to figure out which bonds
        are crossing the boundaries to calculate all angles correctly.
        If all topology is determined together it would also be more efficient since we can
        calculate the angles and bonds at the same time.
        --------------------------------------------------------------------------------------------
        # NOT COMPLETED !!!
        self.bonds, self.periodic_bonds = [], []
        self.angles, self.periodic_angles = [], []
        self.dihedrals, self.periodic_dihedrals = [], []
        self.calculate_vectors()

        atoms_list = list(enumerate(self.atoms))
        max_rad = max([rcov[a[1][3]] for a in atoms_list])
        translations = [[0, 0, 0]]
        for v in self.cell_vectors:
            translations.append(v)
            translations.append([-v[0], -v[1], -v[2]])

        for i in range(0, len(atoms_list)):
            max_cutoff = rcov[atoms_list[i][1][3]] + max_rad + RADIUS_BUFFER
            atom1 = atoms_list[i][1][:3]

            neighbors = dict(bonds=[], atoms=[], bond_lengths=[])
            for j in range(i + 1, len(atoms_list)):
                atom2 = atoms_list[j][1][:3]
                max_bond_distance = (rcov[atoms_list[i][1][3]] + rcov[atoms_list[j][1][3]] + RADIUS_BUFFER)
                for t in translations:
                    atom2_t = [atom2[0] + t[0], atom2[1] + t[1], atom2[2] + t[2]]
                    distance = ((atom1[0] - atom2_t[0]) ** 2 +
                                (atom1[1] - atom2_t[1]) ** 2 +
                                (atom1[2] - atom2_t[2]) ** 2) ** 0.5
                    if distance >= 0.16 and distance <= max_bond_distance:
                        bonded_atoms = sorted((atoms_list[i][0], atoms_list[j][0]))

                        neighbors['atoms'].append((atom2_t[0], atom2_t[1], atom2_t[2], j))
                        neighbors['bonds'].append((bonded_atoms[0], bonded_atoms[1]))
                        neighbors['bond_lengths'].append(distance)
                        self.bonds.append((bonded_atoms[0], bonded_atoms[1], distance))
                        if t = [0, 0, 0]:
                            self.periodic_bonds.append((bonded_atoms[0], bonded_atoms[1], distance))
                        break
            # Get angles
            neighbors['angles'] = get_angles(neighbors['bonds'])
            neighbors['angle_values'] = []
            for angle in neighbors['angles']:
                p1 = [a[3] for a in neighbors['atoms'] if a[3] == angle[0]]  # There must a better way of doing this
                p2 = ...
                p3 = ...
                neighbors['angle_values'].append(calculate_angle(p1, p2, p3))
        """
        return None

    def assign_UFF_atom_types(self, tolerance=10):
        uff = UffTable()
        self.atom_types = []
        self.unique_atom_types = dict()
        for atom_idx, atom in enumerate(self.atoms):
            atom_bonds = [x for x in self.bonds if atom_idx in x]
            n_bonds = len(atom_bonds)
            element = self.atom_names[atom_idx]
            """

            Figure out periodic bonds and calculate angles for those separately
            -> Alternatively first go through each atom and calculate bonds + angles

            """
            atom_angles = [calculate_angle(self.atoms[a[0]][:3], self.atoms[a[1]][:3], self.atoms[a[2]][:3]) for a in get_angles(atom_bonds)]
            avg_angle = sum(atom_angles) / len(atom_angles)
            if all([math.isclose(angle, avg_angle, abs_tol=tolerance) for angle in atom_angles]):
                uff_row = uff.lookup(element, n_bonds, avg_angle)
                self.atom_types.append(uff_row['type'])
                if uff_row['type'] not in self.unique_atom_types:
                    self.unique_atom_types[uff_row['type']] = uff_row
            elif n_bonds == 4:
                # Check square planar case
                ortho = all([math.isclose(angle, avg_angle, abs_tol=tolerance) for angle in sorted(atom_angles)[:4]])
                linear = all([math.isclose(angle, avg_angle, abs_tol=tolerance) for angle in sorted(atom_angles)[4:]])
                if ortho and linear:
                    print('Found square planar')
                    print('Element: %s | NumBonds: %i | AvgAngle: % .3f' % (self.atom_names[atom_idx], len(atom_bonds), avg_angle))
                    print('All angles:', atom_angles)
                    uff_row = uff.lookup(element, n_bonds, avg_angle)
                    self.atom_types.append(uff_row['type'])
                    if uff_row['type'] not in self.unique_atom_types:
                        self.unique_atom_types[uff_row['type']] = uff_row
                else:
                    print('Element: %s | NumBonds: %i | AvgAngle: % .3f' % (self.atom_names[atom_idx], len(atom_bonds), avg_angle))
                    print('All angles:', atom_angles)
                    raise AtomTypingError("%s atom does not fit square planar geometry" % (self.atom_names[atom_idx]))
            else:
                # Throw an error
                print('Element: %s | NumBonds: %i | AvgAngle: % .3f' % (self.atom_names[atom_idx], len(atom_bonds), avg_angle))
                print('All angles:', atom_angles)
                raise AtomTypingError("%s atom cannot be typed" % (self.atom_names[atom_idx]))


class Packing:
    """
    Packing class containing functions used for packing unit cells
    """
    @classmethod
    def factor(cls, uc_size, cut_off):
        """
        Calculate packing factor for given unit cell size and cut off radius
        """
        packing_factor = []
        for uc_length in uc_size:
            packing_factor.append(math.ceil(cut_off / uc_length) * 2 + 1)
        return packing_factor

    @classmethod
    def uc_vectors(cls, uc_size, uc_angle):
        """
        Calculate unit cell vectors for given unit cell size and angles
        """
        a = uc_size[0]
        b = uc_size[1]
        c = uc_size[2]
        alpha = math.radians(uc_angle[0])
        beta = math.radians(uc_angle[1])
        gamma = math.radians(uc_angle[2])

        x_v = [a, 0, 0]
        y_v = [b * math.cos(gamma), b * math.sin(gamma), 0]
        z_v = [0.0] * 3
        z_v[0] = c * math.cos(beta)
        z_v[1] = (c * b * math.cos(alpha) - y_v[0] * z_v[0]) / y_v[1]
        z_v[2] = math.sqrt(c * c - z_v[0] * z_v[0] - z_v[1] * z_v[1])
        uc_vectors = [x_v, y_v, z_v]
        return uc_vectors

    @classmethod
    def translation_vectors(cls, packing_factor, uc_vectors):
        """
        Calculate translation vectors for given packing factor and uc vectors
        """
        packing_amount = []
        for x in range(packing_factor[0]):
            for y in range(packing_factor[1]):
                for z in range(packing_factor[2]):
                    packing_amount.append([x, y, z])

        x_v = uc_vectors[0]
        y_v = uc_vectors[1]
        z_v = uc_vectors[2]
        translation_vectors = []
        for pack in packing_amount:
            x_trans = x_v[0] * pack[0] + y_v[0] * pack[1] + z_v[0] * pack[2]
            y_trans = x_v[1] * pack[0] + y_v[1] * pack[1] + z_v[1] * pack[2]
            z_trans = x_v[2] * pack[0] + y_v[2] * pack[1] + z_v[2] * pack[2]
            translation_vectors.append([x_trans, y_trans, z_trans])

        return translation_vectors

    @classmethod
    def uc_coors(cls, translation_vectors, packing_factors, uc_vectors, atom_coors):
        """
        Calculate packed coordinates for given:
        - translation vectors  - packing factor     - unit cell vectors    - atom coordinates
        """
        x_v = uc_vectors[0]
        y_v = uc_vectors[1]
        z_v = uc_vectors[2]

        packed_coors = [[] for i in range(len(translation_vectors))]
        translation_factor = []
        origin_trans_vec = []
        for dim, factor in enumerate(packing_factors):
            translation_factor.append((factor - 1) / 2)
            origin_translation = (factor - 1) / 2 * x_v[dim] + (factor - 1) / 2 * y_v[dim]
            origin_translation += (factor - 1) / 2 * z_v[dim]
            origin_trans_vec.append(origin_translation)

        packing_index = 0
        for translation in translation_vectors:
            for coor in atom_coors:
                x = coor[0] + translation[0] - origin_trans_vec[0]
                y = coor[1] + translation[1] - origin_trans_vec[1]
                z = coor[2] + translation[2] - origin_trans_vec[2]
                packed_coors[packing_index].append([x, y, z])
            packing_index += 1

        return packed_coors

    @classmethod
    def edge_points(cls, uc_vectors):
        """
        Calculate coordinates of unit cell edges in order of:
        (0, 0, 0) - (a, 0, 0) - (b, 0, 0) - (c, 0, 0)
        (a, b, 0) - (0, b, c) - (a, 0, c) - (a, b, c)
        """
        uc_edges = []
        uc_edges.append([0, 0, 0])

        # (a, 0, 0) - (b, 0, 0) - (c, 0, 0)
        for vec in uc_vectors:
            uc_edges.append(vec)

        # (a, b, 0)
        uc_edges.append([uc_vectors[0][0] + uc_vectors[1][0], uc_vectors[0][1] + uc_vectors[1][1],
                        uc_vectors[0][2] + uc_vectors[1][2]])
        # (0, b, c)
        uc_edges.append([uc_vectors[1][0] + uc_vectors[2][0], uc_vectors[1][1] + uc_vectors[2][1],
                        uc_vectors[1][2] + uc_vectors[2][2]])
        # (a, 0, c)
        uc_edges.append([uc_vectors[0][0] + uc_vectors[2][0], uc_vectors[0][1] + uc_vectors[2][1],
                        uc_vectors[0][2] + uc_vectors[2][2]])
        # (a, b, c)
        uc_edges.append([uc_vectors[0][0] + uc_vectors[1][0] + uc_vectors[2][0],
                        uc_vectors[0][1] + uc_vectors[1][1] + uc_vectors[2][1],
                        uc_vectors[0][2] + uc_vectors[1][2] + uc_vectors[2][2]])
        return uc_edges


class AtomTypingError(Exception):
    """ Raised when atom type could not be determined """
    pass
