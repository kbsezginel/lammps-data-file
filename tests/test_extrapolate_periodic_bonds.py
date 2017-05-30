from lammps_data.bonds import extrapolate_periodic_bonds
from lammps_data.bonds import extrapolate_bonds


def test_zero_cell_vectors_should_give_same_bonds_with_nonperiodic():
    cell_vectors = [[0] * 3] * 3
    atoms = [(0, 0, 0, 1)]
    assert extrapolate_bonds(atoms) == extrapolate_periodic_bonds(atoms, cell_vectors)


def test_periodic_atoms_too_close_should_not_be_bonded():
    cell_vectors = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    atoms = [(0.0, 0.0, 9.9, 1), (0.0, 0.0, 10.059, 1)]
    bonds = extrapolate_periodic_bonds(atoms, cell_vectors)
    assert len(bonds) == 0


def test_periodic_atoms_at_0_16_angstrom_should_be_bonded():
    cell_vectors = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    atoms = [(0.0, 0.0, 9.9, 1), (0.0, 0.0, 10.06, 1)]
    bonds = extrapolate_periodic_bonds(atoms, cell_vectors)
    assert len(bonds) == 1
    assert bonds == [(0, 1)]


def test_perodic_h_atoms_at_lte_1_09_angstrom_should_be_bonded():
    cell_vectors = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    atoms = [(0.0, 0.0, 9.0, 1), (0.0, 0.0, 10.09, 1)]
    bonds = extrapolate_periodic_bonds(atoms, cell_vectors)
    assert len(bonds) == 1
    assert bonds == [(0, 1)]


def test_periodic_h_atoms_at_gt_1_09_angstrom_should_not_be_bonded():
    cell_vectors = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    atoms = [(0.0, 0.0, 9.0, 1), (0.0, 0.0, 10.10, 1)]
    bonds = extrapolate_periodic_bonds(atoms, cell_vectors)
    assert len(bonds) == 0


def test_periodic_si_atoms_at_lte_2_77_angstrom_should_be_bonded():
    cell_vectors = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    atoms = [(0.0, 0.0, 9.0, 14), (0.0, 0.0, 11.77, 14)]
    bonds = extrapolate_periodic_bonds(atoms, cell_vectors)
    assert len(bonds) == 1
    assert bonds == [(0, 1)]


def test_periodic_si_atoms_at_gt_2_77_angstrom_should_not_be_bonded():
    cell_vectors = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    atoms = [(0.0, 0.0, 9.0, 14), (0.0, 0.0, 11.78, 14)]
    bonds = extrapolate_periodic_bonds(atoms, cell_vectors)
    assert len(bonds) == 0


def test_periodic_bond_tuples_should_be_sorted_by_atom_index():
    cell_vectors = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    atoms = [(0.0, 0.0, 9.9, 1), (0.0, 0.0, 10.06, 1)]
    bonds = extrapolate_periodic_bonds(atoms, cell_vectors)
    assert bonds == [(0, 1)]
    atoms = [(0.0, 0.0, 10.06, 1), (0.0, 0.0, 9.9, 1)]
    bonds = extrapolate_periodic_bonds(atoms, cell_vectors)
    assert bonds == [(0, 1)]
