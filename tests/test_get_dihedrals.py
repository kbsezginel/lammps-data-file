from lammps_data.dihedrals import get_dihedrals


def test_triangle_molecule_should_have_no_dihedrals():
    bonds = [(0, 1), (0, 2), (1, 2)]
    assert get_dihedrals(bonds) == []


def test_four_atom_molecule_should_have_one_dihedral():
    bonds = [(0, 1), (1, 2), (2, 3)]
    assert get_dihedrals(bonds) == [(0, 1, 2, 3)]


def test_different_order_of_bond_tuples_should_return_same_dihedral():
    bonds = [(0, 1), (1, 2), (2, 3)]
    assert get_dihedrals(bonds) == [(0, 1, 2, 3)]
    bonds = [(2, 3), (1, 2), (0, 1)]
    assert get_dihedrals(bonds) == [(0, 1, 2, 3)]


def test_different_order_of_bond_tuples_should_return_same_list_of_dihedral_tuples():
    bonds = [(0, 2), (1, 2), (2, 3), (3, 4), (3, 5)]
    assert get_dihedrals(bonds) == [(0, 2, 3, 4), (0, 2, 3, 5), (1, 2, 3, 4), (1, 2, 3, 5)]
    bonds = [(1, 2), (0, 2), (2, 3), (3, 5), (3, 4)]
    assert get_dihedrals(bonds) == [(0, 2, 3, 4), (0, 2, 3, 5), (1, 2, 3, 4), (1, 2, 3, 5)]
