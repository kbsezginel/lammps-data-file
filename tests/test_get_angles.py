from lammps_data.angles import get_angles


def test_separate_diatomic_molecules_should_have_no_angles():
    bonds = [(0, 1), (2, 3)]
    assert get_angles(bonds) == []

def test_molecule_with_two_bonds_should_have_one_angle():
    bonds = [(0, 1), (1, 2)]
    assert get_angles(bonds) == [(0, 1, 2)]

def test_different_order_of_bond_tuples_should_return_same_order_within_angle_tuples():
    bonds = [(0, 1), (1, 2)]
    assert get_angles(bonds) == [(0, 1, 2)]
    bonds = [(1, 2), (0, 1)]
    assert get_angles(bonds) == [(0, 1, 2)]

def test_different_order_of_bond_tuples_should_return_same_order_of_angle_tuples():
    bonds = [(0, 1), (1, 2), (1, 3)]
    assert get_angles(bonds) == [(0, 1, 2), (0, 1, 3), (2, 1, 3)]
    bonds = [(1, 2), (0, 1), (1, 3)]
    assert get_angles(bonds) == [(0, 1, 2), (0, 1, 3), (2, 1, 3)]

def test_tetrahedral_molecule_should_have_six_angles():
    bonds = [(0, 1), (0, 2), (0, 3), (0, 4)]
    assert get_angles(bonds) == [(1, 0, 2),
                                 (1, 0, 3),
                                 (1, 0, 4),
                                 (2, 0, 3),
                                 (2, 0, 4),
                                 (3, 0, 4)]
