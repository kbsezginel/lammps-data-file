from lammps_data.angles import get_angles


def test_separate_diatomic_molecules_should_have_no_angles():
    bonds = [(0, 1), (2, 3)]
    assert get_angles(bonds) == []
