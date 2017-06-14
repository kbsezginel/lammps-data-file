import pytest
from lammps_data.bonds import calculate_bond_length, BondLengthLimitError


def test_C_H_bond_length_should_return_known_distance():
    atoms = [(0.0, 1.396, 0.0, 6), (0.0, 2.479, 0.0, 1)]
    bond = (0, 1)
    assert round(calculate_bond_length(bond, atoms), 3) == 1.083


def test_bond_length_higher_than_limit_should_raise_an_exception():
    with pytest.raises(BondLengthLimitError):
        atoms = [(5.0, 0.0, 0.0, 6), (0.0, 0.0, 0.0, 1)]
        bond = (0, 1)
        calculate_bond_length(bond, atoms, limit=(0.5, 3))


def test_bond_length_lower_than_limit_should_raise_an_exception():
    with pytest.raises(BondLengthLimitError):
        atoms = [(0.4, 0.0, 0.0, 6), (0.0, 0.0, 0.0, 1)]
        bond = (0, 1)
        calculate_bond_length(bond, atoms, limit=(0.5, 3))
