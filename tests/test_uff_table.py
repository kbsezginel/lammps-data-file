import pytest
from lammps_data.uff_table import UffTable, NonExistentElementError, \
                                  NonExistentCoordinationError, NonExistentAngleError, \
                                  LookupMultiplicityError


@pytest.fixture
def uff_table():
    return UffTable()


def test_lookup_of_nonexistent_element_should_raise_an_exception(uff_table):
    with pytest.raises(NonExistentElementError):
        uff_table.lookup('M', 2, 180)


def test_lookup_of_nonexistent_coordination_number_should_raise_an_exception(uff_table):
    with pytest.raises(NonExistentCoordinationError):
        uff_table.lookup('C', 5, 180)


def test_lookup_of_nonexistent_angle_should_raise_an_exception(uff_table):
    with pytest.raises(NonExistentAngleError):
        uff_table.lookup('C', 4, 100, tolerance=0)


def test_lookup_of_matching_angle_should_return_the_row(uff_table):
    expected_row = dict(type='C_1', bond=0.706, angle=180.0, distance=3.851, energy=0.105, scale=12.73, charge=1.912)
    assert uff_table.lookup('C', 2, 180, tolerance=0) == expected_row


def test_lookup_of_angle_within_tolerance_should_return_the_row(uff_table):
    expected_row = dict(type='C_1', bond=0.706, angle=180.0, distance=3.851, energy=0.105, scale=12.73, charge=1.912)
    assert uff_table.lookup('C', 2, 175, tolerance=10) == expected_row


def test_lookup_with_angle_multiplicity_using_low_tolerance_should_return_single_row():
    assert len(uff_table.lookup('O', 4, 104.51, tolerance=0)) == 1


def test_lookup_with_angle_multiplicity_using_high_tolerance_should_raise_an_exception():
    with pytest.raises(LookupMultiplicityError):
        uff_table.lookup('O', 4, 104.51, tolerance=50)
