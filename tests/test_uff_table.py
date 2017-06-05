import pytest
from lammps_data.uff_table import UffTable


@pytest.fixture
def uff_table():
    return UffTable()


def test_lookup_of_nonexistent_element_should_raise_an_exception():
    pass


def test_lookup_of_nonexistent_coordination_number_should_raise_an_exception():
    pass


def test_lookup_of_matching_angle_should_return_the_row():
    pass


def test_lookup_of_nonmatching_angle_should_return_the_row():
    pass


def test_lookup_with_angle_multiplicity_using_low_tolerance_should_return_single_row():
    pass


def test_lookup_with_angle_multiplicity_using_high_tolerance_should_raise_an_exception():
    pass
