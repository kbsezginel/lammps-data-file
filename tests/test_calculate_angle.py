from lammps_data.angles import calculate_angle


def test_linear_points_should_return_180_degrees():
    p1 = [1, 0, 0]
    p2 = [2, 0, 0]
    p3 = [3, 0, 0]
    assert calculate_angle(p1, p2, p3) == 180.0


def test_trigonal_planar_atoms_should_return_120_degrees():
    # Coordinates taken from carbon atoms of a benzene molecule (C-C-C)
    p1 = [1.209, 0.698, 0.000]
    p2 = [1.209, -0.698, 0.000]
    p3 = [0.000, -1.396, 0.000]
    assert round(calculate_angle(p1, p2, p3), 2) == 120.00


def test_tetrahedral_atoms_should_return_109_47_degrees():
    # Coordinates taken from methane molecule (H-C-H)
    p1 = [-0.629118, -0.629118, 0.629118]
    p2 = [0.000000, 0.000000, 0.000000]
    p3 = [-0.629118, 0.629118, -0.629118]
    assert round(calculate_angle(p1, p2, p3), 2) == 109.47


def test_square_planar_atoms_should_return_90_degrees():
    # Coordinates taken from XeF4 molecule (F-Xe-F)
    p1 = [-0.94250, -0.88050, 0.32480]
    p2 = [0.00000, 0.00000, 0.00000]
    p3 = [-0.93700, 0.85680, -0.39600]
    assert round(calculate_angle(p1, p2, p3), 2) == 90.00


def test_water_molecule_should_return_104_5_degrees():
    p1 = [0.757, 0.586, 0.0]
    p2 = [0.000, 0.000, 0.0]
    p3 = [-0.757, 0.586, 0.0]
    assert round(calculate_angle(p1, p2, p3), 2) == 104.51
