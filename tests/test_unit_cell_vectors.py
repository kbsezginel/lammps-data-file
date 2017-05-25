from lammps_data.crystal import Packing


def test_cubic_unit_cell_vectors():
    # Parameters for ECOLEP_clean from CoRE database and expected results from Avogadro
    expected_vectors = [[18.64000, 0.00000, 0.00000],
                        [0.00000, 18.64000, 0.00000],
                        [0.00000, 0.00000, 18.64000]]
    uc_size = [18.64000, 18.64000, 18.64000]
    uc_angle = [90.00000, 90.00000, 90.00000]
    uc_vectors = Packing.uc_vectors(uc_size, uc_angle)
    # Eliminate rounding errors for floats returned by cos, sin functions from math library
    uc_vectors = [[round(i, 5) for i in j] for j in uc_vectors]
    assert uc_vectors == expected_vectors


def test_orthorhombic_unit_cell_vectors():
    # Parameters for MOGYAI_clean from CoRE database and expected results from Avogadro
    expected_vectors = [[10.55690, 0.00000, 0.00000],
                        [0.00000, 14.84590, 0.00000],
                        [0.00000, 0.00000, 19.06630]]
    uc_size = [10.55690, 14.84590, 19.06630]
    uc_angle = [90.00000, 90.00000, 90.00000]
    uc_vectors = Packing.uc_vectors(uc_size, uc_angle)
    uc_vectors = [[round(i, 5) for i in j] for j in uc_vectors]
    assert uc_vectors == expected_vectors


def test_trigonal_unit_cell_vectors():
    # Parameters for MIHHER_clean from CoRE database and expected results from Avogadro
    expected_vectors = [[32.15390, 0.00000, 0.00000],
                        [16.07695, 27.84609, 0.00000],
                        [16.07695, 9.28203, 26.25355]]
    uc_size = [32.15390, 32.15390, 32.15390]
    uc_angle = [60.00000, 60.00000, 60.00000]
    uc_vectors = Packing.uc_vectors(uc_size, uc_angle)
    uc_vectors = [[round(i, 5) for i in j] for j in uc_vectors]
    assert uc_vectors == expected_vectors


def test_tetragonal_unit_cell_vectors():
    # Parameters for LAVTOT_clean from CoRE database and expected results from Avogadro
    expected_vectors = [[10.26700, 0.00000, 0.00000],
                        [0.00000, 10.26700, 0.00000],
                        [0.00000, 0.00000, 14.46200]]
    uc_size = [10.26700, 10.26700, 14.46200]
    uc_angle = [90.00000, 90.00000, 90.00000]
    uc_vectors = Packing.uc_vectors(uc_size, uc_angle)
    uc_vectors = [[round(i, 5) for i in j] for j in uc_vectors]
    assert uc_vectors == expected_vectors


def test_hexagonal_unit_cell_vectors():
    # Parameters for cm301726k_si_004_clean from CoRE database and expected results from Avogadro
    expected_vectors = [[12.57380, 0.00000, 0.00000],
                        [-6.2869, 10.88923, 0.00000],
                        [0.00000, 0.00000, 14.33400]]
    uc_size = [12.57380, 12.57380, 14.33400]
    uc_angle = [90.00000, 90.00000, 120.00000]
    uc_vectors = Packing.uc_vectors(uc_size, uc_angle)
    uc_vectors = [[round(i, 5) for i in j] for j in uc_vectors]
    assert uc_vectors == expected_vectors


def test_monoclinic_unit_cell_vectors():
    # Parameters for OJAKOA_clean from CoRE database and expected results from Avogadro
    expected_vectors = [[9.01670, 0.00000, 0.00000],
                        [-3.61984, 14.69194, 0.00000],
                        [0.00000, 0.00000, 16.82840]]
    uc_size = [9.01670, 15.13130, 16.82840]
    uc_angle = [90.00000, 90.00000, 103.84100]
    uc_vectors = Packing.uc_vectors(uc_size, uc_angle)
    uc_vectors = [[round(i, 5) for i in j] for j in uc_vectors]
    assert uc_vectors == expected_vectors


def test_triclinic_unit_cell_vectors():
    # Parameters for UVARIT_clean from CoRE database and expected results from Avogadro
    expected_vectors = [[8.40900, 0.00000, 0.00000],
                        [2.09071, 13.32498, 0.00000],
                        [1.04535, 6.66249, 12.41697]]
    uc_size = [8.40900, 13.48800, 14.13020]
    uc_angle = [61.49240, 85.75740, 81.08290]
    uc_vectors = Packing.uc_vectors(uc_size, uc_angle)
    uc_vectors = [[round(i, 5) for i in j] for j in uc_vectors]
    assert uc_vectors == expected_vectors
