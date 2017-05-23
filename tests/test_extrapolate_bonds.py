from lammps_data.bonds import extrapolate_bonds


def test_atoms_too_close_should_not_be_bonded():
    atoms = [(0.0, 0.0, 0.0, 1), (0.0, 0.0, 0.159, 1)]
    bonds = extrapolate_bonds(atoms)
    assert len(bonds) == 0

def test_atoms_at_0_16_angstrom_should_be_bonded():
    atoms = [(0.0, 0.0, 0.0, 1), (0.0, 0.0, 0.16, 1)]
    bonds = extrapolate_bonds(atoms)
    assert len(bonds) == 1
    assert bonds == [(0,1)]

def test_h_atoms_at_lte_1_09_angstrom_should_be_bonded():
    atoms = [(0.0, 0.0, 0.0, 1), (0.0, 0.0, 1.09, 1)]
    bonds = extrapolate_bonds(atoms)
    assert len(bonds) == 1
    assert bonds == [(0,1)]

def test_h_atoms_at_gt_1_09_angstrom_should_not_be_bonded():
    atoms = [(0.0, 0.0, 0.0, 1), (0.0, 0.0, 1.10, 1)]
    bonds = extrapolate_bonds(atoms)
    assert len(bonds) == 0

def test_si_atoms_at_lte_2_77_angstrom_should_be_bonded():
    atoms = [(0.0, 0.0, 0.0, 14), (0.0, 0.0, 2.77, 14)]
    bonds = extrapolate_bonds(atoms)
    assert len(bonds) == 1
    assert bonds == [(0,1)]

def test_si_atoms_at_gt_2_77_angstrom_should_not_be_bonded():
    atoms = [(0.0, 0.0, 0.0, 14), (0.0, 0.0, 2.78, 14)]
    bonds = extrapolate_bonds(atoms)
    assert len(bonds) == 0

def test_bond_tuples_should_be_sorted_by_atom_index():
    atoms = [(0.0, 0.0, 0.0, 1), (0.0, 0.0, 0.16, 1)]
    bonds = extrapolate_bonds(atoms)
    assert bonds == [(0,1)]
    atoms = [(0.0, 0.0, 0.16, 1), (0.0, 0.0, 0.0, 1)]
    bonds = extrapolate_bonds(atoms)
    assert bonds == [(0,1)]

def test_ethane_should_have_seven_bonds():

    atoms = [( 1.185080, -0.003838,  0.987524, 1),
             ( 0.751621, -0.022441, -0.020839, 6),
             ( 1.166929,  0.833015, -0.569312, 1),
             ( 1.115519, -0.932892, -0.514525, 1),
             (-0.751587,  0.022496,  0.020891, 6),
             (-1.166882, -0.833372,  0.568699, 1),
             (-1.115691,  0.932608,  0.515082, 1),
             (-1.184988,  0.004424, -0.987522, 1)]

    bonds = extrapolate_bonds(atoms)
    assert len(bonds) == 7

    expected_bonds = [(0, 1),
                      (1, 2),
                      (1, 3),
                      (1, 4),
                      (4, 5),
                      (4, 6),
                      (4, 7)]
    expected_bonds = { frozenset(s) for s in expected_bonds }
    bonds = { frozenset(s) for s in bonds }
    assert len(expected_bonds ^ bonds) == 0
