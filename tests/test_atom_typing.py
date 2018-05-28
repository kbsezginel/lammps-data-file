from itertools import groupby
from lammps_data.crystal import MOF
from lammps_data.bonds import extrapolate_periodic_bonds


def test_hkust1_should_have_correct_atom_types():
    """
    NOT COMPLETED!!!

    hkust1 = MOF('tests/HKUST-1.cif')
    hkust1.calculate_vectors()
    hkust1.bonds = extrapolate_periodic_bonds(hkust1.atoms, hkust1.uc_vectors, RADIUS_BUFFER=0.3)
    hkust1.assign_UFF_atom_types()
    assert hkust1.unique_atom_types.keys == ['C_R', 'Cu', 'H', 'O']
    assert hkust1.atom_types[:5] == ['C', 'C', 'O', 'O', 'Cu']

    n_atom_types = [len(list(group)) for key, group in groupby(sorted(hkust1.atom_types))]
    assert n_atom_types == [100, 50, 150, 124]
    """
    pass
