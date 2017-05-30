from lammps_data.crystal import MOF
from lammps_data.bonds import extrapolate_bonds, extrapolate_periodic_bonds
from lammps_data.angles import get_angles
from lammps_data.dihedrals import get_dihedrals

expected_cell_vectors = [[25.832, 0, 0], [0, 25.832, 0], [0, 0, 25.832]]


def test_irmof1():
    irmof1 = MOF('tests/IRMOF-1.cif')

    irmof1.calculate_vectors()
    cell_vectors = [[round(i, 3) for i in j] for j in irmof1.uc_vectors]
    assert cell_vectors == expected_cell_vectors

    irmof1.bonds = extrapolate_bonds(irmof1.atoms)
    assert len(irmof1.bonds) == 488

    irmof1.pbonds = extrapolate_periodic_bonds(irmof1.atoms, irmof1.uc_vectors)
    assert len(irmof1.pbonds) == 512

    irmof1_222 = irmof1.extend_unit_cell(pack=[2, 2, 2])

    assert len(irmof1_222.atoms) == len(irmof1.atoms) * 8

    packed_cell_vectors = [[round(i, 3) for i in j] for j in irmof1_222.uc_vectors]
    expected_packed_cell_vectors = [[round((i * 2), 3) for i in j] for j in irmof1.uc_vectors]
    assert packed_cell_vectors == expected_packed_cell_vectors

    irmof1_222.bonds = extrapolate_bonds(irmof1_222.atoms)
    assert len(irmof1_222.bonds) == 4000

    irmof1_222.pbonds = extrapolate_periodic_bonds(irmof1_222.atoms, irmof1_222.uc_vectors)
    assert len(irmof1_222.pbonds) == 4096
