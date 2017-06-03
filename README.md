# lammps-data-file
Topology analysis and force field assignment for [Lammps][lammps] simulation package.

## Installation
Python 3.5 is required with additional libraries listed in requirements.txt file.

You can install by cloning the repository and installing requirements as follows:

```bash
git clone https://github.com/WilmerLab/lammps-data-file.git`
cd lammps-data-file
pip install -r requirements.txt
```
## File Input/Output
File loading and saving are handled using [IPMOF][ipmof] library which utilizes [ASE][ase] python library. Any format supported by ASE can be used to load and save files.

## Topology
To determine topology, bonds, angles, and dihedrals can be calculated for a given structure.

### Bonds
Bonds are extrapolated using covalent radii calculated by Pyykkö et al. and a skin distance of 0.45 Å.

For periodic structures bonds that asdasd the boundary can be extrapolated using unit cell vectors. This information is required for Lammps simulations.

### Angles
Angles are determined using bond list. For each bond pair that shares an atom a 3-element tuple is listed with common atom index in the middle, small index first and large index last.

### Dihedrals
Dihedrals are determined using bond list. For a given bond two other bonds that share one of the atoms in that bond are selected. If those bonds also share an atom dihedral is discarded as that means it's a closed structure. Otherwise a 4-element tuple is listed as using two atoms of the original bond in the middle in same order and adding smaller index first and larger index last.

## Force field
A combination of UFF and UFF4MOF force field is used for assigning force field parameters.

## Lammps
Lammps input files and data files can be created with UFF parameters for a given structure.

[lammps]: http://lammps.sandia.gov/
[ipmof]: https://github.com/kbsezginel/IPMOF
[ase]: https://wiki.fysik.dtu.dk/ase/
