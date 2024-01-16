# findTDE
`findTDE` comprises a set of scripts to facilitate easy, high-throughput calculations of threshold displacement energies (TDEs) for materials using ab initio/classical molecular dynamics in `VASP`/`LAMMPS`. The threshold displacement energy is the minimum kinetic energy transfer from incident radiation to a lattice atom that produces a permanent defect. This property is useful for understanding the radiation hardness of a material, and it is a required parameter for binary collision approximation calculations (e.g., SRIM/TRIM).

## Installation
Currently, there is no automatic installation method available. The files may be either downloaded manually or using `git clone`.

## Usage
The find_tde script may be called directly from the command line with several options. The usage may be displayed using the help (-h) option. The convergence mode (-c) determines how subsequent kinetic energy values are chosen for the displacement event, either "standard" (adjust by 5 eV until opposite defect generation is found, then adjust by 1 eV until the TDE is found) or "midpoint" (adjust by 8 eV until opposite defect generation is found, then adjust by half the distance from the current energy to the nearest energy of opposite defect generation). The program selection (-p) chooses whether `VASP` or `LAMMPS` is used for the calculations. If `LAMMPS` is used, the force field file may be chosen (-f).

```bash
find_tde [-h] [-c <standard|midpoint>] [-p <vasp|lammps>] [-f <lmp_ff.type>]
```

The script is currently written to execute via `Slurm` workload manager. This can be adjusted temporarily to execute the appropriate program. 

The script relies on a directory structure. Only the base directory (e.g., "project," can be named anything), main input file ("latt_dirs_to_calc.csv"), inputs directory ("inp"), and perfect supercell directory ("perfect") are required to be made and named as described. `findTDE` should be executed in the "project" directory. Each "displacement" directory, associated "energy" directories, and relevant .csv/.txt files are created by the program.

```
project
│   latt_dirs_to_calc.csv   
│
└───displacement1
│   │   displacement1_data.csv
│   │   displacement1_out.txt
│   │   KE_calcs_list.txt
│   │
│   └───energyA
│   |   │   program_inputs
│   |   │   program_outputs
│   |   │   ...
│   │
│   └───energyB
│   |   │   program_inputs
│   |   │   program_outputs
│   |   │   ...
│   │
│   └───...
│
└───...
│   
└───inp
│   │   INCAR_cgm
│   │   INCAR_md
│   │   KPOINTS
│   │   POSCAR
│   │   POTCAR
│   │   lmp_ff.type
│   │   ...
│   
└───perfect
|   │   OUTCAR
```

The input file "latt_dirs_to_calc.csv" is required to specify the displacement event. This can either be created manually, or by using the `multi_tde.py` accessory script. The heading of this file may be used to describe the file format. The bottom row of the file is read when `findTDE` is executed, and that info is used for that TDE calculation. The first value is a "pseudo" to correspond to the displacement direction, given as an integer number (changes with each unique direction) and either "L" or "S" (describes whether the direction is given as a lattice direction \[u v w\] or spherical direction (rho, phi, theta)). The "atom_type" and "atom_number" detail which atom in the supercell is given the velocity vector to simulate the displacement event (e.g., atom_type: ga and atom_number: 34 corresponds to the 34th Ga atom, as listed in the POSCAR file, being displaced). The initial kinetic energy "ke_i" and cutoff kinetic energy "ke_cut" (stops the program if a defect is not found below this kinetic energy) are then defined. The direction is then defined, either using lattice direction integer notation or spherical coordinate notation (may be floats).

```
########################################
# format of text file
# nL    atom_type    atom_number    ke_i    ke_cut    u    v    w
# n+1S    atom_type    atom_number    ke_i    ke_cut    r    p    t
########################################
```

## Acknowledgements
The findTDE code was developed by Alexander Hauck, Dr. Mia Jin, and Dr. Blair Tuttle at The Pennsylvania State University.
