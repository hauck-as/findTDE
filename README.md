# findTDE
`findTDE` comprises a set of scripts to facilitate easy, high-throughput calculations of threshold displacement energies (TDEs) for materials using ab initio/classical molecular dynamics in `VASP`/`LAMMPS`. The threshold displacement energy is the minimum kinetic energy transfer from incident radiation to a lattice atom that produces a permanent defect. This property is useful for understanding the radiation hardness of a material, and it is a required parameter for binary collision approximation calculations (e.g., SRIM/TRIM).

## Installation
Automatic installation is available with `pip install findtde`. The files may also be either downloaded manually or using `git clone`.

We recommend creating a virtual environment through Conda first, then installing `findTDE`. If you plan to use `LAMMPS` to perform calculations, you can install `LAMMPS` in the same environment with Conda.

## Usage
The find_tde script may be called directly from the command line with several options. The usage may be displayed using the help (-h) option.

```bash
$ find_tde -h
usage: find_tde [-h] [-c <standard|midpoint>] [-p <vasp|lammpsish|lammps>] [-d <nL|nS>]
         [-f <lmp_ff.type>] [-x <execution>]

options:
  -h
  -c <standard|midpoint>      Convergence mode, determines how subsequent kinetic energy values
                                are chosen for the displacement event.
                                "standard" adjusts the previous KE value by 5 eV until opposite
                                defect generation behavior is found (defect formed -> defect not
                                formed, or vice versa), then KE is adjusted by 1 eV until the TDE
                                is found.
                                "midpoint" adjusts the previous KE value by 8 eV until opposite
                                defect generation behavior is found, then KE is adjusted by a
                                bisection method until the TDE is found.
  -p <vasp|lammpsish|lammps>  Program used for displacement event simulations. Dictates how and
                                which input files are created.
                                "vasp" selects the ab initio program VASP and only relies on these
                                associated files. Requires POSCAR, POTCAR, KPOINTS, INCAR_md, and
                                INCAR_cgm files in the "inp" directory. A converged OUTCAR of a
                                pristine supercell calculation matching the structure and
                                energetics-related INCAR settings must exist in the "perfect"
                                directory.
                                "lammpsish" performs setup using VASP files, but uses LAMMPS to
                                perform simulations. Requires POSCAR, POTCAR, input.tde, and
                                lmp_ff.type files in the "inp" directory.
                                "lammps" selects the classical molecular dynamics program LAMMPS.
                                Requires read_data.lmp, input.tde, and lmp_ff.type files in the
                                "inp" directory.
  -d <nL|nS>                  Specify the calculation direction pseudonym instead of using the last
                                line in the latt_dirs_to_calc.csv file. Formatted as either "nL" or
                                "nS", where "n" is an integer corresponding to the n-th unique
                                direction calculated, formatted as either a lattice ("L") or
                                spherical ("S") direction.
  -f <lmp_ff.type>            Specifies the name of the forcefield/interatomic potential used for
                                LAMMPS calculations. Can be any filename+extension (e.g., GaN.sw,
                                ngaal.pb). Assumed to exist in the "inp" directory.
  -x <execution>              Dictates the execution line used for the chosen program. Option
                                should be in quotes to prevent word splitting or other issues.
                                Defaults to "srun vasp_std > vasp.out || mpirun vasp_std > vasp.out"
                                for VASP and "srun lmp -in input.tde || mpiexec lmp -in input.tde"
                                for LAMMPS if not specified.
```

The script may be executed directly, executed via a shell script, or executed by submitting a batch script via a workload manager like `Slurm` or `PBS`. By changing the program execution line (-x option), find_tde can be run in serial and `VASP` or `LAMMPS` may be run in parallel.

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

## Citations
If you use `findTDE` in your research, please cite:
* A. S. Hauck, M. Jin, and B. R. Tuttle, “Atomic displacement threshold energies and defect generation in GaN, AlN, and AlGaN: A high-throughput molecular dynamics investigation,” Applied Physics Letters, vol. 124, no. 15, p. 152107, Apr. 2024, [doi: 10.1063/5.0190371](https://doi.org/10.1063/5.0190371).

## Acknowledgements
The `findTDE` code was developed by Alexander Hauck, Dr. Mia Jin, and Dr. Blair Tuttle at The Pennsylvania State University.
