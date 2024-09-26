# Example for GaN
## Directories
* inp contains the input files for `VASP`/`LAMMPS` that will be copied and edited to perform the displacement simulations
* perfect must contain an OUTCAR file from an energy calculation using the same perfect crystal POSCAR file as in the inp directory. The k-point mesh, potential, and energetically relevant INCAR tags should remain the same between the perfect calculation and the displacement simulations
* 1L\_ga34 contains the directory structure, input files (sans-POTCAR), OUTCAR files, and findTDE output files from a `VASP` displacement simulation of a Ga knockout resulting in a 33 eV TDE

## Base-level Files
* sph\_directions.csv shows a format of spherical directions to calculate using multi\_tde.py
* latt\_dirs\_to\_calc.csv is the general input file for findTDE corresponding the calculation pseudos (e.g., 1L) to the direction as well as specifying the knockout atom and initial/cutoff energies. This can be created manually or using multi\_tde.py
* all\_tde\_data.csv is an output file created with the analysis tools in findTDE combining data from all runs in a directory to be used for further analysis
* gan\_vasp\_psuedo\_keys.csv similarly to latt\_dirs\_to\_calc.csv connects the calculation pseudos to the calculation directions (either spherical or lattice directions) without the additional information is created and used by analysis tools
* find\_tde\_lineplot.png and gan\_ga\_tde.png are output files produced with the analysis tools to show the TDEs for a number of calculations
