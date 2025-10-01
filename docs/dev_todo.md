# `findTDE` upcoming work
* add in VASP ML utility
* change how lammps ff is copied -> don't copy just put in inp or specify path elsewhere
* also use option to specify either lammps or vasp execution line
* add variable energy change between runs
* add VASP standard/gamma point versions
* add JUST LAMMPS version
* add verbose/debug option

## multi_tde.py
* move examples to an example doc
* change how directions are defined

## tde_analysis.py
* remove built in requirement for ga 34/etc., require input for specific things
* 1 atom displacement at a time
* maintain plotting defaults, optional user input
* examples with lines for GaN displacements to copy-paste for auto analysis

## defect_analysis.py
* remove built in dphi dtheta for averaging TDEs
* double check defect type checking
* rework defect location analysis

## cli.py
* build out functionalities for single, multiple TDE determinations

## plotting.py
* set defaults based on autodetection, add choice to set options