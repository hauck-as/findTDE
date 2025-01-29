# `findTDE` upcoming work
* add in VASP ML utility
* change how lammps ff is copied -> don't copy just put in inp or specify path elsewhere
* also use option to specify either lammps or vasp execution line

## multi_tde.py
* move examples to an example doc
* maybe change from screen? dunno
* change from running to submitting jobs probably
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