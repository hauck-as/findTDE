.. _development:

=============================
`findTDE` Development (To-Do)
=============================

-------------
`findTDE` WIP
-------------
* add ability to perform single displacement simulation
* add in VASP ML utility
* change how lammps ff is copied -> don't copy just put in inp or specify path elsewhere
* add variable energy change between runs
* add verbose/debug option
* edit convergence such that kinetic energy is never negative

-----------------
vasp_vel_write.py
-----------------
* change name to just vel_write?

------------
multi_tde.py
------------
* problem appears in LAMMPS execution now on Roar Collab
* move examples to an example doc

---------------
tde_analysis.py
---------------
* remove built in requirement for ga 34/etc., require input for specific things
* 1 atom displacement at a time
* examples with lines for GaN displacements to copy-paste for auto analysis

------------------
defect_analysis.py
------------------
* remove built in dphi dtheta for averaging TDEs
* double check defect type checking
* rework defect location analysis

------
cli.py
------
* build out functionalities for single, multiple TDE determinations
