.. _findtde_module:

findTDE
=======

Base Module
-----------
.. code:: bash

    Bash script to find TDE value for an atom using VASP AIMD/CGM or LAMMPS MD.
    Syntax: find_tde [-c|p|f|h]
    options:
    c     KE convergence mode: standard or midpoint
    p     Simulation program: vasp or lammps
    f     LAMMPS force field: choose file
    h     Print this Help.

    Usage: find_tde [-h] [-c <standard|midpoint>] [-p <vasp|lammps>] [-f <lmp_ff.type>]

Submodules
----------
.. toctree::
    findtde.defect_analysis
    findtde.io
    findtde.multi_tde
    findtde.plotting
    findtde.tde_analysis
    findtde.utils
    findtde.vasp_vel_write