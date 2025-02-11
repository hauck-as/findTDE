.. _findtde_module:

findTDE
=======

Base Module
-----------
.. code-block:: bash
    Help()
    {
    # Display Help
    echo "Bash script to find TDE value for an atom using VASP AIMD/CGM or LAMMPS MD."
    echo
    echo "Syntax: find_tde [-c|p|f|h]"
    echo "options:"
    echo "c     KE convergence mode: standard or midpoint"
    echo "p     Simulation program: vasp or lammps"
    echo "f     LAMMPS force field: choose file"
    echo "h     Print this Help."
    echo
    }

    usage() { echo "Usage: $0 [-h] [-c <standard|midpoint>] [-p <vasp|lammps>]" 1>&2; exit 1; }

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