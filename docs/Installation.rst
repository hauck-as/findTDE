.. _installation:

============
Installation
============

Automatic installation is available via ``pip``.

.. code:: bash

    pip install findtde

We recommend creating a virtual environment through Conda to manage all necessary packages. If you plan to use ``LAMMPS`` to perform calculations, you can `install with Conda <https://docs.lammps.org/Install_conda.html>`_ in the same environment.

.. code:: bash

    conda create -n findtde python=3.12.8
    conda activate findtde
    pip install findtde
    conda install lammps

Alternatively, it can be installed for development purposes by cloning the GitHub repository.

.. code:: bash

    git clone https://github.com/hauck-as/findTDE
    cd findTDE
    pip install -e .
