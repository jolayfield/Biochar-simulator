Installation
============

Requirements
------------

- Python 3.9 or later
- `RDKit <https://www.rdkit.org>`_ ≥ 2023.9
- NumPy ≥ 1.24
- SciPy ≥ 1.10
- NetworkX ≥ 3.1

.. note::
   RDKit is best installed via **conda-forge**.  A pure-PyPI install is
   possible via ``rdkit-pypi`` but the conda route is more reliable and
   gives access to the full toolkit.

Recommended: conda (conda-forge)
---------------------------------

.. code-block:: bash

   conda create -n biochar python=3.11
   conda activate biochar
   conda install -c conda-forge rdkit networkx numpy scipy
   pip install -e ".[dev]"

From PyPI
---------

.. code-block:: bash

   pip install biochar

.. note::
   On some platforms the ``rdkit`` wheel from PyPI may not be available.
   If the install fails, use the conda route above.

Development install
-------------------

Clone the repository and install in editable mode with all dev dependencies:

.. code-block:: bash

   git clone https://github.com/jolayfield/Biochar-simulator.git
   cd Biochar-simulator
   pip install -e ".[dev]"

Verify the installation by running the test suite:

.. code-block:: bash

   pytest tests/ -v
