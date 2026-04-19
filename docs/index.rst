Biochar Simulator
=================

A Python package for generating realistic biochar molecular structures
for GROMACS molecular dynamics simulations.

.. code-block:: python

   from src.biochar_generator import generate_biochar

   mol, coords, gro, top, itp = generate_biochar(
       target_num_carbons=80,
       H_C_ratio=0.5,
       O_C_ratio=0.10,
       functional_groups={"phenolic": 3, "carboxyl": 1},
       output_directory="output",
       seed=42,
   )

**Key features**

- PAH skeletons from 6 to 200+ carbons with exact carbon counts
- H/C and O/C ratio control with configurable tolerances
- Functional groups: phenolic, carboxyl, ether, and more
- Pentagon ring defects for realistic topological disorder
- Slit-pore surface systems (stacked parallel sheets)
- GROMACS-ready ``.gro`` / ``.top`` / ``.itp`` with OPLS-AA force field

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user_guide/single_molecule
   user_guide/functional_groups
   user_guide/surfaces
   user_guide/batch_generation
   user_guide/gromacs_workflow

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/biochar_generator
   api/surface_builder
   api/gromacs_export
   api/opls_typing
   api/constants

.. toctree::
   :maxdepth: 1
   :caption: Development

   changelog

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
