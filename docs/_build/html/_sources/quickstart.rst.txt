Quick Start
===========

Single molecule
---------------

Generate a biochar molecule and write GROMACS files in one call:

.. code-block:: python

   from src.biochar_generator import generate_biochar

   mol, coords, gro, top, itp = generate_biochar(
       target_num_carbons=80,
       H_C_ratio=0.5,
       O_C_ratio=0.10,
       output_directory="output",
       basename="bc400",
       molecule_name="BC400",
       seed=42,
   )

This produces ``output/bc400.gro``, ``output/bc400.top``, and
``output/bc400.itp``.

With functional groups
----------------------

Specify exact functional group counts instead of relying on the O/C ratio:

.. code-block:: python

   mol, coords, gro, top, itp = generate_biochar(
       target_num_carbons=60,
       functional_groups={"phenolic": 3, "carboxyl": 1, "ether": 2},
       output_directory="output",
       basename="bc_fg",
       seed=42,
   )

With pentagon defects
---------------------

Add topological disorder by inserting 5-membered rings during skeleton growth:

.. code-block:: python

   mol, coords, gro, top, itp = generate_biochar(
       target_num_carbons=80,
       H_C_ratio=0.4,
       O_C_ratio=0.05,
       defect_fraction=0.15,   # ~15 % of rings will be pentagons
       output_directory="output",
       seed=42,
   )

Slit-pore surface
-----------------

Stack two parallel sheets separated by a 10 Å pore:

.. code-block:: python

   from src.biochar_generator import generate_surface

   sheets, gro, top, itps = generate_surface(
       target_num_carbons=50,
       functional_groups={"phenolic": 2},
       pore_diameter=10.0,
       num_sheets=2,
       output_directory="output",
       basename="slit_pore",
       seed=42,
   )

Run in GROMACS
--------------

.. code-block:: bash

   gmx grompp -f em.mdp -c output/bc400.gro -p output/bc400.top -o em.tpr
   gmx mdrun -v -deffnm em

See :doc:`user_guide/gromacs_workflow` for a complete energy-minimisation
and MD workflow.
