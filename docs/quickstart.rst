Quick Start
===========

From pyrolysis temperature and feedstock
----------------------------------------

Drive composition (H/C, O/C, aromaticity) directly from pyrolysis conditions
using the UC Davis Biochar Database model:

.. code-block:: python

   from biochar import generate_biochar

   mol, coords, gro, top, itp = generate_biochar(
       target_num_carbons=80,
       temperature=600,          # °C
       feedstock="softwood",     # optional; omit for pooled-feedstock curve
       output_directory="output",
       basename="bc600_sw",
       molecule_name="BC600",
       seed=42,
   )

Explicit ``H_C_ratio`` / ``O_C_ratio`` / ``aromaticity_percent`` kwargs
override the model-derived values when provided.

Query reference properties (surface area, pH, CEC, …) for any temperature:

.. code-block:: python

   from biochar import properties

   props = properties(600, feedstock="softwood")
   print(props["H_C_ratio"], props["surface_area_m2_g"])

From the command line::

   biochar-gen --temperature 600 --feedstock softwood --carbons 80 --name BC600 --seed 42

Single molecule
---------------

Generate a biochar molecule and write GROMACS files in one call:

.. code-block:: python

   from biochar import generate_biochar

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

   from biochar import generate_surface

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
