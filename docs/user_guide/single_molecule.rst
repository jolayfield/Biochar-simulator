Single Molecule Generation
==========================

The primary API for generating a single biochar molecule is
:func:`~biochar.biochar_generator.generate_biochar` (convenience function) or
the lower-level :class:`~biochar.biochar_generator.BiocharGenerator` class.

Generation pipeline
-------------------

.. code-block:: text

   target_num_carbons
         │
         ▼
   Carbon skeleton (PAH library seed + ring-growth)
         │  Parity-aware hex-lattice builder; 100 % aromatic for any size
         ▼
   Oxygen assignment (functional groups dict or O/C-ratio-driven)
         │
         ▼
   Hydrogen assignment (fill valences → target H/C ratio)
         │
         ▼
   3D geometry
         │  ≤ 80 heavy atoms: ETKDGv3 embedding + MMFF94 minimisation
         │  > 80 heavy atoms: 2-D flat sheet + force-field minimisation
         ▼
   OPLS-AA atom typing & partial charges
         │
         ▼
   Validation (composition, geometry, steric clashes)
         │
         ▼
   GROMACS export (.gro / .top / .itp)

Controlling size
----------------

The ``target_num_carbons`` parameter controls skeleton size.
The generator uses the built-in :ref:`PAH library <pah-library>` for sizes
up to 40 C, then grows the skeleton by appending fused hexagonal rings.

.. code-block:: python

   from biochar.biochar_generator import generate_biochar

   # Small molecule — exact match from PAH library
   mol, coords, gro, top, itp = generate_biochar(
       target_num_carbons=24,   # coronene
       output_directory="output",
   )

   # Large sheet — hex-lattice growth path
   mol, coords, gro, top, itp = generate_biochar(
       target_num_carbons=150,
       output_directory="output",
   )

.. _pah-library:

PAH library (seed molecules)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+--------------------------------+-------+------------------+
| Molecule                       | C     | Type             |
+================================+=======+==================+
| benzene                        | 6     | Classic          |
+--------------------------------+-------+------------------+
| naphthalene                    | 10    | Classic          |
+--------------------------------+-------+------------------+
| anthracene / phenanthrene      | 14    | Linear / angular |
+--------------------------------+-------+------------------+
| pyrene                         | 16    | Pericondensed    |
+--------------------------------+-------+------------------+
| chrysene / tetracene           | 18    | Various          |
+--------------------------------+-------+------------------+
| coronene                       | 24    | 7-ring compact   |
+--------------------------------+-------+------------------+
| hex_lattice_28/30/38/40        | 28–40 | Graphene flakes  |
+--------------------------------+-------+------------------+

Controlling composition
-----------------------

H/C and O/C ratios are set as targets with configurable tolerances:

.. code-block:: python

   mol, coords, gro, top, itp = generate_biochar(
       target_num_carbons=80,
       H_C_ratio=0.65,       # high H → lower pyrolysis temperature
       H_C_tolerance=0.10,   # accept ± 10 %
       O_C_ratio=0.20,       # high O → more functional groups
       O_C_tolerance=0.10,
   )

Typical values by pyrolysis temperature:

+-------------+----------+----------+
| Temperature | H/C      | O/C      |
+=============+==========+==========+
| 300–400 °C  | 0.6–0.8  | 0.15–0.25|
+-------------+----------+----------+
| 500–600 °C  | 0.4–0.6  | 0.08–0.15|
+-------------+----------+----------+
| 700–800 °C  | 0.2–0.4  | 0.02–0.08|
+-------------+----------+----------+

Pentagon ring defects
---------------------

Real biochar contains both 6-membered (hexagonal) and 5-membered (pentagonal)
rings due to incomplete graphitization.  The ``defect_fraction`` parameter
controls the probability that each ring added during skeleton growth is a
pentagon instead of a hexagon.

.. code-block:: python

   mol, coords, gro, top, itp = generate_biochar(
       target_num_carbons=80,
       H_C_ratio=0.4,
       O_C_ratio=0.05,
       defect_fraction=0.15,   # ~15 % of rings will be pentagons
       output_directory="output",
       seed=42,
   )

Guidance on ``defect_fraction`` values:

+------------------+------------------------------------------------------+
| Value            | Effect                                               |
+==================+======================================================+
| 0.0 (default)    | Pure hexagonal PAH — flat graphene-like sheet        |
+------------------+------------------------------------------------------+
| 0.05–0.10        | Mild disorder — a few curvature-inducing defects     |
+------------------+------------------------------------------------------+
| 0.10–0.20        | Realistic amorphous biochar topology                 |
+------------------+------------------------------------------------------+
| > 0.30           | Highly disordered — may need more retries to kekulize|
+------------------+------------------------------------------------------+

Pentagon rings introduce curvature and parity constraints.  The generator
automatically handles parity (even/odd node count) to ensure each structure
can be kekulized, retrying with a different sub-seed up to 5 times if needed.

Reproducibility
---------------

Pass an integer ``seed`` to get deterministic output:

.. code-block:: python

   mol, coords, gro, top, itp = generate_biochar(
       target_num_carbons=80,
       seed=42,
   )

Using the class API
-------------------

For more control, use :class:`~biochar.biochar_generator.BiocharGenerator`
directly:

.. code-block:: python

   from biochar.biochar_generator import BiocharGenerator, GeneratorConfig

   config = GeneratorConfig(
       target_num_carbons=100,
       H_C_ratio=0.4,
       O_C_ratio=0.08,
       functional_groups={"phenolic": 4, "carboxyl": 1},
       seed=99,
   )
   gen = BiocharGenerator(config)
   mol, coords, composition = gen.generate()
   gen.print_summary()
   gro, top, itp = gen.export_gromacs("output", "bc_custom")
