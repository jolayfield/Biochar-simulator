Porous Surface Generation
=========================

The surface builder creates slit-pore systems consisting of multiple
parallel graphene-like sheets separated by a controllable pore gap and
exports them as a single GROMACS simulation box.

.. code-block:: text

   SurfaceConfig
         │
         ▼
   Generate N sheets (each via single-molecule pipeline)
         │  Identical sheets: generate once, deep-copy remainder
         │  Distinct sheets:  generate each independently
         ▼
   Flatten each sheet to xy plane (SVD best-fit rotation)
         │
         ▼
   Stack along z: sheet i at z = i × (pore_diameter + 3.4 Å)
         │
         ▼
   Compute periodic box (bounding box + padding)
         │
         ▼
   Centre system in box
         │
         ▼
   GROMACS export
         │  .gro  — all N sheets, contiguous atom numbering
         │  .itp  — one file (identical) or one per unique sheet
         └─ .top  — forcefield + itp includes + [ molecules ] count

Symmetric slit pore
-------------------

Two chemically identical sheets — only one sheet is generated and
deep-copied.  A single ``.itp`` is produced with ``count = 2`` in the
``.top``.

.. code-block:: python

   from biochar.biochar_generator import generate_surface

   sheets, gro, top, itps = generate_surface(
       target_num_carbons=50,
       H_C_ratio=0.3,
       O_C_ratio=0.05,
       functional_groups={"phenolic": 2},
       pore_diameter=10.0,   # Å — gap between vdW surfaces
       num_sheets=2,
       output_directory="output",
       basename="slit_pore",
       seed=42,
   )

The pore diameter is the free gap between the inner van-der-Waals
surfaces of the two sheets.  The centre-to-centre separation is
``pore_diameter + 3.4 Å`` (graphite interlayer spacing).

Asymmetric pore
---------------

Use ``sheet_overrides`` to give each wall different chemistry.  One
``.itp`` is written per unique sheet.

.. code-block:: python

   sheets, gro, top, itps = generate_surface(
       pore_diameter=8.0,
       num_sheets=2,
       sheet_overrides=[
           # Wall 1: phenolic-rich, smaller sheet
           {"target_num_carbons": 40, "functional_groups": {"phenolic": 3}},
           # Wall 2: carboxyl-rich, larger sheet
           {"target_num_carbons": 60, "functional_groups": {"carboxyl": 2}},
       ],
       output_directory="output",
       basename="asymmetric_pore",
   )

Multi-sheet stacking
--------------------

More than two sheets create multiple consecutive pores with equal gaps:

.. code-block:: python

   sheets, gro, top, itps = generate_surface(
       target_num_carbons=50,
       pore_diameter=12.0,
       num_sheets=4,          # three pore gaps
       output_directory="output",
       basename="multilayer",
       seed=7,
   )

Using the class API
-------------------

For programmatic control, use :class:`~biochar.surface_builder.SurfaceBuilder`
directly:

.. code-block:: python

   from biochar.surface_builder import SurfaceBuilder, SurfaceConfig

   config = SurfaceConfig(
       target_num_carbons=60,
       functional_groups={"phenolic": 2, "ether": 1},
       pore_diameter=10.0,
       num_sheets=2,
       box_padding_xy=1.5,   # nm
       box_padding_z=1.0,    # nm
       seed=42,
   )

   builder = SurfaceBuilder(config)
   sheets, box_nm = builder.build()

   print(f"Box: {box_nm} nm")
   for i, sheet in enumerate(sheets):
       z_vals = sheet.coords[:, 2]
       print(f"Sheet {i+1}: z = {z_vals.mean():.2f} Å")

   gro, top, itps = builder.export_gromacs("output", "slit_pore")

Output files
------------

+---------------------------+--------------------------------------------------+
| File                      | Contents                                         |
+===========================+==================================================+
| ``<basename>.gro``        | All sheets; residue numbers 1, 2, ... N          |
+---------------------------+--------------------------------------------------+
| ``<basename>.top``        | Force field include, itp includes, [molecules]   |
+---------------------------+--------------------------------------------------+
| ``<basename>_sheet.itp``  | Identical-sheet topology (count = N in .top)     |
+---------------------------+--------------------------------------------------+
| ``<basename>_sht1.itp``   | Sheet 1 topology (distinct-sheet case)           |
+---------------------------+--------------------------------------------------+
| ``<basename>_sht2.itp``   | Sheet 2 topology (distinct-sheet case)           |
+---------------------------+--------------------------------------------------+
