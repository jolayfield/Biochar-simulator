Batch Generation
================

:func:`~biochar.biochar_generator.generate_biochar_series` generates multiple
biochar structures in one call and optionally writes a combined topology
suitable for running a mixed-molecule GROMACS simulation.

Temperature series
------------------

.. code-block:: python

   from biochar.biochar_generator import generate_biochar_series

   configs = [
       {"molecule_name": "BC400", "target_num_carbons": 80,
        "H_C_ratio": 0.65, "O_C_ratio": 0.20, "seed": 1},
       {"molecule_name": "BC600", "target_num_carbons": 100,
        "H_C_ratio": 0.55, "O_C_ratio": 0.12, "seed": 2},
       {"molecule_name": "BC800", "target_num_carbons": 120,
        "H_C_ratio": 0.40, "O_C_ratio": 0.05, "seed": 3},
   ]

   results = generate_biochar_series(
       configurations=configs,
       output_directory="output/temperature_series",
       create_combined_top=True,
   )

   # results["BC400"] -> (gro_path, top_path, itp_path)

This produces individual ``.gro``/``.top``/``.itp`` for each molecule plus
a ``combined.top`` that references all three.

Configuration dict keys
-----------------------

Each entry in ``configurations`` is a plain dict.  Supported keys:

+------------------------+---------+------------------------------------------+
| Key                    | Default | Description                              |
+========================+=========+==========================================+
| ``molecule_name``      | —       | **Required.** Residue name (≤ 5 chars).  |
+------------------------+---------+------------------------------------------+
| ``target_num_carbons`` | 50      | Carbon count.                            |
+------------------------+---------+------------------------------------------+
| ``H_C_ratio``          | 0.5     | Target H/C ratio.                        |
+------------------------+---------+------------------------------------------+
| ``O_C_ratio``          | 0.1     | Target O/C ratio.                        |
+------------------------+---------+------------------------------------------+
| ``aromaticity_percent``| 90.0    | Target aromatic fraction (%).            |
+------------------------+---------+------------------------------------------+
| ``functional_groups``  | None    | Dict of group → count.                   |
+------------------------+---------+------------------------------------------+
| ``defect_fraction``    | 0.0     | Pentagon probability per ring.           |
+------------------------+---------+------------------------------------------+
| ``seed``               | None    | RNG seed.                                |
+------------------------+---------+------------------------------------------+

Running a mixed GROMACS simulation
-----------------------------------

After batch generation, use the combined topology directly:

.. code-block:: bash

   cd output/temperature_series
   gmx grompp -f md.mdp -c bc400.gro -p combined.top -o topol.tpr
   gmx mdrun -deffnm topol

.. note::
   The combined ``.gro`` is not generated automatically — you would
   typically use ``gmx insert-molecules`` or PACKMOL to position multiple
   biochar molecules in a simulation box before running ``grompp``.
