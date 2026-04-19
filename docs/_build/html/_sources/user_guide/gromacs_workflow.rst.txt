GROMACS Workflow
================

This page describes a complete workflow from structure generation to
production MD using GROMACS.

1. Generate the structure
--------------------------

.. code-block:: python

   from biochar.biochar_generator import generate_biochar

   mol, coords, gro, top, itp = generate_biochar(
       target_num_carbons=80,
       H_C_ratio=0.5,
       O_C_ratio=0.10,
       functional_groups={"phenolic": 3, "carboxyl": 1},
       output_directory="output",
       basename="bc400",
       molecule_name="BC400",
       seed=42,
   )

2. Energy minimisation
-----------------------

Create a minimal ``em.mdp``:

.. code-block:: ini

   ; Steepest descent energy minimisation
   integrator  = steep
   nsteps      = 5000
   emtol       = 100.0
   emstep      = 0.01

   ; Neighbour search
   cutoff-scheme = Verlet
   nstlist       = 10
   ns_type       = grid
   rlist         = 1.0

   ; Electrostatics / VdW
   coulombtype   = PME
   rcoulomb      = 1.0
   vdwtype       = Cut-off
   rvdw          = 1.0

Run minimisation:

.. code-block:: bash

   gmx grompp -f em.mdp -c output/bc400.gro -p output/bc400.top -o em.tpr
   gmx mdrun -v -deffnm em

3. Visualise
-------------

.. code-block:: bash

   # Convert to PDB for VMD or PyMOL
   gmx editconf -f em.gro -o em.pdb

   # Open in VMD
   vmd em.pdb

4. Solvation (optional)
------------------------

Add water for solution-phase simulations:

.. code-block:: bash

   # Expand box
   gmx editconf -f em.gro -o box.gro -c -d 1.0 -bt cubic

   # Solvate
   gmx solvate -cp box.gro -cs spc216.gro -p output/bc400.top -o solvated.gro

   # Add ions if needed
   gmx grompp -f em.mdp -c solvated.gro -p output/bc400.top -o ions.tpr
   gmx genion -s ions.tpr -o ionised.gro -p output/bc400.top -neutral

5. Production MD
-----------------

A basic ``md.mdp`` for NVT at 300 K:

.. code-block:: ini

   integrator  = md
   nsteps      = 500000    ; 1 ns at 2 fs timestep
   dt          = 0.002

   ; Output
   nstxout-compressed = 500
   nstenergy          = 500

   ; Temperature coupling
   tcoupl       = V-rescale
   tc-grps      = System
   tau_t        = 0.1
   ref_t        = 300

   ; Pressure coupling (NPT)
   pcoupl       = Parrinello-Rahman
   tau_p        = 2.0
   ref_p        = 1.0
   compressibility = 4.5e-5

   ; Constraints
   constraints       = h-bonds
   constraint-algorithm = LINCS

.. code-block:: bash

   gmx grompp -f md.mdp -c em.gro -p output/bc400.top -o md.tpr
   gmx mdrun -deffnm md

Force field notes
-----------------

- All files use the **OPLS-AA** force field (``oplsaa.ff/forcefield.itp``).
- Make sure GROMACS can find the force field — it ships with GROMACS and
  is on the default search path.  If ``grompp`` complains, copy the
  ``oplsaa.ff`` directory into your working directory.
- For slit-pore surfaces the ``.top`` uses a ``count = N`` molecule entry
  (identical sheets) or separate entries per sheet type (distinct sheets).
