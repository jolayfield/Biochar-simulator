Changelog
=========

0.1.1 (2026-04-18)
-------------------

- Add ``pyproject.toml`` packaging configuration
- Add MIT ``LICENSE`` file
- Add GitHub Actions CI across Python 3.9–3.12
- Full NumPy-style docstrings for all public API classes and functions
- Sphinx documentation

0.1.0 (2026-04-17)
-------------------

**Slit-pore surface generation**

- New :class:`~biochar_simulator.surface_builder.SurfaceBuilder` and
  :class:`~biochar_simulator.surface_builder.SurfaceConfig` for parallel-sheet slit pores
- New :func:`~biochar_simulator.biochar_generator.generate_surface` convenience function
- :class:`~biochar_simulator.gromacs_export.MultiSheetGROWriter` and
  :class:`~biochar_simulator.gromacs_export.SurfaceTopologyWriter` in ``gromacs_export``
- Identical-sheet optimisation: one ``.itp``, ``count = N`` in ``.top``

**Ether bridge span limit**

- :class:`~biochar_simulator.biochar_generator.GeneratorConfig` gains ``max_ether_span``
  (default 3 → furan-like 5-membered ring)
- Prevents long-range C–O–C bridges that fold the aromatic sheet into a
  nanotube shape

**H-position optimiser**

- O–H hydrogens are rotated around the C–O bond to minimise steric clashes
- C–H bonds on aromatic sp2 carbons are correctly held in-plane

**Bug fixes**

- Interior H atoms in large hex-lattice sheets are placed correctly
- Pentagon-junction C–H contacts reduced via improved geometry pipeline
