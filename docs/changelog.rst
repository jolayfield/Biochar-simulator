Changelog
=========

0.2.0 (2026-06-01)
-------------------

- **Temperature × feedstock composition model** — ``GeneratorConfig`` and
  ``generate_biochar()`` now accept ``temperature`` (°C) and ``feedstock``
  (e.g. ``"softwood"``, ``"grass"``) to derive H/C, O/C, and aromaticity
  targets from the UC Davis Biochar Database.  Explicit ratio kwargs still
  override the derived values.
- **CLI** — ``biochar-gen --temperature 600 --feedstock softwood`` now works;
  ``--hc-ratio`` / ``--oc-ratio`` / ``--aromaticity`` default to ``None``
  (derived from the model when ``--temperature`` is given).
- **New public API** — ``biochar.properties(temperature, feedstock=None)``
  returns the full reference property table (pH, surface area, CEC, …);
  ``biochar.VALID_FEEDSTOCKS`` lists accepted feedstock names;
  ``biochar.TemperatureModel`` is the underlying model class.
- **Documentation** — API reference page for ``temperature_model``;
  temperature/feedstock example added to the Quick Start guide.

0.1.4 (2026-05-31)
-------------------

- **ML-based partial charge refinement** — opt-in ``charge_method="ml"`` using a
  bundled Gaussian-process model trained on OPLS reference charges (issue #4)
- CI and Read the Docs updated to install scikit-learn for the ``ml`` extra

0.1.3 (2026-05-29)
-------------------

- **Amorphous porous packing** — ``pore_type="amorphous"`` (issue #1)
- **S-doping** — thiol and thioether functional groups (issue #3)
- **Ring-substituting nitrogen** — pyridinic / pyrrolic / graphitic (issue #2)
- Expanded test coverage (~83 %)

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

- New :class:`~biochar.surface_builder.SurfaceBuilder` and
  :class:`~biochar.surface_builder.SurfaceConfig` for parallel-sheet slit pores
- New :func:`~biochar.biochar_generator.generate_surface` convenience function
- :class:`~biochar.gromacs_export.MultiSheetGROWriter` and
  :class:`~biochar.gromacs_export.SurfaceTopologyWriter` in ``gromacs_export``
- Identical-sheet optimisation: one ``.itp``, ``count = N`` in ``.top``

**Ether bridge span limit**

- :class:`~biochar.biochar_generator.GeneratorConfig` gains ``max_ether_span``
  (default 3 → furan-like 5-membered ring)
- Prevents long-range C–O–C bridges that fold the aromatic sheet into a
  nanotube shape

**H-position optimiser**

- O–H hydrogens are rotated around the C–O bond to minimise steric clashes
- C–H bonds on aromatic sp2 carbons are correctly held in-plane

**Bug fixes**

- Interior H atoms in large hex-lattice sheets are placed correctly
- Pentagon-junction C–H contacts reduced via improved geometry pipeline
