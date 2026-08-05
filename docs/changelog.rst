Changelog
=========

.. note::

   ``CHANGELOG.md`` at the repository root is the canonical history. This page
   is not a full mirror of it: releases 0.3.0 and 0.4.0 are recorded there and
   are not reproduced here.

0.7.0 (2026-08-05)
------------------

Requirements coverage extended to the four console entry points, which were
excluded on the grounds that argument parsing holds no logic of its own. Writing
``rqm/cli-arguments.md`` showed otherwise: precedence between a loaded config and
the command line, the exit status a pipeline branches on, and consistency between
two individually-valid flags are decided at that layer and nowhere else, and each
had a defect. Fifteen documents, 205 scenarios, every one with a defending test.
``constants.py`` is now the only unspecified module.

.. warning::

   **A config file reloaded with ``--load-config`` was largely discarded. Re-run
   anything built that way.**

   Twelve fields — including ``seed``, ``temperature``, ``feedstock``,
   ``target_num_carbons`` and ``molecule_name`` — were taken from the command
   line unconditionally, and a parsed argument cannot tell a flag the user typed
   from one left at its default. The parser's defaults, not the user's requests,
   overwrote the loaded file. ``H_C_ratio`` and ``O_C_ratio`` survived on a
   different code path, so a reload kept the composition derived from 600 °C
   while dropping the 600 °C that explained it. A run driven entirely by
   command-line flags is unaffected.

.. warning::

   **One change breaks an exit status.** ``biochar-md-setup`` exits 2 when it
   writes no run directories, where it exited 0 — so a pipeline was told the
   stage succeeded and went on to submit an empty directory tree. 2 means here
   what it already means in ``biochar-sweep``. Partial skips stay 0.

- **Explicit functional-group flags replace a loaded group set** rather than
  merging into it. The six counts describe one dictionary between them and
  compete for the same edge sites.
- **``biochar-condense`` refuses a ``--which-repeat`` outside ``--repeats``**,
  before any structure is generated. It previously wrote a ``run_surface.sh``
  pointing at a repeat the condensation run would never create, and nothing
  failed until after the annealing run.
- **``biochar-sweep template`` names only backends that exist.** Its comment
  offered ``gasteiger``, which ``GeneratorConfig`` rejects, and omitted ``ml``.
- **A loaded config's ``box_size`` survives**, being restored to the array the
  field expects rather than left as the list the JSON carries.

0.6.0 (2026-08-04)
------------------

Requirements coverage extended to the five modules that were unspecified:
parameter sweep, pH protonation, valence validation, condensation annealing, and
the QM and ML charge backends. Every module is now specified — fourteen
documents, 188 scenarios, every one with a defending test. Fourteen defects
fixed and five breaking changes; see ``CHANGELOG.md`` for the full list.

.. warning::

   **Two things changed what gets written to disk. Re-run anything whose
   results matter.**

   Structures whose sheets do not kekulise came back from embedding stripped of
   every aromatic flag, so the exported topology carried no improper torsions
   and nothing held those rings planar. A structure that kekulises cleanly is
   unaffected.

   And ``charge_method="ml"`` produces different charges, because the bundled
   model was refitted against the current OPLS typing — at most 0.022 e, RMS
   0.007 e. ``charge_method="opls"`` (the default) and ``"qm"`` are unaffected.

.. warning::

   **Three further changes break existing code or existing outcomes**, all in
   the parameter sweep. A manifest records ``skipped`` where it recorded
   ``failed``; ``on_validation_fail: strict`` raises instead of completing
   silently; and a sweep refuses ``molecule_name`` in ``fixed`` or in an axis,
   since it overrode the per-point name while each point still reported the
   templated one.

- **An aromatic ring heteroatom is no longer reported over-valent.** A furan
  oxygen, a pyrrolic nitrogen and a thiophene sulfur were all judged by summed
  bond orders, so strict mode refused structures containing the ether bridge
  this package builds by design.
- **A titrated structure records its pH**, net charge and sampling seed in the
  files it produces. A protonation state is one draw from an ensemble, and
  nothing on disk said which.
- **Condensation run directories record their provenance**, and their packing
  stage verifies what it placed against what the topology declares.
- **A molecule no longer loses state it was carrying.** Embedding handed back
  the bond-order-rewritten copy it uses internally rather than the caller's
  molecule, and validation answered "does RDKit accept this?" by sanitising the
  caller's molecule — which, for a sheet with no kekulé structure, rewrites its
  aromatic bonds on the way to raising.
- **The bundled ML charge model is no longer a pickle.** It is
  ``charges_gpr_cm5.json``: the fitted hyperparameters and the reference data,
  rebuilt on load and checked against the charges it recorded. This replaces the
  scikit-learn version pin the pickle would otherwise have required.
  ``MLChargeRefinement.model_source`` names which model answered.

0.5.0 (2026-08-03)
------------------

.. warning::

   **Exported topologies changed. Re-run anything whose results matter.**
   Structures exported before this release were missing every 1–4 non-bonded
   interaction and had no term holding aromatic rings planar. Trajectories from
   those topologies are not comparable with trajectories from new ones. The
   structures themselves are unaffected — only the force field applied to them.

- **1–4 interactions are no longer silently dropped.** Topologies declare
  ``nrexcl = 3`` but wrote no ``[ pairs ]`` section; OPLS-AA's ``gen-pairs``
  supplies pair parameters, never the pair list, so every 1–4 interaction was
  excluded and none restored.
- **Aromatic rings carry improper torsions**, so ring carbons are held planar
  rather than free to pyramidalise during dynamics.
- **Atom names above five characters no longer collide** between ``.gro`` and
  ``.itp`` or within a single ``.gro``. Affects systems above ~10,000 atoms of
  one element.
- **The local MD run pipeline could not pass solvation** — ``gmx solvate -p``
  was handed a topology nothing had created.
- **Multi-species ion placement discarded all but the last species.**
- **Annealing follows the pyrolysis temperature** (Wood scaling) instead of
  applying the 400 °C schedule to every structure.
- **Validation reports every geometry error**, not the first three. No pass/fail
  decision changes; ``n_validation_errors`` in sweep manifests becomes the true
  count.
- **Run directories carry** ``run_provenance.json`` recording status, sources,
  net charge, and the annealing schedule used.
- **Strict mode fails on a functional-group shortfall**, not only when a group
  places zero. Sweeps will report more ``fallback`` and fewer ``strict_pass``
  rows; manifest status counts are not comparable with earlier runs.
- **Breaking for direct model callers:**
  ``TemperatureModel.get_valid_range`` now takes the property as its first
  argument. It previously returned H/C's range for every property.
- **New:** ``TemperatureModel.predict_with_evidence`` returns a prediction with
  the observation count, spread, and whether the value was carried in from a
  neighbouring grid point or extrapolated.
- A fallback from a feedstock curve to the pooled curve is reported, and an
  aromaticity extrapolated beyond the fitted H/C range warns.
- **A surface that fell back from amorphous to slit geometry** is no longer
  written with an ``Amorphous surface`` title. ``SurfaceBuilder`` reports the
  realised geometry via ``realised_pore_type`` and ``packing_fell_back``.
- **A carbon skeleton that cannot be built raises** ``SkeletonError`` instead of
  substituting the library's 16-carbon pyrene regardless of the request.
- **Breaking for direct assembler callers:** ``PAHAssembler.generate`` no longer
  accepts ``target_aromaticity``, which was inert. Aromaticity is decided by
  ring topology.
- A ``PAH_LIBRARY`` entry that cannot be parsed or sanitised is reported rather
  than silently dropped.
- Sheets copied by the identical-sheet optimisation no longer share a
  ``composition`` record.
- ``generate_surface`` exposes ``amorphous_fallback``, ``aromaticity_percent``
  and the box paddings.

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

- New :class:`~biochar.workflows.surface_builder.SurfaceBuilder` and
  :class:`~biochar.workflows.surface_builder.SurfaceConfig` for parallel-sheet slit pores
- New :func:`~biochar.workflows.surface_builder.generate_surface` convenience function
- :class:`~biochar.export.gromacs_export.MultiSheetGROWriter` and
  :class:`~biochar.export.gromacs_export.SurfaceTopologyWriter` in ``gromacs_export``
- Identical-sheet optimisation: one ``.itp``, ``count = N`` in ``.top``

**Ether bridge span limit**

- :class:`~biochar.pipeline.biochar_generator.GeneratorConfig` gains ``max_ether_span``
  (default 3 → furan-like 5-membered ring)
- Prevents long-range C–O–C bridges that fold the aromatic sheet into a
  nanotube shape

**H-position optimiser**

- O–H hydrogens are rotated around the C–O bond to minimise steric clashes
- C–H bonds on aromatic sp2 carbons are correctly held in-plane

**Bug fixes**

- Interior H atoms in large hex-lattice sheets are placed correctly
- Pentagon-junction C–H contacts reduced via improved geometry pipeline
