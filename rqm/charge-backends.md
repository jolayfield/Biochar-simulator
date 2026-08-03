# Feature: Charge Backends <!-- rq-fb7f157f -->

Three things can decide a structure's partial charges. `charges/qm_charges.py` runs an AM1
calculation through an external MOPAC binary and maps it to LigParGen's 1.14×CM1A;
`charges/ml_charges.py` predicts them from a bundled Gaussian-process model; and the static OPLS
table in `pipeline/opls_typing.py` is what happens when neither is asked for.

They are not interchangeable, and the differences are the point. One needs a compiled Fortran
program that is not pip-installable. One needs a pickle trained by a particular scikit-learn on
particular reference data. One is a lookup table. A caller picking between them is picking a
provenance, and the failure mode throughout is a structure whose charges came from somewhere other
than where the caller believes.

So the requirements here are about **saying which backend answered, and under what conditions**:
that an absent binary is a refusal rather than a fallback, that a substituted model announces
itself, that a model trained under a different library says so, and that an empirical factor
derived for neutral molecules says so when it is applied to an ion.

## Feature API <!-- rq-f03c2518 -->

- `QMChargeAssigner.assign(mol, coords, ...) -> Dict[int, float]` <!-- rq-4e3bf53e -->
  - AM1 through MOPAC, mapped to CM1A, scaled by 1.14, corrected to the formal charge.

- `cm1a_from_am1(charges, bond_orders, atomic_numbers) -> List[float]` <!-- rq-a27a7e8c -->
  - The Class IV mapping from AM1 Mulliken charges and bond orders, transcribed from AMSOL 7.1.

- `scale_and_neutralize(cm1a_charges, total_charge=0.0, scale=1.14) -> List[float]` <!-- rq-f88930c4 -->
  - Applies LigParGen's condensed-phase factor and forces the sum to the formal charge.

- `parse_net_atomic_charges(out_text)`, `parse_bond_orders(out_text)` <!-- rq-dc7908c7 -->
  - Readers for the two MOPAC output tables the mapping needs.

- `QMChargeError` <!-- rq-0be99a80 -->
  - Raised when the QM path cannot run — most often because MOPAC is not installed.

- `MLChargeRefinement(model_path=None)` <!-- rq-1594d21a -->
  - `refine(mol, atom_types) -> Dict[int, float]` predicts a charge per atom.
  - Reports which model answered, and whether the library loading it is the one that wrote it.

- `MLChargeRefinement.train_and_save(...)`, `build_and_save_bundled_model(...)` <!-- rq-0c2933ce -->
  - Retraining on DFT-derived references, and rebuilding the bundled artifact.

## A Backend That Cannot Run Refuses <!-- rq-b8eee7a8 -->

MOPAC is a compiled Fortran program, not a Python dependency, so `charge_method="qm"` is a request
that can fail on a perfectly good install. It fails loudly: `QMChargeError` naming the binary it
looked for and how to get it.

Falling back to the OPLS table would be worse than the error. The caller asked for QM charges
because the table was not good enough for what they are doing, and a structure exported with table
charges under a QM request is indistinguishable from one that got them.

```gherkin
Feature: Refuse rather than substitute a different backend

  @rq-6b61adce
  Scenario: A missing MOPAC binary is reported
    Given no MOPAC binary on the path
    When QM charges are requested
    Then the request is refused, naming the binary and how to install it

  @rq-44546ffc
  Scenario: An unrecognised charge method is refused at construction
    Given a charge method that is not one of the three
    When a generator configuration is built
    Then the configuration is refused
```

## The Charges Sum to the Formal Charge <!-- rq-c2f9b1a9 -->

Whatever produced them, the per-atom charges must sum to the molecule's formal charge. A topology
whose charges sum to something else is a system with a net charge nobody declared, and the
neutralising step downstream balances against the declared figure.

CM1A conserves charge by construction, so after scaling the sum is the scale times the formal
charge; the residual correction only cleans up the finite precision MOPAC prints at.

```gherkin
Feature: Conserve the molecular charge

  @rq-ff6dbfa6
  Scenario: The CM1A mapping conserves the total charge
    Given AM1 charges and bond orders for a molecule
    When they are mapped to CM1A
    Then the mapped charges sum to the same total

  @rq-c49c52b0
  Scenario: Scaling ends at the formal charge
    Given a set of CM1A charges and a formal charge
    When they are scaled and corrected
    Then they sum to the formal charge
```

## The 1.14 Factor Was Derived for Neutral Molecules <!-- rq-22fa0ea0 -->

LigParGen's 1.14 is an empirical condensed-phase correction fitted against neutral organic liquids.
Applied to an ion it is an extrapolation, and the arithmetic makes that visible: scaling a molecule
of formal charge −3 gives a sum of −3.42, and the correction that brings it back to −3 does so by
shifting **every** atom by the same amount. That uniform smear is neither CM1A nor LigParGen; it is
an artefact of applying a neutral-molecule factor to something that is not one.

This is not a corner case. Combining a requested pH with the ML backend is refused, and the refusal
recommends the QM backend by name — so the charged-molecule path is the one callers are sent down.

The behaviour stays: changing an ion's charges silently would be a worse answer than an
extrapolated one. What changes is that it says so, naming the formal charge, the factor, and what
the correction did.

```gherkin
Feature: Say when the factor is being extrapolated

  @rq-817541fc
  Scenario: Scaling a charged molecule reports the extrapolation
    Given CM1A charges for a molecule with a non-zero formal charge
    When they are scaled and corrected
    Then a warning names the formal charge and that the factor is a neutral-molecule parameterisation

  @rq-2d0042a6
  Scenario: Scaling a neutral molecule reports nothing
    Given CM1A charges for a molecule with no formal charge
    When they are scaled and corrected
    Then no extrapolation warning is issued
```

## The Refiner Says Which Model Answered <!-- rq-b1e463e8 -->

The ML backend has two models: the bundled artifact trained on OPLS reference charges, and a
fallback built at run time when that artifact cannot be found. They are not the same model and do
not give the same charges.

The substitution is currently a log line, which is to say invisible: a caller who asked for the
bundled model and got the fallback has nothing in their hands that differs. The refiner reports
which one answered, as a value rather than as a side effect.

```gherkin
Feature: Name the model behind the charges

  @rq-c202e1dc
  Scenario: The refiner reports the bundled model when it loaded
    Given the bundled model artifact
    When a refiner is constructed over it
    Then it reports that the bundled model answered

  @rq-3cf63133
  Scenario: A substituted fallback model announces itself
    Given a model path that does not exist
    When a refiner is constructed over it
    Then it reports the fallback and warns that the charges are not the bundled model's
```

## A Model Written by a Different Library Says So <!-- rq-39798517 -->

The bundled model is a pickled scikit-learn pipeline, and a pickle is only as portable as the
library that wrote it. Loading one under a different scikit-learn produces that library's own
`InconsistentVersionWarning` — which says, correctly, that it "might lead to breaking code or
invalid results" — and then the charges are used anyway.

That warning is the wrong shape for this. It names two version strings and a class, not the model
file, not what the model is for, and not what a reader should do about it; and it arrives from a
dependency, so it reads as library noise rather than as a statement about these charges. The
package constraint is a floor with no ceiling, so the mismatch is expected rather than exceptional.

The version that wrote the artifact is recorded beside it and checked on load, and a mismatch is
reported in terms of the charges it affects.

```gherkin
Feature: Check the library against the artifact

  @rq-1e274fe6
  Scenario: A model written by a different scikit-learn is reported
    Given a bundled model recorded as written by a different scikit-learn version
    When a refiner is constructed
    Then a warning names the model, both versions, and that the charges may be invalid

  @rq-75d3624b
  Scenario: A matching version is not reported
    Given a bundled model recorded as written by the running scikit-learn version
    When a refiner is constructed
    Then no version warning is issued
```

## Every Atom Gets a Charge <!-- rq-b7d22b5f -->

A refinement returns one charge per atom, hydrogens included. A missing entry is not a neutral
atom; it is an atom the exporter writes as zero, which is a different molecule.

```gherkin
Feature: Cover the whole molecule

  @rq-d3ccbbe3
  Scenario: Refinement returns a charge for every atom
    Given a molecule and its atom types
    When its charges are refined
    Then every atom index appears in the result
```

## QM Charges Belong to the Exported Geometry <!-- rq-17279822 -->

LigParGen optimises at the AM1 level before computing charges. This module does not, by default:
it runs a single-point calculation on the coordinates the generator produced, so the charges belong
to the geometry actually written to the GROMACS files rather than to a re-optimised one the user
never sees. Re-optimisation is available and is a caller's decision.

The difference from LigParGen is deliberate and is recorded, because the two produce different
numbers and a reader comparing against the server needs to know which they have.

```gherkin
Feature: Charge the structure that was exported

  @rq-1e6d149c
  Scenario: Charges are computed on the given coordinates by default
    Given a molecule with coordinates
    When QM charges are assigned without requesting optimisation
    Then the calculation is a single point rather than a geometry optimisation
```

## Cross-references <!-- rq-47323ea3 -->

- The refusal of a pH request combined with the ML backend, and the choice between the three
  methods, are specified in `generation-config.md`.
- The static OPLS charges these two back onto, and the atom types they key from, are specified in
  `opls-typing.md`.
- The total-charge statement a topology makes, and the step that neutralises it, are specified in
  `gromacs-export.md`.
