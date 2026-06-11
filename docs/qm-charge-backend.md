# QM (1.14\*CM1A) Charge Backend — Implementation & Process Notes

This document records **what was built** (an opt-in LigParGen-style QM charge
method) and **how it was arrived at**, including the research trail and the
judgment calls made along the way. It is meant to be enough for a future
maintainer to trust, extend, or correct the implementation.

## 1. Goal

Add a *secondary, opt-in* partial-charge method to the biochar generator that
uses **the same methodology as LigParGen** — AM1 → CM1A → ×1.14 scaling — while
keeping the static OPLS-AA lookup table as the default behavior. Ship it so the
package still installs from both **pip** and **conda**.

## 2. How charges worked before this change

The generator had two charge paths, selected by `GeneratorConfig.charge_method`:

| method | source | summary |
|--------|--------|---------|
| `"opls"` (default) | `opls_typing.ChargeAssigner` + `constants.OPLS_ATOM_TYPES` | static per-atom-type table lookup, then neutralize |
| `"ml"` | `ml_charges.MLChargeRefinement` | GPR model (trained on OPLS charges) predicts per-atom charges |

This change adds a third: `"qm"`.

## 3. Decision history

1. **Which backend reproduces LigParGen?** LigParGen runs AM1 in BOSS, maps to
   CM1A, and scales neutral-molecule charges by 1.14. Three faithful-ish options
   were considered: (a) the LigParGen web API, (b) a local semiempirical engine
   (MOPAC) + a CM1A implementation, (c) GFN2-xTB (rejected — not CM1A).
2. **Packaging drove the choice.** The user wants pip **and** conda distribution.
   MOPAC is on conda-forge (real per-platform builds, currently 23.x) but is
   **not** a usable PyPI package (PyPI only has an unrelated `mopac 0.9.1`). The
   user chose the **MOPAC route** knowing this asymmetry.
3. **CM1A parameters had to be sourced authoritatively.** The CM1A coefficients
   live in the paywalled Storer et al. 1995 paper and the AMSOL manual. Rather
   than reconstruct them from memory (which risks silently-wrong charges), the
   user downloaded the open-source **AMSOL 7.1** distribution, and the exact
   algorithm + parameters were transcribed from its Fortran source.

## 4. What was implemented

### New module: `biochar/qm_charges.py`

- `QMChargeAssigner` — discovers a `mopac` executable on `PATH` (or `$MOPAC_CMD`),
  writes an AM1 input deck from the molecule + 3D coordinates, runs MOPAC, parses
  the output, applies CM1A, scales by 1.14, and returns `{atom_idx: charge}`.
- `cm1a_from_am1()` — the pure-Python CM1A mapping (no external deps; unit-tested
  in isolation).
- `scale_and_neutralize()` — applies the 1.14 factor and cleans up rounding so the
  total charge equals the molecular formal charge exactly.
- `parse_net_atomic_charges()` / `parse_bond_orders()` — MOPAC `.out` parsers,
  validated against real OpenMOPAC 23.2.5 output.

### Wiring

- `GeneratorConfig.charge_method` now accepts `"qm"` (validation + docstring).
- `BiocharGenerator._assign_opls_properties(mol, coords)` gained the `coords`
  argument (QM needs the geometry) and a `"qm"` branch.
- CLI: `--charge-method` now offers `qm`.

### Packaging

- `pyproject.toml`: a `qm` extra (intentionally empty — MOPAC is not a pip
  package; the extra documents the requirement and lets `pip install biochar[qm]`
  succeed). The backend finds `mopac` on `PATH` at runtime.
- `conda-recipe/meta.yaml`: `mopac >=23` added under `run_constrained` (optional —
  installed only when the user opts in, not forced on every install).

### Tests: `tests/test_qm_charges.py`

CM1A mapping, scaling, and parser tests run with canned data (no binary needed).
An end-to-end test is gated on `shutil.which("mopac")`, and config-validation /
missing-binary error paths are covered.

## 5. The CM1A algorithm (as implemented)

Transcribed from AMSOL 7.1 `new/chgmp1.f` (the CM1 forward mapping) and
`include/PARAMS.i` (the AM1 parameter `DATA` blocks). For each atom *i* with AM1
net (Mulliken) charge `QM(i)` and Mayer bond orders `B_ij`:

```
QC(i) = SCALA1[Z_i]                       (scale)
QD(i) = OFSTA1[Z_i] + Σ_j  B(i,j)·FD(j,i) (offset; B in AMSOL sign convention)
Δ(i)  = QC(i)·QM(i) + QD(i)
Q(i)  = QM(i) + Σ_j  B(i,j)·Δ(j)          ⇒  charge-conserving
```

with `B(i,j) = -B_ij` off-diagonal and `B(i,i) = Σ_j B_ij`, so the result
satisfies `Σ_i Q(i) = Σ_i QM(i)` algebraically. `FD(j,i)` carries the
element-pair offsets: H-on-partner (`HOFFA1`), O-on-partner (`OXOFF1`),
N-on-partner (`ANTOF1`). The AM1 N–C bond uses a special `tanh` switch that
activates for high-order (nitrile) bonds and vanishes for aromatic/amine N–C.

### AM1 parameter values (atomic-number-keyed), verbatim from AMSOL 7.1

```
SCALA1 (scale):   N 0.3846   F 0.1468   P 0.0     S -0.1311  Cl 0.0405  Br 0.1761  I 0.2380
OFSTA1 (offset):  O -0.0283  F 0.0399   S -0.0956 Cl -0.0276 Br -0.0802 I -0.1819
HOFFA1 (H on X):  N 0.0850   O 0.1447   Si 0.0640
OXOFF1 (O on X):  S -0.0600
ANTOF1 (N on X):  C -0.0880 (special N–C term)   O -0.0630
scale factor:     1.14
```

For the biochar element set (C, H, N, O, S) the relevant non-zero parameters are
N/S scales, O/S offsets, H-on-N and H-on-O offsets, and the N–C / N–O terms.
Carbon and H-on-C have no CM1A correction, which is why pure PAHs map almost
directly to scaled AM1 charges.

## 6. Validation performed

- **Parsers** match real MOPAC 23.2.5 output for methane (charges and bond
  orders), and raise clearly on malformed input.
- **Null case:** methane has no applicable CM1A corrections, so CM1A == AM1 net
  charges (confirmed to 1e-9).
- **Hand-trace:** water, computed by hand from real MOPAC AM1 numbers
  (`O = -0.3854`, `H = +0.1927`, `B_OH = 0.963`), gives unscaled CM1A
  `O = -0.7083`, `H = +0.3542`; the code reproduces this exactly. ×1.14 →
  `O = -0.807`, `H = +0.403`, sum 0.
- **Charge conservation** holds for arbitrary inputs (property test).
- **End-to-end:** a 24-carbon sheet with functional groups generated with
  `charge_method="qm"` yields 36 charges summing to 0.0000 with a sensible range
  (≈ ±0.15 e).
- **Known-behavior check:** anilines and phenols come out strongly polar — the
  exact systematic CM1A errors documented in the literature (and the reason the
  later LBCC correction exists). Reproducing these *increases* confidence that
  the mapping is genuinely CM1A.

## 7. Judgment calls (made autonomously, worth knowing)

1. **Single-point AM1 by default (`optimize=False`).** LigParGen optimizes the
   geometry at AM1 first. We instead run a single-point on the generator's
   geometry so charges match the exact coordinates written to GROMACS, it stays
   fast for large sheets, and AM1 re-optimization can't distort a deliberately
   built structure. Empirically the difference is tiny (water: −0.807 vs −0.804).
   `optimize=True` is available for strict LigParGen-style behavior on small
   fragments.
2. **"NET ATOMIC CHARGES", not "MULLIKEN POPULATIONS AND CHARGES".** In the ZDO
   (NDDO) framework AMSOL's `QM` is the diagonal-density charge, which is MOPAC's
   *NET ATOMIC CHARGES* section — not MOPAC's deorthogonalized "Mulliken" section.
   Verified on methane (`-0.2659` = 4 − diagonal population).
3. **Final neutralization.** CM1A conserves charge analytically, but MOPAC prints
   charges/bond orders to finite precision, so a tiny equal-distribution
   correction forces the sum to the exact formal charge.
4. **Base 1.14\*CM1A only; LBCC deferred.** LigParGen's default export is plain
   1.14\*CM1A. The 19 localized bond-charge corrections (1.14\*CM1A-LBCC,
   Vilseck et al. 2017) were extracted but **not** applied; a future `lbcc=True`
   option can add them after scaling.
5. **MOPAC as an optional dependency, not a hard one.** Most users want the
   default OPLS path; forcing a Fortran QM engine on every install is wrong. Hence
   the empty pip extra + conda `run_constrained` + runtime discovery.

## 8. Fidelity caveats

- Charges come from **MOPAC's** AM1; LigParGen uses **BOSS's** AM1. Same
  Hamiltonian, independent codes — values agree closely, not bit-for-bit.
- Bond orders are read from MOPAC's printed table (3 decimals). The resulting
  charge error is ≲ 0.001 e and is cleaned up by neutralization.
- CM1A has documented systematic errors (anilines, phenols, nitro groups); this
  implementation faithfully reproduces them. Use `-LBCC` (future) for those cases.

## 9. Using it

```bash
conda install -c conda-forge mopac          # provide the QM backend
biochar-gen --carbons 24 --charge-method qm --name BCqm
```

```python
from biochar import generate_biochar
r = generate_biochar(target_num_carbons=24, charge_method="qm", molecule_name="BCqm")
```

## 10. References

- Dodda, Cabeza de Vaca, Tirado-Rives & Jorgensen, *Nucleic Acids Res.* **45**,
  W331 (2017) — LigParGen.
- Storer, Giesen, Cramer & Truhlar, *J. Comput.-Aided Mol. Des.* **9**, 87 (1995)
  — Class IV charge models (CM1).
- Vilseck, Tirado-Rives & Jorgensen, *J. Comput. Chem.* (2017) — 1.14\*CM1A-LBCC.
- AMSOL 7.1 (Hawkins, Giesen, Lynch, Chambers, Rossi, Storer, Li, Zhu, Thompson,
  Winget, Lynch, Rinaldi, Liotard, Cramer, Truhlar) — `chgmp1.f`, `PARAMS.i`.
```
