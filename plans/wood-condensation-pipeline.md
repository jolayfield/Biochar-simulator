# Plan: Wood et al. 2024 condensation pipeline (parallel construction mode)

**Goal:** add a second construction pipeline to biochar-simulator that reproduces
the Wood, Mašek & Erastova (2024, *Cell Reports Physical Science* 5, 102037)
methodology: pack a library of island-type building blocks into a box and
condense them into an amorphous bulk solid via HTT-scaled simulated annealing,
then expand into an exposed surface, validated on true density / TEM / SASA.

## Locked decisions

1. **Setup-only** — the pipeline writes the packed box + exact Wood `.mdp` files
   + run scripts, ready to run on the user's own GROMACS/HPC. It never invokes
   `gmx`. (Wood's finished models require running ~45 ns × 3 repeats per model;
   that run is the user's step.)
2. **Generate our own building blocks** — do not vendor Wood's published library.
   The existing single-sheet generator (PAH core + aliphatic/O decoration)
   *is* the block source; we select/mix blocks to hit targets.
3. **Parallel mode** — the existing single-sheet generator, `sweep`, `md_setup`,
   and biochar-pfas keep working unchanged. This is additive (new module[s]).

## What "match exactly" means here

- **Exactly matched** (documented protocol): box packing intent, the annealing
  schedule + engine settings, the HTT→temperature/timestep scaling, surface
  creation by z-expansion + semi-isotropic relaxation, and the validation
  observables (true density via free volume, SASA, simulated TEM).
- **Ours, not Wood's** (by decision 2): the actual building-block structures.
  The *architecture* (island: polyaromatic core + alkyl/aryl/O arms) matches;
  the specific molecules are generated here.

## Architecture (layers)

```
 Layer 1  existing single-sheet generator  ─▶ one island building block
          (PAH core + aliphatic arms + O-groups; already built)
 Layer 2  block-library selection          ─▶ a mix of blocks that hits
   NEW    (compose to target H/C, O/C, aromaticity, and true density via
           aromatic-domain-size mix)
 Layer 3  box packing                       ─▶ loose cubic box of N blocks
   NEW    (initial low-density system + combined .top; setup-only)
 Layer 4  Wood annealing setup              ─▶ EM → NVT(hot) → NPT(anneal) →
   NEW    (exact .mdp + run_pipeline.sh,       NPT(final), HTT-scaled, ×3 seeds
           HTT-scaled)
 Layer 5  surface creation setup            ─▶ z-expand + EM + NPT(semi-iso)
   NEW    (second-stage setup)                 .mdp + script
 Layer 6  validation helpers                ─▶ density / SASA / (TEM) analysis
   NEW                                          command scripts for the user's runs
```

## The Wood protocol to reproduce (reference numbers)

**Per model:** pack blocks in a cubic box → EM → NVT → NPT anneal → NPT final,
**×3 repeats** (different start configs) for ergodicity.

HTT-scaled annealing (their Tables 6 & 7):

| Target HTT | NVT 10 ns @ T | NPT anneal 25 ns (T at t = 0 / 10 / 20 / 25 ns) | timestep |
|---|---|---|---|
| 400 °C | 1000 K | 1000 / 1000 / 300 / 300 K | 1 fs |
| 600 °C | 2000 K | 2000 / 2000 / 300 / 300 K | 0.5 fs |
| 800 °C | 3000 K | 3000 / 3000 / 300 / 300 K | 0.5 fs |

(hold high T 0–10 ns → cool to 300 K 10–20 ns → hold 300 K 20–25 ns), then a
**10 ns NPT final** at 300 K / 1 bar.

Engine settings (GROMACS, OPLS-AA; from Ungerer et al. kerogen method):
- **EM**: steepest descent, F_max < 500 kJ/mol/nm, PME, **1.4 nm** vdW, PBC xyz.
- **NVT**: PME, **1 nm** vdW, v-rescale thermostat, τ_t = 0.1 ps.
- **NPT anneal**: PME, **1 nm** vdW, v-rescale 0.1 ps, **Berendsen barostat,
  τ_p = 1 ps, 100 bar, isotropic**, 4-point annealing.
- **NPT final**: PME, **1.4 nm** vdW, v-rescale 0.1 ps, Berendsen 1 ps,
  **1 bar**, 300 K, **2 fs**.

Surface: elongate box in z so layers are separated by **~10 nm** vacuum, EM,
**10 ns NPT semi-isotropic** (xy decoupled from z).

Validation:
- **True density**: GROMACS free-volume tool, He probe **0.13 nm**, averaged
  over the converged 300 K portion; ρ_solid = ρ_system / (1 − V_free,0.13).
- **SASA**: GROMACS `sasa`, N₂ probe **0.18 nm**, normalized nSASA = SASA / (2·A_xy).
- **Simulated TEM**: computeM `ctem`, 200 keV, 0.02 rad aperture, 20 nm defocus.
- **Convergence**: plateaus in RMS, density, box dimensions.

## Phased execution

- **Phase 1 — annealing `.mdp` templates + HTT scaling** (self-contained, exact).
  New module `biochar/condensation.py`: the 4 dry `.mdp` templates parameterized
  by (peak_T, timestep, anneal times), a `HTT → (peak_T, timestep)` mapping
  (anchored on 400/600/800 °C; interpolate/clamp between), and a
  `run_pipeline.sh` renderer chaining EM→NVT→NPT→NPT with ×3 seed repeats.
- **Phase 2 — building-block library selection.** Given targets (H/C, O/C,
  aromaticity, true density), pick core sizes + arm decoration (reuse the
  existing generator) and a block *mix*. Encode the aromatic-domain-size ↔
  true-density relationship as a tuning knob.
- **Phase 3 — box packing (setup-only).** Emit a loose cubic box of the selected
  blocks (via `gmx insert-molecules` staging or a packmol-style input) + the
  combined `.top`. Reuse ideas from the `PreSolvationStage` insertion renderer.
- **Phase 4 — surface creation setup.** z-expansion + EM + 10 ns NPT semi-iso
  `.mdp` + script (second stage, run after condensation).
- **Phase 5 — validation helpers.** Render density / SASA / TEM analysis command
  scripts for the user's produced trajectories.
- **Phase 6 — config + CLI + docs + tests.** A `CondensationConfig`, a
  `biochar-condense` CLI, README section, tests (all setup-only / file-content).

## Open sub-decisions (resolve as we go)

1. **Block parameterization**: reuse biochar-simulator's built-in OPLS-AA typer
   (self-contained, already emits `.itp`) vs. LigParGen/PolyParGen like Wood
   (more faithful, external, offline). Recommend the built-in typer for a
   self-contained setup-only tool; LigParGen as an opt-in.
2. **Curvature**: include heptagons + Wood's ratios (~1 pentagon : 5 hexagons,
   ~1 heptagon : 10 hexagons) in the blocks. Needs a heptagon fuse in
   `carbon_skeleton.py` (pentagons already exist).
3. **HTT interpolation**: Wood only anchored 3 HTTs; decide how to map arbitrary
   requested temperatures (linear interp of peak_T between anchors + clamp).
4. **Non-periodic vs periodic packing**: Wood's bulk is periodic in x,y,z; the
   packed box must be set up for full PBC condensation.

## Refinement (single-molecule block)

Per follow-up decision, the building block is a **single biochar-simulator
molecule** (not a Wood-style mix of many different blocks): generate ONE molecule
at the target composition, pack **N copies** of it into the box, and anneal. This
collapses Phase 2 (library selection) into "get one molecule" and makes Phase 3 a
single `gmx insert-molecules` packing (seeded per repeat for the 3 different
starting configs). Both molecule sources are supported: **generate fresh**
(`generate_and_condense`) or **use an existing `.gro`/`.itp`** (`setup_condensation`).

## Progress

- [x] Methodology extracted from the paper; decisions locked (above).
- [x] **Phase 1 — annealing `.mdp` + HTT scaling** → PR #15 (merged). Exact Wood
      Tables 6/7 numbers, 16 tests.
- [x] **Phase 2+3 — single-molecule packing** → `setup_condensation` /
      `generate_and_condense`: dry `.top` (molecule × N), loose box estimate,
      per-repeat seeded `insert-molecules` packing wired into the run script.
      (Phase 2 "library selection" simplified to the single-molecule block.)
- [x] **Phase 4 — surface creation** → `setup_surface` / `render_surface_script`:
      z-expand the condensed bulk by ~10 nm vacuum (editconf), EM, then 10 ns
      semi-isotropic NPT (xy at 1 bar, z frozen to preserve the gap).
- [x] **Phase 5 — validation helpers** → `render_validation_script` /
      `write_validation_setup`: an `analyze.sh` with true density (freevolume,
      He 0.13 nm; rho_solid = rho_sys/(1-Vf)), normalized SASA (N2 0.18 nm, last
      2 ns), convergence (energy/rms), and the TEM command.
- [x] **Phase 6 — CLI / exports / docs** → `biochar-condense` (`generate` /
      `from-files`), `add_surface_and_validation` convenience, `__init__` exports,
      pyproject entry point, README section. 36 condensation tests total.
- [x] **Curvature (sub-decision 2)** → heptagon fusion (`_fuse_heptagon`) +
      `_choose_ring_size`; `heptagon_fraction` threaded through `PAHAssembler` →
      `GeneratorConfig` → `generate_biochar`; Wood's ratios exported as
      `WOOD_PENTAGON_FRACTION` (2/13) / `WOOD_HEPTAGON_FRACTION` (1/13) and wired
      into `biochar-condense generate --wood-curvature`. 9 heptagon tests added.
