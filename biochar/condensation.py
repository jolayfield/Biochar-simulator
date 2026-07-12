"""
Wood et al. 2024 condensation-annealing setup (parallel construction mode).

Reproduces the simulated-annealing protocol from Wood, Mašek & Erastova,
*Cell Reports Physical Science* 5, 102037 (2024): island-type building blocks
are packed into a periodic box and condensed into an amorphous bulk solid via an
HTT-scaled heat/cool cycle, then (later phases) expanded into an exposed surface.

This module is **Phase 1**: the exact GROMACS ``.mdp`` templates, the
HTT → (peak temperature, timestep) mapping, and a run-script renderer that
chains the four dry stages with the ×3 repeats Wood used for ergodicity.

**Setup-only:** it writes files; it never invokes ``gmx``. Producing a finished
condensed model means running these ~45 ns (×3) simulations on your own GROMACS
build, which is the user's step.

Protocol per model (Wood Tables 6 & 7 + experimental procedures):

    EM (steepest descent) → NVT 10 ns @ peak_T → NPT anneal 25 ns
    (hold peak_T 0–10 ns, cool to 300 K 10–20 ns, hold 300 K 20–25 ns,
     Berendsen 100 bar) → NPT final 10 ns (300 K, 1 bar, 2 fs)

HTT scaling (their anchors):

    400 °C → 1000 K, 1.0 fs      600 °C → 2000 K, 0.5 fs      800 °C → 3000 K, 0.5 fs

The NVT + anneal stages use the HTT-scaled timestep; the final equilibration is
always 2 fs (Wood: "for all systems ... a 2 fs timestep").

Note: the paper does not state the constraint scheme explicitly. The 0.5–1 fs
annealing timesteps imply unconstrained bonds (small dt resolves bond vibrations
at high T), so the annealing stages here use ``constraints = none`` and the 2 fs
final stage uses ``constraints = h-bonds``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


class CondensationError(Exception):
    """Raised for invalid condensation-setup inputs."""


# --------------------------------------------------------------------------- #
# HTT → (peak temperature, timestep) scaling
# --------------------------------------------------------------------------- #
# Wood anchored three HTTs; peak_T is interpolated linearly between/clamped
# outside these, and the timestep drops to 0.5 fs once the peak exceeds ~1500 K
# (i.e. above the 400 °C anchor), matching their Table 6/7 choices.
_HTT_ANCHORS = {
    400.0: (1000.0, 1.0),
    600.0: (2000.0, 0.5),
    800.0: (3000.0, 0.5),
}


@dataclass(frozen=True)
class AnnealSpec:
    """The two HTT-dependent knobs of the Wood annealing protocol."""

    peak_T_K: float       # NVT + anneal high temperature
    timestep_fs: float    # NVT + anneal timestep (final stage is always 2 fs)

    def __post_init__(self):
        if self.peak_T_K < 300.0:
            raise CondensationError("peak_T_K must be >= 300 K")
        if self.timestep_fs <= 0.0:
            raise CondensationError("timestep_fs must be > 0")


def anneal_spec_for_htt(htt_c: float) -> AnnealSpec:
    """Map a pyrolysis HTT (°C) to Wood's annealing peak temperature + timestep.

    Anchored on 400/600/800 °C; peak_T is linearly interpolated between anchors
    and clamped outside them. Timestep is 1.0 fs at/below the 400 °C anchor and
    0.5 fs above it (Wood used 0.5 fs for the hotter 600/800 °C systems).
    """
    anchors = sorted(_HTT_ANCHORS.items())
    lo_htt, (lo_T, _) = anchors[0]
    hi_htt, (hi_T, _) = anchors[-1]

    if htt_c <= lo_htt:
        peak_T = lo_T
    elif htt_c >= hi_htt:
        peak_T = hi_T
    else:
        # piecewise-linear interpolation of peak_T between the bracketing anchors
        peak_T = lo_T
        for (h0, (t0, _)), (h1, (t1, _)) in zip(anchors, anchors[1:]):
            if h0 <= htt_c <= h1:
                frac = (htt_c - h0) / (h1 - h0)
                peak_T = t0 + frac * (t1 - t0)
                break

    timestep = 1.0 if htt_c <= lo_htt else 0.5
    return AnnealSpec(peak_T_K=round(peak_T, 1), timestep_fs=timestep)


# --------------------------------------------------------------------------- #
# Fixed protocol durations (ns) — Wood experimental procedures
# --------------------------------------------------------------------------- #
_NVT_NS = 10.0
_ANNEAL_NS = 25.0
_FINAL_NS = 10.0
_FINAL_DT_FS = 2.0
# annealing schedule: hold peak (0–10 ns), cool to 300 K (10–20 ns), hold (20–25 ns)
_ANNEAL_TIMES_PS = (0, 10000, 20000, 25000)


def _nsteps(duration_ns: float, timestep_fs: float) -> int:
    return int(round(duration_ns * 1e6 / timestep_fs))


# --------------------------------------------------------------------------- #
# .mdp templates (parameterized by peak_T, timestep, nsteps)
# --------------------------------------------------------------------------- #
def _em_mdp() -> str:
    return """; Wood et al. 2024 condensation — energy minimization (steepest descent)
integrator      = steep
nsteps          = 50000
emtol           = 500.0
emstep          = 0.01

cutoff-scheme   = Verlet
nstlist         = 10
rlist           = 1.4

; Electrostatics — PME
coulombtype     = PME
rcoulomb        = 1.4
fourierspacing  = 0.12
pme-order       = 4

; Van der Waals — 1.4 nm
vdwtype         = cut-off
rvdw            = 1.4

constraints     = none
pbc             = xyz
"""


def _nvt_mdp(spec: AnnealSpec) -> str:
    dt = spec.timestep_fs / 1000.0
    nsteps = _nsteps(_NVT_NS, spec.timestep_fs)
    return f"""; Wood et al. 2024 condensation — NVT, {_NVT_NS:g} ns @ {spec.peak_T_K:g} K
integrator      = md
nsteps          = {nsteps}
dt              = {dt:.4f}

; Output
nstlog          = 5000
nstenergy       = 5000
nstxout-compressed = 50000
compressed-x-grps = System

; Neighbour searching
cutoff-scheme   = Verlet
nstlist         = 10
rlist           = 1.0

; Electrostatics — PME
coulombtype     = PME
rcoulomb        = 1.0
fourierspacing  = 0.12
pme-order       = 4

; Van der Waals — 1 nm
vdwtype         = cut-off
rvdw            = 1.0
DispCorr        = EnerPres

; Thermostat — velocity-rescale, 0.1 ps
tcoupl          = V-rescale
tc-grps         = System
tau_t           = 0.1
ref_t           = {spec.peak_T_K:g}

; NVT — no pressure coupling
pcoupl          = no

; Velocities generated at the peak temperature
gen_vel         = yes
gen_temp        = {spec.peak_T_K:g}
gen_seed        = -1

; Bonds — unconstrained (small timestep resolves bond vibrations at high T)
constraints     = none
pbc             = xyz
continuation    = no
"""


def _npt_anneal_mdp(spec: AnnealSpec) -> str:
    dt = spec.timestep_fs / 1000.0
    nsteps = _nsteps(_ANNEAL_NS, spec.timestep_fs)
    times = " ".join(str(t) for t in _ANNEAL_TIMES_PS)
    temps = f"{spec.peak_T_K:g} {spec.peak_T_K:g} 300 300"
    return f"""; Wood et al. 2024 condensation — NPT simulated annealing, {_ANNEAL_NS:g} ns
; hold {spec.peak_T_K:g} K (0–10 ns) -> cool to 300 K (10–20 ns) -> hold 300 K (20–25 ns)
integrator      = md
nsteps          = {nsteps}
dt              = {dt:.4f}

; Output
nstlog          = 5000
nstenergy       = 5000
nstxout-compressed = 50000
compressed-x-grps = System

; Neighbour searching
cutoff-scheme   = Verlet
nstlist         = 10
rlist           = 1.0

; Electrostatics — PME
coulombtype     = PME
rcoulomb        = 1.0
fourierspacing  = 0.12
pme-order       = 4

; Van der Waals — 1 nm
vdwtype         = cut-off
rvdw            = 1.0
DispCorr        = EnerPres

; Thermostat — velocity-rescale, 0.1 ps
tcoupl          = V-rescale
tc-grps         = System
tau_t           = 0.1
ref_t           = {spec.peak_T_K:g}

; Simulated annealing
annealing           = single
annealing-npoints   = 4
annealing-time      = {times}
annealing-temp      = {temps}

; Pressure — Berendsen isotropic, 100 bar, 1 ps
pcoupl          = Berendsen
pcoupltype      = isotropic
tau_p           = 1.0
ref_p           = 100.0
compressibility = 4.5e-5

gen_vel         = no
constraints     = none
pbc             = xyz
continuation    = yes
"""


def _npt_final_mdp() -> str:
    nsteps = _nsteps(_FINAL_NS, _FINAL_DT_FS)
    return f"""; Wood et al. 2024 condensation — final NPT equilibration, {_FINAL_NS:g} ns, 300 K, 1 bar
integrator      = md
nsteps          = {nsteps}
dt              = 0.002

; Output
nstlog          = 5000
nstenergy       = 5000
nstxout-compressed = 50000
compressed-x-grps = System

; Neighbour searching
cutoff-scheme   = Verlet
nstlist         = 10
rlist           = 1.4

; Electrostatics — PME
coulombtype     = PME
rcoulomb        = 1.4
fourierspacing  = 0.12
pme-order       = 4

; Van der Waals — 1.4 nm
vdwtype         = cut-off
rvdw            = 1.4
DispCorr        = EnerPres

; Thermostat — velocity-rescale, 0.1 ps
tcoupl          = V-rescale
tc-grps         = System
tau_t           = 0.1
ref_t           = 300

; Pressure — Berendsen isotropic, 1 bar, 1 ps
pcoupl          = Berendsen
pcoupltype      = isotropic
tau_p           = 1.0
ref_p           = 1.0
compressibility = 4.5e-5

gen_vel         = no

; 2 fs timestep — h-bond constraints
constraints     = h-bonds
constraint-algorithm = LINCS
lincs-iter      = 1
lincs-order     = 4
pbc             = xyz
continuation    = yes
"""


def render_mdp_set(spec: AnnealSpec) -> dict[str, str]:
    """The four dry `.mdp` files for one Wood condensation run."""
    return {
        "em.mdp": _em_mdp(),
        "nvt.mdp": _nvt_mdp(spec),
        "npt_anneal.mdp": _npt_anneal_mdp(spec),
        "npt_final.mdp": _npt_final_mdp(),
    }


# --------------------------------------------------------------------------- #
# Run-script renderer (chains the 4 stages, ×N repeats)
# --------------------------------------------------------------------------- #
def render_condensation_script(
    gro_name: str,
    top_name: str,
    spec: AnnealSpec,
    n_repeats: int = 3,
    gmx_bin: str = "gmx",
    ntomp: int = 8,
) -> str:
    """`run_condensation.sh` — chain EM → NVT → NPT anneal → NPT final for each
    of *n_repeats* independent repeats (Wood ran 3 with different start configs
    for ergodicity). Setup-only: this script is reviewed and run by the user.
    """
    if n_repeats < 1:
        raise CondensationError("n_repeats must be >= 1")

    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "# ============================================================",
        "# Wood et al. 2024 condensation annealing",
        f"# peak T = {spec.peak_T_K:g} K, timestep = {spec.timestep_fs:g} fs (final: 2 fs)",
        "# EM -> NVT 10 ns -> NPT anneal 25 ns (heat/cool) -> NPT final 10 ns",
        f"# {n_repeats} independent repeat(s) for ergodicity.",
        "# Generated by biochar.condensation -- review, then run on your GROMACS build.",
        "# ============================================================",
        "",
        f'GMX="{gmx_bin}"',
        'SIM="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        f"NTOMP={ntomp}",
        f'GRO="$SIM/{gro_name}"',
        f'TOP="$SIM/{top_name}"',
        "",
        f"for rep in $(seq 1 {n_repeats}); do",
        '  R="$SIM/rep_$rep"',
        '  mkdir -p "$R"',
        '  echo "=== repeat $rep ==="',
        "",
        "  # [1/4] Energy minimization",
        '  "$GMX" grompp -f "$SIM/em.mdp" -c "$GRO" -p "$TOP" -o "$R/em.tpr" -maxwarn 2',
        '  "$GMX" mdrun -v -ntmpi 1 -ntomp "$NTOMP" -s "$R/em.tpr" -deffnm "$R/em"',
        "",
        "  # [2/4] NVT (hot) — fresh velocities per repeat (gen_seed = -1)",
        '  "$GMX" grompp -f "$SIM/nvt.mdp" -c "$R/em.gro" -p "$TOP" -o "$R/nvt.tpr" -maxwarn 2',
        '  "$GMX" mdrun -v -ntmpi 1 -ntomp "$NTOMP" -s "$R/nvt.tpr" -deffnm "$R/nvt"',
        "",
        "  # [3/4] NPT simulated annealing (heat/cool)",
        '  "$GMX" grompp -f "$SIM/npt_anneal.mdp" -c "$R/nvt.gro" -t "$R/nvt.cpt" -p "$TOP" -o "$R/anneal.tpr" -maxwarn 2',
        '  "$GMX" mdrun -v -ntmpi 1 -ntomp "$NTOMP" -s "$R/anneal.tpr" -deffnm "$R/anneal"',
        "",
        "  # [4/4] NPT final equilibration (300 K, 1 bar)",
        '  "$GMX" grompp -f "$SIM/npt_final.mdp" -c "$R/anneal.gro" -t "$R/anneal.cpt" -p "$TOP" -o "$R/final.tpr" -maxwarn 2',
        '  "$GMX" mdrun -v -ntmpi 1 -ntomp "$NTOMP" -s "$R/final.tpr" -deffnm "$R/final"',
        '  echo "repeat $rep condensed structure: $R/final.gro"',
        "done",
        "",
        'echo "All repeats complete. Check RMS/density/box plateaus for convergence."',
        "",
    ]
    return "\n".join(lines) + "\n"


def write_condensation_setup(
    output_dir: str | Path,
    gro_name: str,
    top_name: str,
    htt_c: float | None = None,
    spec: AnnealSpec | None = None,
    n_repeats: int = 3,
    gmx_bin: str = "gmx",
    ntomp: int = 8,
) -> Path:
    """Write the four `.mdp` files + `run_condensation.sh` into *output_dir*.

    Provide either *htt_c* (mapped via :func:`anneal_spec_for_htt`) or an explicit
    *spec*. *gro_name*/*top_name* are the packed building-block system that a
    later phase produces and that this script consumes. Setup-only.
    """
    if spec is None:
        if htt_c is None:
            raise CondensationError("provide either htt_c or spec")
        spec = anneal_spec_for_htt(htt_c)

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, content in render_mdp_set(spec).items():
        (out / name).write_text(content)

    script = render_condensation_script(
        gro_name, top_name, spec, n_repeats=n_repeats, gmx_bin=gmx_bin, ntomp=ntomp
    )
    script_path = out / "run_condensation.sh"
    script_path.write_text(script)
    import stat as _stat
    script_path.chmod(script_path.stat().st_mode | _stat.S_IEXEC | _stat.S_IXGRP | _stat.S_IXOTH)
    return out
