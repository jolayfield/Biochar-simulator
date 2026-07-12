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

import math
import re
import shutil
import stat
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


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
# Packing: build a bulk box from N copies of a single biochar molecule
# --------------------------------------------------------------------------- #
_SECTION_RE = re.compile(r"^\s*\[\s*(.+?)\s*\]")


def moleculetype_name(itp_text: str) -> str:
    """Return the ``[ moleculetype ]`` name from a biochar molecule `.itp`."""
    lines = itp_text.splitlines()
    for i, ln in enumerate(lines):
        m = _SECTION_RE.match(ln)
        if m and m.group(1).lower() == "moleculetype":
            for data in lines[i + 1:]:
                s = data.strip()
                if s and not s.startswith(";"):
                    return s.split()[0]
    raise CondensationError("no [ moleculetype ] found in the molecule .itp")


def _gro_extent_nm(gro_text: str) -> float:
    """Largest xyz span (nm) of the molecule in a GROMACS `.gro` file."""
    lines = gro_text.splitlines()
    n = int(lines[1].strip())
    xs, ys, zs = [], [], []
    for ln in lines[2:2 + n]:
        # .gro columns are fixed-width: x,y,z are 8.3f at cols 20-28,28-36,36-44
        xs.append(float(ln[20:28])); ys.append(float(ln[28:36])); zs.append(float(ln[36:44]))
    return max(max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))


def estimate_box_nm(gro_text: str, n_copies: int, loose_factor: float = 1.6) -> float:
    """A generous cubic box side (nm) so *n_copies* pack loosely without overlap.

    ~one molecule per (extent x loose_factor)^3 cell, arranged in an n^(1/3) grid.
    Deliberately loose: Wood packs loosely, then the Berendsen 100-bar anneal
    contracts the box. Override with an explicit box if you prefer.
    """
    if n_copies < 1:
        raise CondensationError("n_copies must be >= 1")
    cell = max(_gro_extent_nm(gro_text), 0.5) * loose_factor
    side = cell * math.ceil(n_copies ** (1.0 / 3.0))
    return round(side, 3)


def render_condensation_top(itp_name: str, moltype: str, n_copies: int) -> str:
    """A dry (no water) `.top` for a bulk of *n_copies* of one molecule type."""
    if n_copies < 1:
        raise CondensationError("n_copies must be >= 1")
    return (
        '#include "oplsaa.ff/forcefield.itp"\n'
        f'#include "{itp_name}"\n\n'
        "[ system ]\n"
        f"Condensation bulk: {n_copies} x {moltype}\n\n"
        "[ molecules ]\n"
        f"{moltype:<20s} {n_copies}\n"
    )


@dataclass
class PackSpec:
    """How to pack the bulk box from one molecule's coordinates."""

    molecule_gro: str    # single-molecule .gro filename (in the run dir)
    n_copies: int
    box_nm: float        # cubic box side for gmx insert-molecules
    insertion_try: int = 200


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
    pack: Optional[PackSpec] = None,
) -> str:
    """`run_condensation.sh` — chain EM → NVT → NPT anneal → NPT final for each
    of *n_repeats* independent repeats (Wood ran 3 with different start configs
    for ergodicity). If *pack* is given, each repeat first packs N copies of a
    single molecule into a loose box with `gmx insert-molecules` (a distinct
    ``-seed`` per repeat gives the different starting configurations). Setup-only:
    this script is reviewed and run by the user.
    """
    if n_repeats < 1:
        raise CondensationError("n_repeats must be >= 1")

    header = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "# ============================================================",
        "# Wood et al. 2024 condensation annealing",
        f"# peak T = {spec.peak_T_K:g} K, timestep = {spec.timestep_fs:g} fs (final: 2 fs)",
    ]
    if pack is not None:
        header.append(
            f"# Pack {pack.n_copies} x molecule into a {pack.box_nm:g} nm box, then:"
        )
    header += [
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
    ]
    lines = header
    if pack is not None:
        b = pack.box_nm
        lines += [
            f"  # [0/4] Pack {pack.n_copies} copies into a loose {b:g} nm cubic box (seeded per repeat)",
            f'  "$GMX" insert-molecules -ci "$SIM/{pack.molecule_gro}" -nmol {pack.n_copies} '
            f'-box {b:g} {b:g} {b:g} -try {pack.insertion_try} -seed "$rep" -o "$GRO"',
            "",
        ]
    lines += [
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
    script_path.chmod(script_path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return out


def setup_condensation(
    output_dir: str | Path,
    molecule_gro: str | Path,
    molecule_itp: str | Path,
    n_copies: int,
    htt_c: float | None = None,
    spec: AnnealSpec | None = None,
    box_nm: float | None = None,
    n_repeats: int = 3,
    gmx_bin: str = "gmx",
    ntomp: int = 8,
) -> Path:
    """Set up a Wood condensation run that packs *n_copies* of ONE biochar
    molecule and anneals them into a bulk.

    Takes an existing single-molecule `.gro` + `.itp` (e.g. from
    ``BiocharGenerator.export_gromacs``). Writes into *output_dir*: the molecule
    files, a dry `.top` (molecule × n_copies), the four `.mdp`s, and
    `run_condensation.sh` which packs (``gmx insert-molecules``, seeded per
    repeat) then runs EM → NVT → NPT anneal → NPT final. Setup-only.

    Anneal temperature comes from *htt_c* (mapped via :func:`anneal_spec_for_htt`)
    or an explicit *spec*. The packing box defaults to a loose estimate
    (:func:`estimate_box_nm`); override with *box_nm*.
    """
    if spec is None:
        if htt_c is None:
            raise CondensationError("provide either htt_c or spec")
        spec = anneal_spec_for_htt(htt_c)
    if n_copies < 1:
        raise CondensationError("n_copies must be >= 1")

    molecule_gro, molecule_itp = Path(molecule_gro), Path(molecule_itp)
    if not molecule_gro.exists():
        raise CondensationError(f"molecule .gro not found: {molecule_gro}")
    if not molecule_itp.exists():
        raise CondensationError(f"molecule .itp not found: {molecule_itp}")

    gro_text = molecule_gro.read_text()
    itp_text = molecule_itp.read_text()
    moltype = moleculetype_name(itp_text)
    if box_nm is None:
        box_nm = estimate_box_nm(gro_text, n_copies)

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    shutil.copy(molecule_gro, out / molecule_gro.name)
    shutil.copy(molecule_itp, out / molecule_itp.name)

    (out / "system.top").write_text(
        render_condensation_top(molecule_itp.name, moltype, n_copies)
    )
    for name, content in render_mdp_set(spec).items():
        (out / name).write_text(content)

    pack = PackSpec(molecule_gro=molecule_gro.name, n_copies=n_copies, box_nm=box_nm)
    script = render_condensation_script(
        "packed.gro", "system.top", spec,
        n_repeats=n_repeats, gmx_bin=gmx_bin, ntomp=ntomp, pack=pack,
    )
    script_path = out / "run_condensation.sh"
    script_path.write_text(script)
    script_path.chmod(script_path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return out


def generate_and_condense(
    output_dir: str | Path,
    n_copies: int,
    generator_config=None,
    htt_c: float | None = None,
    spec: AnnealSpec | None = None,
    box_nm: float | None = None,
    n_repeats: int = 3,
    gmx_bin: str = "gmx",
    ntomp: int = 8,
) -> Path:
    """Generate one biochar molecule *fresh* and set up its condensation run.

    Thin wrapper over :func:`setup_condensation`: builds a molecule with
    :class:`biochar.biochar_generator.BiocharGenerator` (using *generator_config*,
    or a default :class:`GeneratorConfig`), exports its `.gro`/`.itp`, then packs
    *n_copies* and anneals. Requires RDKit (imported lazily so the rest of this
    module stays dependency-free).
    """
    from .biochar_generator import BiocharGenerator, GeneratorConfig

    cfg = generator_config or GeneratorConfig()
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    gen = BiocharGenerator(cfg)
    gen.generate()
    gro, _top, itp = gen.export_gromacs(str(out / "_molecule"), basename="molecule")

    return setup_condensation(
        out, gro, itp, n_copies, htt_c=htt_c, spec=spec, box_nm=box_nm,
        n_repeats=n_repeats, gmx_bin=gmx_bin, ntomp=ntomp,
    )


# --------------------------------------------------------------------------- #
# Phase 4: turn a condensed bulk into an exposed surface
# --------------------------------------------------------------------------- #
# Wood created two exposed surfaces per bulk by elongating the box in z so each
# periodic layer is separated by ~10 nm of vacuum, then re-equilibrating with
# semi-isotropic pressure coupling ("the surface-forming xy plane decoupled from
# the interlayer z direction"). Interpretation note: to *preserve* the vacuum
# gap, z must not be pressure-scaled toward 1 bar (that would collapse it), so
# the xy plane is coupled at 1 bar (compressibility 4.5e-5) while z is frozen
# (compressibility 0). The models are dry — solvation for adsorption is separate.
_DEFAULT_GAP_NM = 10.0


def _surface_npt_mdp() -> str:
    nsteps = _nsteps(_FINAL_NS, _FINAL_DT_FS)  # 10 ns @ 2 fs
    return f"""; Wood et al. 2024 — surface equilibration NPT, {_FINAL_NS:g} ns, 300 K, semi-isotropic
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

; Pressure — semi-isotropic: xy coupled at 1 bar, z frozen to preserve the
; vacuum gap (compressibility order is "xy z")
pcoupl          = Berendsen
pcoupltype      = semiisotropic
tau_p           = 1.0
ref_p           = 1.0 1.0
compressibility = 4.5e-5 0

gen_vel         = no

; 2 fs timestep — h-bond constraints
constraints     = h-bonds
constraint-algorithm = LINCS
lincs-iter      = 1
lincs-order     = 4
pbc             = xyz
continuation    = yes
"""


def render_surface_script(
    bulk_gro: str,
    top_name: str,
    gap_nm: float = _DEFAULT_GAP_NM,
    gmx_bin: str = "gmx",
    ntomp: int = 8,
) -> str:
    """`run_surface.sh` — z-expand a condensed bulk by *gap_nm* of vacuum, then
    EM + semi-isotropic NPT to relax the two exposed surfaces. Setup-only.
    """
    if gap_nm <= 0:
        raise CondensationError("gap_nm must be > 0")
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "# ============================================================",
        "# Wood et al. 2024 — create exposed surfaces from a condensed bulk",
        f"# z-expand by ~{gap_nm:g} nm vacuum, EM, then 10 ns semi-isotropic NPT",
        "# (xy at 1 bar, z frozen to keep the gap). Generated by biochar.condensation.",
        "# ============================================================",
        "",
        f'GMX="{gmx_bin}"',
        'SIM="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        f"NTOMP={ntomp}",
        f'BULK="$SIM/{bulk_gro}"',
        f'TOP="$SIM/{top_name}"',
        f"GAP={gap_nm:g}",
        "",
        "# [1/3] Elongate the box in z (read the current box, add the vacuum gap, re-center)",
        'read BX BY BZ _ <<< "$(tail -n 1 "$BULK")"',
        'NEWZ="$(awk "BEGIN{printf \\"%.5f\\", $BZ + $GAP}")"',
        'echo "box z: $BZ -> $NEWZ nm"',
        '"$GMX" editconf -f "$BULK" -o "$SIM/expanded.gro" -box "$BX" "$BY" "$NEWZ" -c',
        "",
        "# [2/3] Energy minimization of the expanded system",
        '"$GMX" grompp -f "$SIM/em.mdp" -c "$SIM/expanded.gro" -p "$TOP" -o "$SIM/surf_em.tpr" -maxwarn 2',
        '"$GMX" mdrun -v -ntmpi 1 -ntomp "$NTOMP" -s "$SIM/surf_em.tpr" -deffnm "$SIM/surf_em"',
        "",
        "# [3/3] Surface equilibration — 10 ns semi-isotropic NPT",
        '"$GMX" grompp -f "$SIM/surf_npt.mdp" -c "$SIM/surf_em.gro" -p "$TOP" -o "$SIM/surf_npt.tpr" -maxwarn 2',
        '"$GMX" mdrun -v -ntmpi 1 -ntomp "$NTOMP" -s "$SIM/surf_npt.tpr" -deffnm "$SIM/surf_npt"',
        "",
        'echo "Exposed-surface structure: $SIM/surf_npt.gro"',
        'echo "Analyse SASA over the last 2 ns for surface roughness (see validation helpers)."',
        "",
    ]
    return "\n".join(lines) + "\n"


def setup_surface(
    output_dir: str | Path,
    bulk_gro: str | Path,
    top: str | Path,
    itp: str | Path | None = None,
    gap_nm: float = _DEFAULT_GAP_NM,
    gmx_bin: str = "gmx",
    ntomp: int = 8,
) -> Path:
    """Set up the surface-creation run from a condensed bulk.

    *bulk_gro* is a converged condensed structure (e.g. a repeat's ``final.gro``);
    *top* is its topology (the ``system.top`` from :func:`setup_condensation`);
    *itp* is the molecule include the topology references (copied alongside so
    ``grompp`` resolves it). Writes ``em.mdp`` + ``surf_npt.mdp`` +
    ``run_surface.sh`` into *output_dir*. Setup-only.
    """
    bulk_gro, top = Path(bulk_gro), Path(top)
    for p in (bulk_gro, top):
        if not p.exists():
            raise CondensationError(f"input not found: {p}")

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    shutil.copy(bulk_gro, out / bulk_gro.name)
    shutil.copy(top, out / top.name)
    if itp is not None:
        itp = Path(itp)
        if not itp.exists():
            raise CondensationError(f"itp not found: {itp}")
        shutil.copy(itp, out / itp.name)

    (out / "em.mdp").write_text(_em_mdp())
    (out / "surf_npt.mdp").write_text(_surface_npt_mdp())

    script = render_surface_script(
        bulk_gro.name, top.name, gap_nm=gap_nm, gmx_bin=gmx_bin, ntomp=ntomp
    )
    script_path = out / "run_surface.sh"
    script_path.write_text(script)
    script_path.chmod(script_path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return out
