"""Generate GROMACS run inputs (.mdp templates + a driver script) for every
structure produced by a `biochar.sweep` run.

This is the stage *after* `biochar-sweep`: it turns each row of a sweep
manifest — a built `.gro`/`.top`/`.itp` triple — into a ready-to-submit
GROMACS directory following the dry-anneal / solvate / wet-equilibrate
protocol used by hand in this repo's `sim_200_v2` campaign (Wood et al. 2024
annealing schedule). GROMACS itself is not invoked here (no `gmx` binary in
this environment); this module only writes files, so it runs anywhere.

Pipeline per structure:
    dry:  EM -> NVT (anneal pre-eq) -> NPT (anneal, heat/cool) -> NPT (final, 1 bar)
    wet:  solvate + add ions -> EM -> NVT -> NPT (semiisotropic)

Usage
-----
    from biochar.md_setup import setup_md_from_manifest

    setup_md_from_manifest(
        "sweep_out/temperature_grid/manifest.csv",
        output_root="sweep_out/temperature_grid/md_runs",
        ion_profile="mn_calcareous_default",
    )

Or via CLI: ``biochar-md-setup sweep_out/.../manifest.csv --output-root ...``
"""

from __future__ import annotations

import csv
import dataclasses
import os
import shutil
import stat
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


class MDSetupError(Exception):
    """Raised for malformed manifests or unresolvable structure files."""


# ---------------------------------------------------------------------------
# Ion profiles — Minnesota water-chemistry scenarios for the wet stage.
#
# These are illustrative starting points (calcium-bicarbonate-type glacial
# till/drift groundwater is the dominant regional pattern in Minnesota), NOT
# a site-specific measurement. Concentrations vary widely by aquifer and
# county. Override with monitoring-well data for a specific site — see the
# MN DNR (Bulletin 26, "Natural Quality of Minnesota's Ground Water") and
# MPCA ambient groundwater monitoring reports for real values.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class IonProfile:
    """Bulk ion concentrations (mM) for `gmx genion`-style solvation.

    `counter_ion` / `counter_conc_mM` is whatever anion balances charge
    (typically Cl-); GROMACS topology neutrality is still enforced by
    `genion -neutral` at run time regardless of these targets.
    """

    name: str
    ca_mM: float = 0.0
    mg_mM: float = 0.0
    na_mM: float = 0.0
    k_mM: float = 0.0
    counter_ion: str = "CL"
    description: str = ""


ION_PROFILES: dict[str, IonProfile] = {
    "pure_water": IonProfile(
        name="pure_water",
        description="No background electrolyte — solvation-only control.",
    ),
    "mn_calcareous_default": IonProfile(
        name="mn_calcareous_default",
        ca_mM=1.5,
        mg_mM=0.6,
        na_mM=0.5,
        k_mM=0.1,
        description=(
            "Illustrative Ca-HCO3-type Minnesota glacial-till groundwater "
            "(hardness-dominated by Ca2+/Mg2+, low Na+/K+). NOT a measured "
            "value for any specific site -- replace with monitoring-well "
            "data (MN DNR / MPCA) for production runs."
        ),
    ),
    "na_dominated": IonProfile(
        name="na_dominated",
        na_mM=5.0,
        k_mM=0.2,
        description="Sodium-dominated scenario, e.g. softened or road-salt-impacted water.",
    ),
}


def get_ion_profile(name_or_profile) -> IonProfile:
    if isinstance(name_or_profile, IonProfile):
        return name_or_profile
    try:
        return ION_PROFILES[name_or_profile]
    except KeyError as exc:
        raise MDSetupError(
            f"Unknown ion_profile {name_or_profile!r}; choices: {sorted(ION_PROFILES)}"
        ) from exc


# ---------------------------------------------------------------------------
# .mdp templates — transcribed from sim_200_v2/*.mdp (Wood et al. 2024
# annealing protocol: dry EM -> NVT pre-eq -> NPT anneal (300->1000->300 K) ->
# NPT final; then wet EM -> NVT -> NPT semiisotropic).
# ---------------------------------------------------------------------------

DRY_EM_MDP = """; Dry Energy Minimization — steepest descent (Wood et al. 2024)
integrator      = steep
nsteps          = 50000
emtol           = 500.0
emstep          = 0.01

; Output
nstlog          = 500
nstenergy       = 500
nstxout-compressed = 0

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 10
rlist           = 1.4

; Electrostatics — PME
coulombtype     = PME
rcoulomb        = 1.4
fourierspacing  = 0.12
pme-order       = 4

; Van der Waals
vdwtype         = cut-off
rvdw            = 1.4
DispCorr        = no

; Bonds — no constraints for EM
constraints     = none

; Periodic boundary conditions
pbc             = xyz
"""

ANNEAL_NVT_MDP = """; NVT pre-equilibration — 0.5 ns at 300 K to safely initialize velocities
integrator      = md
nsteps          = 500000
dt              = 0.001

; Output
nstlog          = 10000
nstenergy       = 10000
nstxout         = 0
nstvout         = 0
nstfout         = 0
nstxout-compressed = 100000
compressed-x-grps = System

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 10
rlist           = 1.0

; Electrostatics — PME
coulombtype     = PME
rcoulomb        = 1.0
fourierspacing  = 0.12
pme-order       = 4

; Van der Waals
vdwtype         = cut-off
rvdw            = 1.0
DispCorr        = no

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = System
tau_t           = 0.1
ref_t           = 300

; Velocity generation
gen_vel         = yes
gen_temp        = 300
gen_seed        = -1

; Bonds
constraints     = h-bonds
constraint-algorithm = LINCS
lincs-iter      = 1
lincs-order     = 4

; Periodic boundary conditions
pbc             = xyz

; Continuation
continuation    = no
"""

ANNEAL_NPT_MDP = """; Simulated annealing — NPT, 5.5 ns, Berendsen 100 bar (Wood et al. 2024, scaled)
; Schedule: hold 300K (0.5ns) -> heat to 1000K (1ns) -> hold 1000K (1ns) -> cool to 300K (2ns) -> hold 300K (1ns)
integrator      = md
nsteps          = 5500000
dt              = 0.001

; Output
nstlog          = 10000
nstenergy       = 10000
nstxout         = 0
nstvout         = 0
nstfout         = 0
nstxout-compressed = 100000
compressed-x-grps = System

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 10
rlist           = 1.0

; Electrostatics — PME
coulombtype     = PME
rcoulomb        = 1.0
fourierspacing  = 0.12
pme-order       = 4

; Van der Waals
vdwtype         = cut-off
rvdw            = 1.0
DispCorr        = no

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = System
tau_t           = 0.1
ref_t           = 300

; Simulated annealing schedule — gradual heat-up then cool-down
annealing           = single
annealing-npoints   = 6
annealing-time      = 0 500 1500 2500 4500 5500
annealing-temp      = 300 300 1000 1000 300 300

; Pressure coupling — Berendsen, 100 bar, isotropic
pcoupl          = Berendsen
pcoupltype      = isotropic
tau_p           = 1.0
ref_p           = 100.0
compressibility = 4.5e-5
refcoord-scaling = com

; Velocity generation
gen_vel         = no

; Bonds — no constraints at high T; 1 fs timestep resolves bond vibrations
constraints     = none

; Periodic boundary conditions
pbc             = xyz

; Continuation
continuation    = yes
"""

FINAL_NPT_MDP = """; Final dry equilibration — NPT, 2 ns, Berendsen 1 bar, 300 K (Wood et al. 2024, scaled)
integrator      = md
nsteps          = 1000000
dt              = 0.002

; Output
nstlog          = 5000
nstenergy       = 5000
nstxout         = 0
nstvout         = 0
nstfout         = 0
nstxout-compressed = 50000
compressed-x-grps = System

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 10
rlist           = 1.4

; Electrostatics — PME
coulombtype     = PME
rcoulomb        = 1.4
fourierspacing  = 0.12
pme-order       = 4

; Van der Waals
vdwtype         = cut-off
rvdw            = 1.4
rvdw-switch     = 1.2
DispCorr        = EnerPres

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = System
tau_t           = 0.1
ref_t           = 300

; Pressure coupling — Berendsen, 1 bar, isotropic
pcoupl          = Berendsen
pcoupltype      = isotropic
tau_p           = 1.0
ref_p           = 1.0
compressibility = 4.5e-5
refcoord-scaling = com

; Velocity generation
gen_vel         = no

; Bonds
constraints     = h-bonds
constraint-algorithm = LINCS
lincs-iter      = 1
lincs-order     = 4

; Periodic boundary conditions
pbc             = xyz

; Continuation
continuation    = yes
"""

WET_EM_MDP = """; Wet Energy Minimization — steepest descent
integrator      = steep
nsteps          = 50000
emtol           = 1000.0
emstep          = 0.01

; Output
nstlog          = 500
nstenergy       = 500
nstxout-compressed = 0

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 10
rlist           = 1.2

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.2
fourierspacing  = 0.12
pme-order       = 4

; Van der Waals
vdwtype         = cut-off
rvdw            = 1.2
rvdw-switch     = 1.0
DispCorr        = EnerPres

; Bonds — no constraints for EM
constraints     = none

; Periodic boundary conditions
pbc             = xyz
"""

WET_NVT_MDP = """; Wet NVT Equilibration — 100 ps
integrator      = md
nsteps          = 50000
dt              = 0.002

; Output
nstlog          = 500
nstenergy       = 500
nstxout         = 0
nstvout         = 0
nstfout         = 0
nstxout-compressed = 5000
compressed-x-grps = System

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 10
rlist           = 1.2

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.2
fourierspacing  = 0.12
pme-order       = 4

; Van der Waals
vdwtype         = cut-off
rvdw            = 1.2
rvdw-switch     = 1.0
DispCorr        = EnerPres

; Temperature coupling — separate biochar and water
tcoupl          = V-rescale
tc-grps         = non-Water Water
tau_t           = 0.1 0.1
ref_t           = 300 300

; Pressure coupling
pcoupl          = no

; Velocity generation
gen_vel         = yes
gen_temp        = 300
gen_seed        = -1

; Bonds
constraints     = h-bonds
constraint-algorithm = LINCS
lincs-iter      = 1
lincs-order     = 4

; Periodic boundary conditions
pbc             = xyz

; Continuation
continuation    = no
"""

WET_NPT_MDP = """; Wet NPT Equilibration — 100 ps, semiisotropic C-rescale (XY fixed, Z free)
integrator      = md
nsteps          = 50000
dt              = 0.002

; Output
nstlog          = 500
nstenergy       = 500
nstxout         = 0
nstvout         = 0
nstfout         = 0
nstxout-compressed = 5000
compressed-x-grps = System

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 10
rlist           = 1.2

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.2
fourierspacing  = 0.12
pme-order       = 4

; Van der Waals
vdwtype         = cut-off
rvdw            = 1.2
rvdw-switch     = 1.0
DispCorr        = EnerPres

; Temperature coupling — separate biochar and water
tcoupl          = V-rescale
tc-grps         = non-Water Water
tau_t           = 0.1 0.1
ref_t           = 300 300

; Pressure coupling — C-rescale semiisotropic (XY frozen, Z free)
pcoupl          = C-rescale
pcoupltype      = semiisotropic
tau_p           = 2.0
ref_p           = 1.0 1.0
compressibility = 0 4.5e-5
refcoord-scaling = com

; Velocity generation
gen_vel         = no

; Bonds
constraints     = h-bonds
constraint-algorithm = LINCS
lincs-iter      = 1
lincs-order     = 4

; Periodic boundary conditions
pbc             = xyz

; Continuation
continuation    = yes
"""

_MDP_FILES = {
    "dry_em.mdp": DRY_EM_MDP,
    "anneal_nvt.mdp": ANNEAL_NVT_MDP,
    "anneal_npt.mdp": ANNEAL_NPT_MDP,
    "final_npt.mdp": FINAL_NPT_MDP,
    "wet_em.mdp": WET_EM_MDP,
    "wet_nvt.mdp": WET_NVT_MDP,
    "wet_npt.mdp": WET_NPT_MDP,
}


@dataclass
class MoleculeInsertion:
    """One `gmx insert-molecules` call: put *n_copies* of *gro_file* into the box."""

    gro_file: str            # a .gro filename staged into the run dir (see PreSolvationStage.extra_files)
    n_copies: int
    n_try: int = 500         # gmx insert-molecules -try


@dataclass
class PreSolvationStage:
    """
    A generic extra stage spliced into the pipeline **after dry-anneal, before
    box padding/solvation**: insert some molecules into the equilibrated slab's
    own box, then switch to a topology that already accounts for them.

    This is a domain-neutral extension point -- ``md_setup`` knows nothing about
    what the inserted molecules are. A downstream workflow (e.g. ligand/PFAS
    sorption) builds the merged topology + coordinate files itself and hands
    them in here; ``md_setup`` renders the ``insert-molecules`` chain, copies
    ``extra_files`` into the run dir, and uses ``solvation_top`` from solvation
    onward.
    """

    name: str                                   # banner label in the generated script
    insertions: list[MoleculeInsertion]         # rendered into the insert-molecules chain
    solvation_top: str                          # topology filename to use from solvation onward
    extra_files: list[str] = field(default_factory=list)  # source paths copied into the run dir


@dataclass
class MDSetupConfig:
    """Options for turning one manifest row into a run directory."""

    solvent_pad_nm: float = 1.2       # min solute-to-box-edge distance for editconf
    water_model: str = "spce"        # spce | tip3p | tip4p (must match force field .itp)
    ion_profile: str = "mn_calcareous_default"
    ntomp: int = 8
    gmx_bin: str = "gmx"

    # --- cluster (SLURM) submission, optional ---
    # When `cluster="slurm"`, `setup_one_structure` ALSO writes a
    # `submit_chain.sh` launcher that submits the 7 MD stages as dependent
    # sbatch jobs through a site GROMACS wrapper script, instead of (or in
    # addition to) the local `run_pipeline.sh` (see `_render_slurm_chain_script`).
    # The three site-specific fields below have no universal default -- set them
    # to your cluster's values; `setup_one_structure` raises if they are empty
    # when `cluster="slurm"`.
    cluster: Optional[str] = None
    slurm_submit_script: str = ""   # path to the site's grompp/mdrun sbatch wrapper
    slurm_sif_file: str = ""        # path to the GROMACS Singularity/Apptainer image
    slurm_partition: str = ""       # SLURM partition (add slurm_gpu=True for a GPU partition)
    slurm_cpus_per_task: int = 8
    slurm_gpu: bool = False
    slurm_preprocess_time: str = "00:30:00"   # wall time for the editconf/solvate/genion job

    # --- Optional extra molecule-insertion stage before solvation ---
    # When set, an `insert-molecules` stage is spliced into the dry structure
    # right after the final anneal and before box padding/solvation, and every
    # stage from that point on uses `pre_solvation_stage.solvation_top`. This is
    # domain-neutral: a downstream workflow supplies the coordinates + merged
    # topology (see `PreSolvationStage`). `md_setup` renders and stages only.
    pre_solvation_stage: Optional[PreSolvationStage] = None


def _read_gro_dims(gro_path: Path) -> tuple[int, tuple[float, float, float]]:
    """Return (n_atoms, box_vector) from a .gro file without invoking gmx."""
    lines = gro_path.read_text().splitlines()
    if len(lines) < 3:
        raise MDSetupError(f"{gro_path} does not look like a valid .gro file")
    n_atoms = int(lines[1].strip())
    box_line = lines[2 + n_atoms].split()
    box = tuple(float(x) for x in box_line[:3])
    return n_atoms, box


def _write_mdp_templates(run_dir: Path) -> None:
    for filename, content in _MDP_FILES.items():
        (run_dir / filename).write_text(content)


def _render_insertion_lines(
    stage: PreSolvationStage,
    input_gro_expr: str,
    output_gro_expr: str,
    gmx_expr: str = '"$GMX"',
    sim_prefix: str = "$SIM/",
    indent: str = "",
) -> list[str]:
    """Chain one `gmx insert-molecules -ci <file>.gro -nmol n -try k` call per
    entry in `stage.insertions`, inside the existing (pre-solvation) box. Each
    call's `-o` feeds the next call's `-f`, so N insertions compose into one
    decorated structure; the first call reads `input_gro_expr` and the last
    call writes `output_gro_expr`. `sim_prefix` is prepended to the staged
    `.gro` filenames ("$SIM/" for the local script; "" when the script's cwd is
    already the run dir, e.g. the SLURM scratch copy).
    """
    lines = []
    prev_expr = input_gro_expr
    n = len(stage.insertions)
    for i, ins in enumerate(stage.insertions):
        is_last = i == n - 1
        stem = Path(ins.gro_file).stem
        out_expr = output_gro_expr if is_last else f'"{sim_prefix}_pre_{stem}.gro"'
        lines.append(
            f'{indent}{gmx_expr} insert-molecules -f {prev_expr} '
            f'-ci "{sim_prefix}{ins.gro_file}" -nmol {ins.n_copies} -try {ins.n_try} '
            f'-o {out_expr}'
        )
        prev_expr = out_expr
    return lines


def _render_pipeline_script(
    run_dir: Path,
    label: str,
    gro_name: str,
    top_name: str,
    cfg: MDSetupConfig,
    ion: IonProfile,
    pre_solvation_stage: Optional[PreSolvationStage] = None,
) -> str:
    """Build a run_pipeline.sh for one structure: dry anneal, then solvate + wet stages.

    Box sizing and ion counts are computed by `gmx` at submit time from the
    *actual* equilibrated structure (via `editconf`/`solvate`/`genion -conc`)
    rather than hand-tuned per run — the failure mode visible in this repo's
    `sim_200_v2/#slab_for_solvation.gro.N#` backups, where the pad/box had to
    be manually re-guessed after each dry run finished at a different size.

    If `pre_solvation_stage` is given, an extra `insert-molecules` stage is
    spliced in right after the dry annealing and before box padding/solvation:
    the surface is annealed bare first, then the requested molecules are
    inserted into the already-equilibrated slab's own box, and every stage from
    that point on (box+solvate, genion, wet EM/NVT/NPT) uses
    `pre_solvation_stage.solvation_top` instead of the bare biochar topology.
    """
    gmx = cfg.gmx_bin
    # genion only supports one -pname/-nname pair per invocation in most
    # GROMACS builds; multi-cation profiles need sequential genion calls.
    # We emit one genion stage per cation actually requested (>0 mM), each
    # followed by -neutral so the final topology is always charge-neutral.
    # Valence: Ca2+/Mg2+ are divalent (-pq 2); Na+/K+ are monovalent (-pq 1).
    cations = [
        (sym, conc, pq)
        for sym, conc, pq in (
            ("CA", ion.ca_mM, 2),
            ("MG", ion.mg_mM, 2),
            ("NA", ion.na_mM, 1),
            ("K", ion.k_mM, 1),
        )
        if conc > 0
    ]

    script_lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "# ============================================================",
        f"# Biochar MD pipeline — {label}",
        "# Generated by biochar.md_setup — do not hand-edit box/ion values;",
        "# re-run biochar-md-setup if the input structure changes.",
        "# Dry: EM -> NVT(anneal pre-eq) -> NPT(anneal 300->1000->300K) -> NPT(final, 1 bar)",
        "# Wet: box+solvate+ions -> EM -> NVT -> NPT (semiisotropic)",
        "# ============================================================",
        "",
        f'GMX="{gmx}"',
        'SIM="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        f'NTOMP={cfg.ntomp}',
        "",
        'LOGFILE="$SIM/pipeline.log"',
        'exec > >(tee -a "$LOGFILE") 2>&1',
        f'echo "Pipeline started for {label}: $(date)"',
        "",
        "# ---- [1/8] Dry Energy Minimization ----",
        'mkdir -p "$SIM/dry_em"',
        f'"$GMX" grompp -f "$SIM/dry_em.mdp" -c "$SIM/{gro_name}" -p "$SIM/{top_name}" -o "$SIM/dry_em/em.tpr" -maxwarn 2',
        '"$GMX" mdrun -v -ntmpi 1 -ntomp "$NTOMP" -s "$SIM/dry_em/em.tpr" -deffnm "$SIM/dry_em/em"',
        "",
        "# ---- [2/8] Anneal NVT pre-equilibration ----",
        'mkdir -p "$SIM/anneal_nvt"',
        f'"$GMX" grompp -f "$SIM/anneal_nvt.mdp" -c "$SIM/dry_em/em.gro" -p "$SIM/{top_name}" -o "$SIM/anneal_nvt/nvt.tpr" -maxwarn 2',
        '"$GMX" mdrun -v -ntmpi 1 -ntomp "$NTOMP" -s "$SIM/anneal_nvt/nvt.tpr" -deffnm "$SIM/anneal_nvt/nvt"',
        "",
        "# ---- [3/8] Simulated annealing NPT (300 -> 1000 -> 300 K) ----",
        'mkdir -p "$SIM/anneal_npt"',
        f'"$GMX" grompp -f "$SIM/anneal_npt.mdp" -c "$SIM/anneal_nvt/nvt.gro" -p "$SIM/{top_name}" -o "$SIM/anneal_npt/npt.tpr" -maxwarn 2',
        '"$GMX" mdrun -v -ntmpi 1 -ntomp "$NTOMP" -s "$SIM/anneal_npt/npt.tpr" -deffnm "$SIM/anneal_npt/npt"',
        "",
        "# ---- [4/8] Final dry NPT (1 bar, 300 K) ----",
        'mkdir -p "$SIM/final_npt"',
        f'"$GMX" grompp -f "$SIM/final_npt.mdp" -c "$SIM/anneal_npt/npt.gro" -p "$SIM/{top_name}" -o "$SIM/final_npt/npt.tpr" -maxwarn 2',
        '"$GMX" mdrun -v -ntmpi 1 -ntomp "$NTOMP" -s "$SIM/final_npt/npt.tpr" -deffnm "$SIM/final_npt/npt"',
        'echo "Dry annealing complete: $(date)"',
        "",
    ]

    # Solvation-stage topology: switch to the stage's topology if an insertion
    # stage was provided (it already accounts for the inserted molecules);
    # otherwise keep using the bare biochar one.
    solv_top_name = pre_solvation_stage.solvation_top if pre_solvation_stage else top_name
    dry_structure_expr = '"$SIM/final_npt/npt.gro"'

    if pre_solvation_stage:
        insertion_lines = _render_insertion_lines(
            pre_solvation_stage,
            input_gro_expr=dry_structure_expr,
            output_gro_expr='"$SIM/pre_solvation.gro"',
        )
        script_lines += [
            f"# ---- [4b/8] {pre_solvation_stage.name} ----",
            f"# Inserts: {', '.join(i.gro_file for i in pre_solvation_stage.insertions)} "
            f"(solvation topology: {pre_solvation_stage.solvation_top})",
        ]
        script_lines += insertion_lines
        script_lines += [""]
        dry_structure_expr = '"$SIM/pre_solvation.gro"'

    script_lines += [
        "# ---- [5/8] Box + solvate ----",
        "# Box padding computed from the final dry structure, not hand-tuned:",
        f'"$GMX" editconf -f {dry_structure_expr} -o "$SIM/boxed.gro" -c -d {cfg.solvent_pad_nm} -bt cubic',
        f'"$GMX" solvate -cp "$SIM/boxed.gro" -cs {cfg.water_model}.gro -p "$SIM/wet.top" -o "$SIM/solvated.gro"',
        f'cp "$SIM/{solv_top_name}" "$SIM/wet.top.base"',
        f'grep -q "{cfg.water_model}" "$SIM/wet.top" || echo \'#include "oplsaa.ff/{cfg.water_model}.itp"\' >> "$SIM/wet.top"',
        "",
        "# ---- [6/8] Add ions to reach the target water chemistry ----",
        f'# Ion profile: {ion.name} ({ion.description})',
        f'"$GMX" grompp -f "$SIM/wet_em.mdp" -c "$SIM/solvated.gro" -p "$SIM/wet.top" -o "$SIM/genion.tpr" -maxwarn 2',
    ]
    if cations:
        for sym, conc_mM, pq in cations:
            script_lines.append(
                f'echo SOL | "$GMX" genion -s "$SIM/genion.tpr" -p "$SIM/wet.top" '
                f'-o "$SIM/solvated.gro" -pname {sym} -pq {pq} -nname {ion.counter_ion} '
                f'-conc {conc_mM / 1000.0:.6f} -neutral'
            )
    else:
        script_lines.append(
            f'echo SOL | "$GMX" genion -s "$SIM/genion.tpr" -p "$SIM/wet.top" '
            f'-o "$SIM/solvated.gro" -pname NA -nname {ion.counter_ion} -neutral'
        )
    script_lines += [
        "",
        "# ---- [7/8] Wet EM ----",
        'mkdir -p "$SIM/wet_em"',
        f'"$GMX" grompp -f "$SIM/wet_em.mdp" -c "$SIM/solvated.gro" -p "$SIM/wet.top" -o "$SIM/wet_em/em.tpr" -maxwarn 2',
        '"$GMX" mdrun -v -ntmpi 1 -ntomp "$NTOMP" -s "$SIM/wet_em/em.tpr" -deffnm "$SIM/wet_em/em"',
        "",
        "# ---- [8/8] Wet NVT then NPT (semiisotropic) ----",
        'mkdir -p "$SIM/wet_nvt" "$SIM/wet_npt"',
        f'"$GMX" grompp -f "$SIM/wet_nvt.mdp" -c "$SIM/wet_em/em.gro" -p "$SIM/wet.top" -o "$SIM/wet_nvt/nvt.tpr" -maxwarn 2',
        '"$GMX" mdrun -v -ntmpi 1 -ntomp "$NTOMP" -s "$SIM/wet_nvt/nvt.tpr" -deffnm "$SIM/wet_nvt/nvt"',
        f'"$GMX" grompp -f "$SIM/wet_npt.mdp" -c "$SIM/wet_nvt/nvt.gro" -p "$SIM/wet.top" -o "$SIM/wet_npt/npt.tpr" -maxwarn 2',
        '"$GMX" mdrun -v -ntmpi 1 -ntomp "$NTOMP" -s "$SIM/wet_npt/npt.tpr" -deffnm "$SIM/wet_npt/npt"',
        "",
        'echo "Pipeline complete for {label}: $(date)"'.replace("{label}", label),
        'echo "Production-ready structure: $SIM/wet_npt/npt.gro (with $SIM/wet.top)"',
        "",
    ]
    return "\n".join(script_lines) + "\n"


def _render_solvate_ions_slurm(
    cfg: MDSetupConfig, ion: IonProfile,
    pre_solvation_stage: Optional[PreSolvationStage] = None,
) -> str:
    """A small custom sbatch script for the (optional insertion +)
    box+solvate+genion step.

    `gromacs_plumed.slurm` (the site's GROMACS wrapper) only runs
    `grompp`+`mdrun`; it has no `editconf`/`solvate`/`genion` support. This
    script fills that one gap, following the same container/scratch/results
    conventions as the wrapper so it drops into the same dependency chain.
    Positional args: $1 = input .gro (final dry structure), $2 = topology
    (.top) to extend with water + ions -- pass the stage's `solvation_top`
    here when `pre_solvation_stage` is set, so the atom count already matches
    the decorated coordinates this script produces -- $3 = results folder name.

    When `pre_solvation_stage` is set, its molecule(s) are inserted into the dry
    structure's own box (`gmx insert-molecules`) before `editconf`/`solvate`,
    mirroring the local `run_pipeline.sh`'s insertion stage. The staged `.gro`
    files are expected alongside this script in the run directory (put there by
    `setup_one_structure` via `PreSolvationStage.extra_files`) -- the `rsync`
    below already pulls the whole run directory into scratch, so they arrive
    for free.

    All `gmx` calls run directly via `singularity exec` (no generated inner
    script / no sed rewriting) so this file's content is exactly what runs.
    """
    cations = [
        (sym, conc, pq)
        for sym, conc, pq in (
            ("CA", ion.ca_mM, 2),
            ("MG", ion.mg_mM, 2),
            ("NA", ion.na_mM, 1),
            ("K", ion.k_mM, 1),
        )
        if conc > 0
    ]
    if cations:
        genion_lines = [
            f'echo SOL | singularity exec "$sif_file" gmx genion -s genion.tpr -p wet.top '
            f'-o solvated.gro -pname {sym} -pq {pq} -nname {ion.counter_ion} '
            f'-conc {conc_mM / 1000.0:.6f} -neutral'
            for sym, conc_mM, pq in cations
        ]
    else:
        genion_lines = [
            f'echo SOL | singularity exec "$sif_file" gmx genion -s genion.tpr -p wet.top '
            f'-o solvated.gro -pname NA -nname {ion.counter_ion} -neutral'
        ]
    genion_block = "\n".join(genion_lines)

    if pre_solvation_stage:
        insertion_lines = _render_insertion_lines(
            pre_solvation_stage,
            input_gro_expr="_in.gro",
            output_gro_expr="pre_solvation.gro",
            gmx_expr='singularity exec "$sif_file" gmx',
            sim_prefix="",   # script cwd is already the scratch run dir
        )
        insertion_block = "\n".join(insertion_lines) + "\n"
        dry_structure_var = "pre_solvation.gro"
    else:
        insertion_block = ""
        dry_structure_var = "_in.gro"

    return f"""#!/bin/bash
#SBATCH --job-name=solvate_ions
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition={cfg.slurm_partition}
#SBATCH --time={cfg.slurm_preprocess_time}

# ============================================================================
# Box + solvate + genion -- the one gap in gromacs_plumed.slurm (grompp/mdrun
# only). Generated by biochar.md_setup; mirrors that script's scratch/results
# conventions so it drops into the same --dependency=afterok chain.
#
# Args: $1 = input .gro (final dry structure), $2 = input .top, $3 = results name
# Ion profile: {ion.name} ({ion.description})
# ============================================================================

set -euo pipefail
sif_file="{cfg.slurm_sif_file}"
scratch_path=/scratch/$USER.$SLURM_JOBID
storage_path=$SLURM_SUBMIT_DIR
in_gro=$1
in_top=$2
results_name=${{3:-solvate_ions}}
results_path=$storage_path/$results_name.$SLURM_JOBID

mkdir -p $scratch_path
rsync -a "$storage_path"/ "$scratch_path"/
cd $scratch_path

cp "$in_gro" ./_in.gro
cp "$in_top" ./wet.top
touch .marker

{insertion_block}singularity exec "$sif_file" gmx editconf -f {dry_structure_var} -o boxed.gro -c -d {cfg.solvent_pad_nm} -bt cubic
singularity exec "$sif_file" gmx solvate -cp boxed.gro -cs {cfg.water_model}.gro -p wet.top -o solvated.gro
grep -q "{cfg.water_model}" wet.top || echo '#include "oplsaa.ff/{cfg.water_model}.itp"' >> wet.top
singularity exec "$sif_file" gmx grompp -f wet_em.mdp -c solvated.gro -p wet.top -o genion.tpr -maxwarn 2
{genion_block}
run_status=$?

mkdir -p $results_path
find . -newer .marker -type f -print0 | rsync -a --files-from=- --from0 ./ "$results_path"/
cd $storage_path
rm -rf $scratch_path
exit $run_status
"""


def _render_slurm_chain_script(
    run_dir: Path,
    label: str,
    gro_name: str,
    top_name: str,
    cfg: MDSetupConfig,
    pre_solvation_stage: Optional[PreSolvationStage] = None,
) -> str:
    """Build `submit_chain.sh`: sbatch the 7 (or 8, with an insertion stage) MD
    stages as a SLURM dependency chain through the site's `gromacs_plumed.slurm`
    wrapper, with a small custom job (`solvate_ions.slurm`) filling the
    solvate/genion gap between the dry and wet halves -- and, when
    `pre_solvation_stage` is given, inserting its molecule(s) inside that same
    job before box/solvate (see `_render_solvate_ions_slurm`).

    The dry stages always run against the bare biochar topology (`top_name`)
    -- insertion only changes the *coordinates*, not the biochar's own bonded
    terms, so there is nothing gained by annealing extra atoms through 3 dry
    stages. `solvate_ions.slurm` is handed `pre_solvation_stage.solvation_top`
    instead of `top_name` precisely at the point molecules are inserted, so
    every wet-stage grompp call downstream sees an atom count that matches the
    decorated coordinates.

    Each `sbatch --parsable` call prints a jobid immediately, before the job
    runs -- so the *next* stage's `-results-name.<jobid>/` input path is
    known at submit time even though the job that will create it hasn't
    executed yet. `--dependency=afterok:<jobid>` is what makes that safe.
    """
    script = cfg.slurm_submit_script
    part = cfg.slurm_partition
    cpus = cfg.slurm_cpus_per_task
    gpu_sbatch = "" if cfg.slurm_gpu else " --gpus=0"
    cpu_flag = "" if cfg.slurm_gpu else " --cpu"
    wet_top_expr = '"$SIM/solvate_ions.$J_solvate_ions/wet.top"'

    def gmx_stage(mdp: str, conf_expr: str, top_expr: str, out: str, name: str,
                  dep_var: Optional[str]) -> list[str]:
        dep = f" --dependency=afterok:${dep_var}" if dep_var else ""
        return [
            f'J_{name}=$(sbatch --parsable -p {part}{gpu_sbatch} --cpus-per-task={cpus}{dep} \\',
            f'  "{script}"{cpu_flag} -f "$SIM/{mdp}" -c {conf_expr} -p {top_expr} \\',
            f'  -o {out} -maxwarn 2 -results-dir "$SIM" -results-name {name})',
            f'echo "[{name}] submitted as job $J_{name}"',
            "",
        ]

    dry_top_expr = f'"$SIM/{top_name}"'
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "# ============================================================",
        f"# SLURM submission chain -- {label}",
        "# Generated by biochar.md_setup. Submits the 7 MD stages as a",
        "# dependency chain through the site GROMACS wrapper (slurm_submit_script),",
        "# with one custom job (solvate_ions.slurm) filling the editconf/solvate/genion gap.",
        "# Run from this directory: ./submit_chain.sh",
        "# ============================================================",
        "",
        'SIM="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        'cd "$SIM"',
        "",
        "# ---- [1/8] Dry Energy Minimization ----",
    ]
    lines += gmx_stage("dry_em.mdp", f'"$SIM/{gro_name}"', dry_top_expr, "em", "dry_em", None)
    lines += ["# ---- [2/8] Anneal NVT pre-equilibration ----"]
    lines += gmx_stage("anneal_nvt.mdp", '"$SIM/dry_em.$J_dry_em/em.gro"', dry_top_expr,
                        "nvt", "anneal_nvt", "J_dry_em")
    lines += ["# ---- [3/8] Simulated annealing NPT (300 -> 1000 -> 300 K) ----"]
    lines += gmx_stage("anneal_npt.mdp", '"$SIM/anneal_nvt.$J_anneal_nvt/nvt.gro"', dry_top_expr,
                        "npt", "anneal_npt", "J_anneal_nvt")
    lines += ["# ---- [4/8] Final dry NPT (1 bar, 300 K) ----"]
    lines += gmx_stage("final_npt.mdp", '"$SIM/anneal_npt.$J_anneal_npt/npt.gro"', dry_top_expr,
                        "npt", "final_npt", "J_anneal_npt")

    solvate_ions_top_name = (
        pre_solvation_stage.solvation_top if pre_solvation_stage else top_name
    )
    lines += [
        "# ---- [5/8] Box + solvate + ions (custom job; wrapper has no editconf/solvate;",
        "#      also does molecule insertion here if this structure has a stage) ----",
        f'J_solvate_ions=$(sbatch --parsable -p {part} --dependency=afterok:$J_final_npt \\',
        f'  "$SIM/solvate_ions.slurm" "$SIM/final_npt.$J_final_npt/npt.gro" "$SIM/{solvate_ions_top_name}" solvate_ions)',
        'echo "[solvate_ions] submitted as job $J_solvate_ions"',
        "",
        "# ---- [6/8] Wet EM ----",
    ]
    lines += gmx_stage("wet_em.mdp", '"$SIM/solvate_ions.$J_solvate_ions/solvated.gro"', wet_top_expr,
                        "em", "wet_em", "J_solvate_ions")
    lines += ["# ---- [7/8] Wet NVT ----"]
    lines += gmx_stage("wet_nvt.mdp", '"$SIM/wet_em.$J_wet_em/em.gro"', wet_top_expr,
                        "nvt", "wet_nvt", "J_wet_em")
    lines += ["# ---- [8/8] Wet NPT (semiisotropic) ----"]
    lines += gmx_stage("wet_npt.mdp", '"$SIM/wet_nvt.$J_wet_nvt/nvt.gro"', wet_top_expr,
                        "npt", "wet_npt", "J_wet_nvt")

    lines += [
        f'echo "Chain submitted for {label}. Track with: squeue -u $USER"',
        'echo "Final production structure will land in wet_npt.$J_wet_npt/npt.gro"',
        "",
    ]
    return "\n".join(lines) + "\n"


def setup_one_structure(
    gro_path: str | Path,
    top_path: str | Path,
    output_dir: str | Path,
    label: str = "structure",
    config: Optional[MDSetupConfig] = None,
) -> Path:
    """Write .mdp templates + run_pipeline.sh + copies of gro/top/itp into `output_dir`.

    If `config.pre_solvation_stage` is set, its `extra_files` (e.g. a merged
    topology + inserted-molecule `.gro`/`.itp`) are copied into `output_dir`,
    and an `insert-molecules` stage plus its `solvation_top` are threaded
    through both the local `run_pipeline.sh` and (if `config.cluster == "slurm"`)
    the SLURM chain. `md_setup` stays domain-neutral: whatever built those files
    (e.g. a ligand/PFAS workflow) is the caller's concern.
    """
    cfg = config or MDSetupConfig()
    gro_path, top_path = Path(gro_path), Path(top_path)
    if not gro_path.exists():
        raise MDSetupError(f"gro file not found: {gro_path}")
    if not top_path.exists():
        raise MDSetupError(f"top file not found: {top_path}")

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    shutil.copy(gro_path, out / gro_path.name)
    shutil.copy(top_path, out / top_path.name)
    itp_path = top_path.with_suffix(".itp")
    if itp_path.exists():
        shutil.copy(itp_path, out / itp_path.name)

    _write_mdp_templates(out)

    stage = cfg.pre_solvation_stage
    if stage:
        for src in stage.extra_files:
            src_path = Path(src)
            if not src_path.exists():
                raise MDSetupError(
                    f"pre_solvation_stage.extra_files entry not found: {src_path}"
                )
            shutil.copy(src_path, out / src_path.name)

    ion = get_ion_profile(cfg.ion_profile)
    script = _render_pipeline_script(
        out, label, gro_path.name, top_path.name, cfg, ion,
        pre_solvation_stage=stage,
    )
    script_path = out / "run_pipeline.sh"
    script_path.write_text(script)
    script_path.chmod(script_path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)

    if cfg.cluster == "slurm":
        missing = [
            name for name, val in (
                ("slurm_submit_script", cfg.slurm_submit_script),
                ("slurm_sif_file", cfg.slurm_sif_file),
                ("slurm_partition", cfg.slurm_partition),
            ) if not val
        ]
        if missing:
            raise MDSetupError(
                "cluster='slurm' requires these site-specific MDSetupConfig fields "
                f"to be set for your cluster: {', '.join(missing)}."
            )
        solvate_ions_script = _render_solvate_ions_slurm(cfg, ion, pre_solvation_stage=stage)
        solvate_path = out / "solvate_ions.slurm"
        solvate_path.write_text(solvate_ions_script)
        solvate_path.chmod(solvate_path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)

        chain_script = _render_slurm_chain_script(
            out, label, gro_path.name, top_path.name, cfg, pre_solvation_stage=stage,
        )
        chain_path = out / "submit_chain.sh"
        chain_path.write_text(chain_script)
        chain_path.chmod(chain_path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    elif cfg.cluster is not None:
        raise MDSetupError(
            f"Unknown cluster {cfg.cluster!r}; only 'slurm' is currently implemented."
        )

    return out


def setup_md_from_manifest(
    manifest_csv: str | Path,
    output_root: str | Path,
    ion_profile: str = "mn_calcareous_default",
    config: Optional[MDSetupConfig] = None,
) -> list[dict]:
    """Read a `biochar.sweep` manifest.csv and write one MD run directory per row.

    Only rows whose `status` produced files (`strict_pass` or `fallback`) are
    processed; `skipped`/`failed` rows are reported but not written.

    Returns a list of per-row result dicts (also usable to build a summary
    table): {label, status, run_dir, gro_path, top_path, skipped_reason}.
    """
    manifest_csv = Path(manifest_csv)
    if not manifest_csv.exists():
        raise MDSetupError(f"manifest not found: {manifest_csv}")

    cfg = config or MDSetupConfig(ion_profile=ion_profile)
    # `sweep.run_sweep` records gro/top paths as given by `output_directory`
    # in the sweep config, which is typically relative to the *cwd the sweep
    # was run from* rather than to manifest.csv's own folder (manifest.csv
    # lives inside that same output_directory, so the two often coincide,
    # but not always -- e.g. a manifest copied elsewhere for archiving).
    # Try the path as-is (cwd-relative) first, then fall back to resolving
    # it against locations derived from the manifest's own path.
    candidate_bases = [Path("."), manifest_csv.parent, manifest_csv.parent.parent]
    out_root = Path(output_root)
    out_root.mkdir(parents=True, exist_ok=True)

    def _resolve(rel_path: str) -> Optional[Path]:
        if not rel_path:
            return None
        p = Path(rel_path)
        if p.is_absolute() and p.exists():
            return p
        for base in candidate_bases:
            candidate = base / rel_path
            if candidate.exists():
                return candidate
        return None

    results = []
    with manifest_csv.open(newline="") as fh:
        for row in csv.DictReader(fh):
            label = row.get("label") or f"row{row.get('index', '?')}"
            status = row.get("status", "")
            if status not in ("strict_pass", "fallback"):
                results.append(
                    {"label": label, "status": status, "run_dir": None,
                     "skipped_reason": f"no structure files (status={status})"}
                )
                continue

            gro_abs = _resolve(row.get("gro_path", ""))
            top_abs = _resolve(row.get("top_path", ""))
            if not gro_abs or not top_abs:
                results.append(
                    {"label": label, "status": status, "run_dir": None,
                     "skipped_reason": "gro/top file missing on disk"}
                )
                continue

            run_dir = out_root / label
            setup_one_structure(gro_abs, top_abs, run_dir, label=label, config=cfg)
            results.append(
                {"label": label, "status": status, "run_dir": str(run_dir),
                 "gro_path": str(gro_abs), "top_path": str(top_abs), "skipped_reason": None}
            )

    return results
