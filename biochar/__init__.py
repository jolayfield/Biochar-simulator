"""
Biochar Simulator Package

Tools for generating biochar molecular structures for GROMACS simulations.
"""

import logging

# Standard library best practice: add NullHandler so the package never emits
# "No handler found" warnings when used without logging configured.
logging.getLogger(__name__).addHandler(logging.NullHandler())

from .biochar_generator import (
    BiocharGenerator,
    BiocharResult,
    GeneratorConfig,
    ValidationError,
    generate_biochar,
    generate_biochar_series,
    generate_surface,
)
from .heteroatom_assignment import CompositionResult, CompositionInfo
from .surface_builder import SurfaceBuilder, SurfaceConfig, SheetResult
from .temperature_model import TemperatureModel, properties, VALID_FEEDSTOCKS
from .qm_charges import QMChargeError
from .sweep import (
    run_sweep,
    expand_grid,
    build_point,
    load_sweep_config,
    GridPoint,
    PointResult,
    SweepError,
)
from .md_setup import (
    setup_md_from_manifest,
    setup_one_structure,
    MDSetupConfig,
    MDSetupError,
    IonProfile,
    ION_PROFILES,
    get_ion_profile,
    PreSolvationStage,
    MoleculeInsertion,
)
from .condensation import (
    AnnealSpec,
    CondensationError,
    anneal_spec_for_htt,
    setup_condensation,
    generate_and_condense,
    setup_surface,
    add_surface_and_validation,
    write_validation_setup,
    estimate_box_nm,
)

__version__ = "0.4.0"
__all__ = [
    "BiocharGenerator",
    "BiocharResult",
    "GeneratorConfig",
    "ValidationError",
    "QMChargeError",
    "CompositionResult",
    "CompositionInfo",
    "generate_biochar",
    "generate_biochar_series",
    "generate_surface",
    "SurfaceBuilder",
    "SurfaceConfig",
    "SheetResult",
    "TemperatureModel",
    "properties",
    "VALID_FEEDSTOCKS",
    "run_sweep",
    "expand_grid",
    "build_point",
    "load_sweep_config",
    "GridPoint",
    "PointResult",
    "SweepError",
    "setup_md_from_manifest",
    "setup_one_structure",
    "MDSetupConfig",
    "MDSetupError",
    "IonProfile",
    "ION_PROFILES",
    "get_ion_profile",
    "PreSolvationStage",
    "MoleculeInsertion",
    "AnnealSpec",
    "CondensationError",
    "anneal_spec_for_htt",
    "setup_condensation",
    "generate_and_condense",
    "setup_surface",
    "add_surface_and_validation",
    "write_validation_setup",
    "estimate_box_nm",
]
