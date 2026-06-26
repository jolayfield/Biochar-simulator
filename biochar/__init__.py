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

__version__ = "0.3.0"
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
]
