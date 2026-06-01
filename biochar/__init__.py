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
    GeneratorConfig,
    ValidationError,
    generate_biochar,
    generate_biochar_series,
    generate_surface,
    generate_biochar_series,
)
from .heteroatom_assignment import CompositionResult, CompositionInfo
from .surface_builder import SurfaceBuilder, SurfaceConfig, SheetResult

__version__ = "0.1.4"
__all__ = [
    "BiocharGenerator",
    "GeneratorConfig",
    "ValidationError",
    "CompositionResult",
    "CompositionInfo",
    "generate_biochar",
    "generate_biochar_series",
    "generate_surface",
    "generate_biochar_series",
    "SurfaceBuilder",
    "SurfaceConfig",
    "SheetResult",
]
