"""
Biochar Simulator Package

Tools for generating biochar molecular structures for GROMACS simulations.
"""

from .biochar_generator import (
    BiocharGenerator,
    GeneratorConfig,
    ValidationError,
    generate_biochar,
    generate_surface,
    generate_biochar_series,
)
from .heteroatom_assignment import CompositionResult, CompositionInfo
from .surface_builder import SurfaceBuilder, SurfaceConfig, SheetResult

__version__ = "0.1.2"
__all__ = [
    "BiocharGenerator",
    "GeneratorConfig",
    "ValidationError",
    "CompositionResult",
    "CompositionInfo",
    "generate_biochar",
    "generate_surface",
    "generate_biochar_series",
    "SurfaceBuilder",
    "SurfaceConfig",
    "SheetResult",
]
