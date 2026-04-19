"""
Biochar Simulator Package

Tools for generating biochar molecular structures for GROMACS simulations.
"""

from .biochar_generator import (
    BiocharGenerator,
    GeneratorConfig,
    generate_biochar,
    generate_surface,
)
from .surface_builder import SurfaceBuilder, SurfaceConfig, SheetResult

__version__ = "0.1.1"
__all__ = [
    "BiocharGenerator",
    "GeneratorConfig",
    "generate_biochar",
    "generate_surface",
    "SurfaceBuilder",
    "SurfaceConfig",
    "SheetResult",
]
