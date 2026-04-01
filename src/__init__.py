"""
Biochar Simulator Package

Tools for generating biochar molecular structures for GROMACS simulations.
"""

from .biochar_generator import BiocharGenerator, GeneratorConfig, generate_biochar

__version__ = "0.1.0"
__all__ = ["BiocharGenerator", "GeneratorConfig", "generate_biochar"]
