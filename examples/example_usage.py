"""
Example: Generate Biochar Structures

Demonstrates how to use the BiocharGenerator to create biochar building blocks
for GROMACS molecular dynamics simulations.
"""

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from biochar_simulator.biochar_generator import BiocharGenerator, GeneratorConfig, generate_biochar


def example_1_basic_generation():
    """Example 1: Basic biochar generation with default parameters."""
    print("\n" + "=" * 70)
    print("EXAMPLE 1: Basic Biochar Generation")
    print("=" * 70)

    config = GeneratorConfig(
        target_num_carbons=50,
        H_C_ratio=0.5,
        O_C_ratio=0.1,
        aromaticity_percent=90.0,
        seed=42,  # For reproducibility
    )

    generator = BiocharGenerator(config)
    mol, coords, composition = generator.generate()
    generator.print_summary()

    # Export to GROMACS files
    generator.export_gromacs(
        output_directory="output",
        basename="biochar_50C",
    )


def example_2_different_compositions():
    """Example 2: Generate biochar with different compositions."""
    print("\n" + "=" * 70)
    print("EXAMPLE 2: Different Compositions")
    print("=" * 70)

    compositions = [
        {
            "name": "Low_O",
            "num_carbons": 75,
            "H_C": 0.4,
            "O_C": 0.05,
        },
        {
            "name": "Medium_O",
            "num_carbons": 100,
            "H_C": 0.6,
            "O_C": 0.15,
        },
        {
            "name": "High_O",
            "num_carbons": 60,
            "H_C": 0.8,
            "O_C": 0.30,
        },
    ]

    for comp in compositions:
        print(f"\nGenerating {comp['name']}...")

        config = GeneratorConfig(
            target_num_carbons=comp["num_carbons"],
            H_C_ratio=comp["H_C"],
            O_C_ratio=comp["O_C"],
            aromaticity_percent=85.0,
            seed=42,
        )

        generator = BiocharGenerator(config)
        mol, coords, composition = generator.generate()
        generator.print_summary()

        generator.export_gromacs(
            output_directory="output",
            basename=f"biochar_{comp['name']}",
        )


def example_3_convenience_function():
    """Example 3: Using the convenience function."""
    print("\n" + "=" * 70)
    print("EXAMPLE 3: Using Convenience Function")
    print("=" * 70)

    mol, coords, gro_path, top_path, itp_path = generate_biochar(
        target_num_carbons=100,
        H_C_ratio=0.55,
        O_C_ratio=0.12,
        aromaticity_percent=92.0,
        functional_groups={"phenolic": 3, "carboxyl": 1, "ether": 2},
        output_directory="output",
        basename="biochar_convenience",
        seed=123,
    )

    print(f"\nGenerated structure with {mol.GetNumAtoms()} atoms")
    print(f"Output files created successfully")


def example_4_large_structure():
    """Example 4: Generate a larger biochar structure."""
    print("\n" + "=" * 70)
    print("EXAMPLE 4: Large Biochar Structure")
    print("=" * 70)

    config = GeneratorConfig(
        target_num_carbons=500,
        H_C_ratio=0.45,
        O_C_ratio=0.08,
        aromaticity_percent=95.0,
        functional_groups={"phenolic": 5, "ether": 2},
        seed=999,
    )

    generator = BiocharGenerator(config)
    print("Generating large structure (may take a moment)...")
    mol, coords, composition = generator.generate()
    generator.print_summary()

    generator.export_gromacs(
        output_directory="output",
        basename="biochar_large_500C",
    )


def example_5_pentagon_defects():
    """Example 5: Biochar with pentagon ring defects."""
    print("\n" + "=" * 70)
    print("EXAMPLE 5: Biochar with Pentagon Ring Defects")
    print("=" * 70)

    config = GeneratorConfig(
        target_num_carbons=80,
        H_C_ratio=0.45,
        O_C_ratio=0.10,
        defect_fraction=0.15,  # 15% of rings will be pentagons
        aromaticity_percent=95.0,
        seed=555,
    )

    generator = BiocharGenerator(config)
    mol, coords, composition = generator.generate()
    print(f"Generated defective biochar with ~15% pentagonal rings")
    generator.print_summary()

    generator.export_gromacs(
        output_directory="output",
        basename="biochar_defects_15pct",
    )


def example_6_porous_surface():
    """Example 6: Porous slit-pore surface."""
    print("\n" + "=" * 70)
    print("EXAMPLE 6: Porous Slit-Pore Surface")
    print("=" * 70)

    from biochar_simulator.biochar_generator import generate_surface

    sheets, gro_path, top_path, itp_paths = generate_surface(
        target_num_carbons=40,
        H_C_ratio=0.3,
        O_C_ratio=0.05,
        functional_groups={"phenolic": 2, "ether": 1},
        pore_diameter=10.0,  # Angstroms
        num_sheets=2,
        output_directory="output",
        basename="surface_slit",
        seed=666,
    )

    print(f"\nGenerated {len(sheets)} sheets for slit-pore system")
    print(f"Pore diameter: 10.0 Å")
    print(f"Files: {gro_path.name}, {top_path.name}")
    for itp in itp_paths:
        print(f"  {itp.name}")


if __name__ == "__main__":
    # Create output directory
    Path("output").mkdir(exist_ok=True)

    # Run examples
    examples = [
        ("Example 1", example_1_basic_generation),
        ("Example 2", example_2_different_compositions),
        ("Example 3", example_3_convenience_function),
        ("Example 4", example_4_large_structure),
        ("Example 5", example_5_pentagon_defects),
        ("Example 6", example_6_porous_surface),
    ]

    for name, example_func in examples:
        try:
            example_func()
            print(f"\n✓ {name} completed successfully")
        except Exception as e:
            print(f"\n✗ {name} failed: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "=" * 70)
    print("All examples completed!")
    print("Check the 'output' directory for generated GROMACS files.")
    print("=" * 70)
