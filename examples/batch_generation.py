"""
Batch Generation Example: Mixed Biochar Simulations

Demonstrates how to generate multiple biochar structures for mixed simulations.
Useful for temperature series, composition studies, or combined systems.
"""

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.biochar_generator import generate_biochar_series


def example_1_temperature_series():
    """Example 1: Temperature series (400, 600, 800°C)"""
    print("\n" + "=" * 70)
    print("EXAMPLE 1: Temperature Series Biochar")
    print("=" * 70)

    # Typical composition changes with pyrolysis temperature
    configs = [
        {
            "molecule_name": "BC400",
            "target_num_carbons": 80,
            "H_C_ratio": 0.65,  # More hydrogenated at lower temp
            "O_C_ratio": 0.20,  # More oxygen at lower temp
            "seed": 101,
        },
        {
            "molecule_name": "BC600",
            "target_num_carbons": 85,
            "H_C_ratio": 0.55,  # Moderate H content
            "O_C_ratio": 0.12,  # Lower O content
            "seed": 102,
        },
        {
            "molecule_name": "BC800",
            "target_num_carbons": 90,
            "H_C_ratio": 0.40,  # Less H at higher temp
            "O_C_ratio": 0.05,  # Much less O at high temp
            "seed": 103,
        },
    ]

    results = generate_biochar_series(
        configurations=configs,
        output_directory="output/temperature_series",
        create_combined_top=True,
        verbose=True,
    )

    print("\nGenerated files:")
    for mol_name, (gro, top, itp) in results.items():
        print(f"  {mol_name}: {gro.name}, {top.name}, {itp.name}")


def example_2_composition_series():
    """Example 2: Composition series (varying H/C ratio)"""
    print("\n" + "=" * 70)
    print("EXAMPLE 2: Composition Series (H/C Ratio)")
    print("=" * 70)

    configs = [
        {
            "molecule_name": "BCH04",
            "target_num_carbons": 70,
            "H_C_ratio": 0.4,
            "O_C_ratio": 0.1,
            "seed": 201,
        },
        {
            "molecule_name": "BCH06",
            "target_num_carbons": 70,
            "H_C_ratio": 0.6,
            "O_C_ratio": 0.1,
            "seed": 202,
        },
        {
            "molecule_name": "BCH08",
            "target_num_carbons": 70,
            "H_C_ratio": 0.8,
            "O_C_ratio": 0.1,
            "seed": 203,
        },
    ]

    results = generate_biochar_series(
        configurations=configs,
        output_directory="output/composition_series",
        create_combined_top=True,
        verbose=True,
    )

    print("\nGenerated files:")
    for mol_name, (gro, top, itp) in results.items():
        print(f"  {mol_name}: {gro.name}")


def example_3_oxygen_series():
    """Example 3: Oxygen content series"""
    print("\n" + "=" * 70)
    print("EXAMPLE 3: Oxygen Content Series")
    print("=" * 70)

    configs = [
        {
            "molecule_name": "BCO05",
            "target_num_carbons": 75,
            "H_C_ratio": 0.5,
            "O_C_ratio": 0.05,
            "seed": 301,
        },
        {
            "molecule_name": "BCO12",
            "target_num_carbons": 75,
            "H_C_ratio": 0.5,
            "O_C_ratio": 0.12,
            "seed": 302,
        },
        {
            "molecule_name": "BCO20",
            "target_num_carbons": 75,
            "H_C_ratio": 0.5,
            "O_C_ratio": 0.20,
            "seed": 303,
        },
    ]

    results = generate_biochar_series(
        configurations=configs,
        output_directory="output/oxygen_series",
        create_combined_top=True,
        verbose=True,
    )

    print("\nGenerated files:")
    for mol_name, (gro, top, itp) in results.items():
        print(f"  {mol_name}: {gro.name}")


def example_4_size_series():
    """Example 4: Size series (different carbon counts)"""
    print("\n" + "=" * 70)
    print("EXAMPLE 4: Size Series (Carbon Count)")
    print("=" * 70)

    configs = [
        {
            "molecule_name": "BC10",
            "target_num_carbons": 20,
            "H_C_ratio": 0.5,
            "O_C_ratio": 0.1,
            "seed": 401,
        },
        {
            "molecule_name": "BC20",
            "target_num_carbons": 50,
            "H_C_ratio": 0.5,
            "O_C_ratio": 0.1,
            "seed": 402,
        },
        {
            "molecule_name": "BC30",
            "target_num_carbons": 100,
            "H_C_ratio": 0.5,
            "O_C_ratio": 0.1,
            "seed": 403,
        },
    ]

    results = generate_biochar_series(
        configurations=configs,
        output_directory="output/size_series",
        create_combined_top=True,
        verbose=True,
    )

    print("\nGenerated files:")
    for mol_name, (gro, top, itp) in results.items():
        print(f"  {mol_name}: {gro.name}")


def example_5_defect_structures():
    """Example 5: Defect structures (pentagon ring insertion)"""
    print("\n" + "=" * 70)
    print("EXAMPLE 5: Defect Structures (Pentagon Ring Insertion)")
    print("=" * 70)
    print("Comparing pure hexagonal vs defective biochar:")
    print("- BC800: Pure hexagonal (defect_fraction=0.0)")
    print("- BD805: 5% pentagonal defects (defect_fraction=0.05)")
    print("- BD815: 15% pentagonal defects (defect_fraction=0.15)")

    configs = [
        {
            "molecule_name": "BC800",
            "target_num_carbons": 80,
            "H_C_ratio": 0.40,
            "O_C_ratio": 0.05,
            "defect_fraction": 0.0,
            "seed": 601,
        },
        {
            "molecule_name": "BD805",
            "target_num_carbons": 80,
            "H_C_ratio": 0.40,
            "O_C_ratio": 0.05,
            "defect_fraction": 0.05,  # 5% pentagons
            "seed": 602,
        },
        {
            "molecule_name": "BD815",
            "target_num_carbons": 80,
            "H_C_ratio": 0.40,
            "O_C_ratio": 0.05,
            "defect_fraction": 0.15,  # 15% pentagons
            "seed": 603,
        },
    ]

    results = generate_biochar_series(
        configurations=configs,
        output_directory="output/defect_study",
        create_combined_top=True,
        verbose=True,
    )

    print("\nGenerated files:")
    for mol_name, (gro, top, itp) in results.items():
        print(f"  {mol_name}: {gro.name}")

    print("\n✓ Ready for studying effect of topological disorder!")


def example_6_custom_mixed_system():
    """Example 6: Custom mixed system (different biochar types together)"""
    print("\n" + "=" * 70)
    print("EXAMPLE 6: Custom Mixed System")
    print("=" * 70)
    print("Simulating a system with multiple biochar types:")
    print("- BC400: Low-temperature biochar (high O/C)")
    print("- BC800: High-temperature biochar (low O/C)")
    print("- BCHWP: Chemically oxidized biochar")

    configs = [
        {
            "molecule_name": "BC400",
            "target_num_carbons": 60,
            "H_C_ratio": 0.65,
            "O_C_ratio": 0.25,  # Highly oxidized
            "seed": 701,
        },
        {
            "molecule_name": "BC800",
            "target_num_carbons": 100,
            "H_C_ratio": 0.35,
            "O_C_ratio": 0.03,  # Minimal oxidation
            "seed": 702,
        },
        {
            "molecule_name": "BCHWP",
            "target_num_carbons": 80,
            "H_C_ratio": 0.5,
            "O_C_ratio": 0.30,  # Chemically oxidized
            "seed": 703,
        },
    ]

    results = generate_biochar_series(
        configurations=configs,
        output_directory="output/mixed_system",
        create_combined_top=True,
        verbose=True,
    )

    print("\nGenerated files:")
    for mol_name, (gro, top, itp) in results.items():
        print(f"  {mol_name}: {gro.name}")

    print("\n✓ Ready for GROMACS simulation with mixed biochar types!")
    print("  Use: gmx grompp -f md.mdp -c combined.gro -p combined.top")


if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("BATCH GENERATION EXAMPLES - Mixed Biochar Simulations")
    print("=" * 70)

    examples = [
        ("Example 1: Temperature Series", example_1_temperature_series),
        ("Example 2: Composition Series", example_2_composition_series),
        ("Example 3: Oxygen Series", example_3_oxygen_series),
        ("Example 4: Size Series", example_4_size_series),
        ("Example 5: Defect Structures", example_5_defect_structures),
        ("Example 6: Custom Mixed System", example_6_custom_mixed_system),
    ]

    completed = 0
    for name, example_func in examples:
        try:
            example_func()
            completed += 1
        except Exception as e:
            print(f"\n✗ {name} failed: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "=" * 70)
    print(f"COMPLETED: {completed}/{len(examples)} examples")
    print("=" * 70)
    print("\nOutput directory structure:")
    print("output/")
    print("├── temperature_series/  (BC400, BC600, BC800)")
    print("├── composition_series/  (BCH04, BCH06, BCH08)")
    print("├── oxygen_series/       (BCO05, BCO12, BCO20)")
    print("├── size_series/         (BC10, BC20, BC30)")
    print("├── defect_study/        (BC800, BD805, BD815 + combined.top)")
    print("└── mixed_system/        (BC400, BC800, BCHWP + combined.top)")
    print("\nEach directory contains:")
    print("  - Individual .gro structure files")
    print("  - Individual .top topology files")
    print("  - Individual .itp molecule definitions")
    print("  - combined.top (for running mixed simulations)")
    print("=" * 70 + "\n")
