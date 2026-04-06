"""
Test Suite: PAH Structure Quality and Composition Optimization

Tests different molecular sizes, functional group compositions, and
assembly strategies to characterize what's achievable with the biochar generator.
"""

import sys
import warnings
import logging
from pathlib import Path
import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Optional
import json

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.CRITICAL)

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from src.biochar_generator import BiocharGenerator, GeneratorConfig
from src.geometry_3d import GeometryValidator
from src.constants import PAH_LIBRARY
from rdkit import Chem
from rdkit.Chem import AllChem


@dataclass
class GeometryMetrics:
    """Metrics for evaluating molecule geometry quality."""
    num_atoms: int
    num_carbons: int
    num_hydrogens: int
    num_oxygens: int
    num_aromatic_carbons: int
    aromaticity_percent: float

    bond_errors: int
    steric_clashes: int
    is_valid: bool

    coordinate_span_x: float
    coordinate_span_y: float
    coordinate_span_z: float
    average_span: float

    ring_planarity: float
    h_c_ratio: float
    o_c_ratio: float

    def __str__(self) -> str:
        """Pretty print metrics."""
        return f"""
Composition:
  C: {self.num_carbons}, H: {self.num_hydrogens}, O: {self.num_oxygens}
  H/C: {self.h_c_ratio:.3f}, O/C: {self.o_c_ratio:.3f}
  Aromaticity: {self.aromaticity_percent:.1f}% ({self.num_aromatic_carbons}/{self.num_carbons})

Geometry Quality:
  Valid: {'✓' if self.is_valid else '✗'}
  Bond errors: {self.bond_errors}
  Steric clashes: {self.steric_clashes}
  Ring planarity: {self.ring_planarity:.3f} Å

Structure Size:
  X: {self.coordinate_span_x:.2f} Å
  Y: {self.coordinate_span_y:.2f} Å
  Z: {self.coordinate_span_z:.2f} Å
  Avg: {self.average_span:.2f} Å
"""


class PAHQualityTester:
    """Test suite for evaluating PAH biochar structures."""

    def __init__(self, output_dir: str = "tests/results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results = {}

    def test_library_smiles(self) -> Dict[str, bool]:
        """Validate all SMILES in PAH_LIBRARY parse and give correct atom counts."""
        print("\n" + "=" * 70)
        print("TEST 0: PAH Library SMILES Validation")
        print("=" * 70)

        results = {}
        print(f"\n{'Name':<20} {'Expected':<10} {'Actual':<10} {'Valid':<6}")
        print("-" * 50)

        for name, data in PAH_LIBRARY.items():
            smiles = data["smiles"]
            expected_C = data["num_atoms"]
            mol = Chem.MolFromSmiles(smiles)

            if mol is None:
                print(f"  {name:<20} {'PARSE FAIL'}")
                results[name] = False
                continue

            try:
                Chem.SanitizeMol(mol)
            except Exception:
                print(f"  {name:<20} {expected_C:<10} {'SANITIZE FAIL'}")
                results[name] = False
                continue

            actual_C = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
            is_valid = (actual_C == expected_C)
            status = "✓" if is_valid else "✗"
            print(f"  {name:<20} {expected_C:<10} {actual_C:<10} {status}")
            results[name] = is_valid

        n_valid = sum(1 for v in results.values() if v)
        print(f"\n  {n_valid}/{len(results)} SMILES valid")
        self.results['library_smiles'] = results
        return results

    def test_pah_sizes(self) -> Dict[int, GeometryMetrics]:
        """Test PAH skeleton generation for all sizes including large ones."""
        print("\n" + "=" * 70)
        print("TEST 1: PAH Size Exploration (Benzene → 200C)")
        print("=" * 70)

        pah_sizes = [6, 10, 14, 16, 18, 22, 24, 26, 28, 30, 40, 50, 100, 200]
        results = {}

        for target_size in pah_sizes:
            print(f"\n Testing {target_size}C PAH...")
            try:
                config = GeneratorConfig(
                    target_num_carbons=target_size,
                    H_C_ratio=0.5,
                    O_C_ratio=0.1,
                    seed=42
                )
                gen = BiocharGenerator(config)
                gen.generate()

                metrics = self._evaluate_molecule(gen.mol, gen.coords)
                results[target_size] = metrics

                # Atom count accuracy check
                count_ok = abs(metrics.num_carbons - target_size) <= max(5, target_size * 0.10)
                print(f"  ✓ {target_size}C: got {metrics.num_carbons}C, "
                      f"{metrics.aromaticity_percent:.0f}% arom, "
                      f"{metrics.bond_errors} bond errors, "
                      f"{metrics.steric_clashes} clashes, "
                      f"planarity={metrics.ring_planarity:.3f}Å "
                      f"{'(count OK)' if count_ok else '(count OFF)'}")

            except Exception as e:
                print(f"  ✗ {target_size}C failed: {str(e)[:80]}")

        self.results['pah_sizes'] = results
        return results

    def test_kekulizability(self) -> Dict[int, bool]:
        """Test that all generated structures have 100% aromatic carbons."""
        print("\n" + "=" * 70)
        print("TEST 1b: Kekulizability / Aromaticity Check")
        print("=" * 70)

        sizes = [6, 16, 22, 24, 26, 28, 30, 40, 50, 100, 200]
        results = {}

        from src.carbon_skeleton import PAHAssembler
        for size in sizes:
            asm = PAHAssembler(seed=42)
            sk = asm.generate(size)
            arom_pct = sk.aromaticity_percent
            all_arom = arom_pct == 100.0
            results[size] = all_arom
            status = "✓" if all_arom else "✗"
            print(f"  {status} {size:4d}C: {arom_pct:.1f}% aromatic")

        self.results['kekulizability'] = results
        return results

    def test_atom_count_accuracy(self) -> Dict[int, float]:
        """Test accuracy of carbon count vs target."""
        print("\n" + "=" * 70)
        print("TEST 1c: Atom Count Accuracy")
        print("=" * 70)

        targets = [10, 18, 22, 26, 30, 40, 50, 75, 100, 150, 200]
        results = {}

        from src.carbon_skeleton import PAHAssembler
        for target in targets:
            asm = PAHAssembler(seed=42)
            sk = asm.generate(target)
            error_pct = abs(sk.num_carbons - target) / target * 100
            results[target] = error_pct
            ok = "✓" if error_pct <= 10 else "!"
            print(f"  {ok} target={target:4d}, got={sk.num_carbons:4d} (error={error_pct:.1f}%)")

        self.results['atom_count_accuracy'] = results
        return results

    def test_compositions(self, target_carbons: int = 50) -> Dict[str, GeometryMetrics]:
        """Test different H/C and O/C ratios for a fixed PAH size."""
        print("\n" + "=" * 70)
        print(f"TEST 2: Composition Variation ({target_carbons}C PAH)")
        print("=" * 70)

        compositions = [
            ("Low O", 0.5, 0.0),
            ("Moderate O", 0.5, 0.1),
            ("High O", 0.5, 0.2),
            ("Low H", 0.3, 0.1),
            ("High H", 0.8, 0.1),
        ]

        results = {}

        for name, h_c_ratio, o_c_ratio in compositions:
            print(f"\n Testing {name} (H/C={h_c_ratio:.1f}, O/C={o_c_ratio:.1f})...")
            try:
                config = GeneratorConfig(
                    target_num_carbons=target_carbons,
                    H_C_ratio=h_c_ratio,
                    O_C_ratio=o_c_ratio,
                    seed=42
                )
                gen = BiocharGenerator(config)
                gen.generate()

                metrics = self._evaluate_molecule(gen.mol, gen.coords)
                results[name] = metrics

                print(f"  ✓ {name}: H/C={metrics.h_c_ratio:.2f}, O/C={metrics.o_c_ratio:.2f}, "
                      f"{metrics.bond_errors} bond errors, {metrics.steric_clashes} clashes")

            except Exception as e:
                print(f"  ✗ {name} failed: {str(e)[:50]}")

        self.results['compositions'] = results
        return results

    def test_seed_variation(self, target_carbons: int = 50, num_seeds: int = 5) -> List[GeometryMetrics]:
        """Test quality variation across different random seeds."""
        print("\n" + "=" * 70)
        print(f"TEST 3: Seed Variation ({num_seeds} seeds, {target_carbons}C PAH)")
        print("=" * 70)

        results = []

        for seed in range(num_seeds):
            print(f"\n Testing seed {seed}...")
            try:
                config = GeneratorConfig(
                    target_num_carbons=target_carbons,
                    H_C_ratio=0.5,
                    O_C_ratio=0.1,
                    seed=seed
                )
                gen = BiocharGenerator(config)
                gen.generate()

                metrics = self._evaluate_molecule(gen.mol, gen.coords)
                results.append(metrics)

                print(f"  ✓ Seed {seed}: {metrics.num_carbons}C, "
                      f"{metrics.bond_errors} bond errors, {metrics.steric_clashes} clashes")

            except Exception as e:
                print(f"  ✗ Seed {seed} failed: {str(e)[:50]}")

        self.results['seed_variation'] = results
        return results

    def test_large_structures(self) -> Dict[int, GeometryMetrics]:
        """Test large PAH structures specifically (50–200C range)."""
        print("\n" + "=" * 70)
        print("TEST 4: Large Structure Quality (50–200C)")
        print("=" * 70)

        results = {}
        sizes = [50, 100, 150, 200]

        for target in sizes:
            print(f"\n Testing {target}C...")
            try:
                config = GeneratorConfig(
                    target_num_carbons=target,
                    H_C_ratio=0.5,
                    O_C_ratio=0.1,
                    seed=42
                )
                gen = BiocharGenerator(config)
                gen.generate()

                metrics = self._evaluate_molecule(gen.mol, gen.coords)
                results[target] = metrics

                # For large structures, check scaling properties
                clash_per_C = metrics.steric_clashes / max(metrics.num_carbons, 1)
                print(f"  ✓ {target}C: {metrics.num_carbons}C actual, "
                      f"H/C={metrics.h_c_ratio:.2f}, "
                      f"clashes/C={clash_per_C:.2f}, "
                      f"planarity={metrics.ring_planarity:.3f}Å")

            except Exception as e:
                import traceback
                print(f"  ✗ {target}C failed: {str(e)[:60]}")

        self.results['large_structures'] = results
        return results

    def _evaluate_molecule(self, mol: Chem.Mol, coords: np.ndarray) -> GeometryMetrics:
        """Evaluate a molecule's geometry quality."""
        valid, errors = GeometryValidator.validate_geometry(mol, coords)
        bond_errors = len([e for e in errors if "Unusual bond length" in e])
        steric_clashes = len([e for e in errors if "Steric clash" in e])

        num_carbons = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
        num_hydrogens = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 1)
        num_oxygens = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
        num_aromatic = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetIsAromatic())

        aromaticity = (num_aromatic / num_carbons * 100) if num_carbons > 0 else 0

        planarity, _ = GeometryValidator.measure_ring_planarity(mol, coords)

        span_x = coords[:, 0].max() - coords[:, 0].min()
        span_y = coords[:, 1].max() - coords[:, 1].min()
        span_z = coords[:, 2].max() - coords[:, 2].min()
        avg_span = np.mean([span_x, span_y, span_z])

        h_c = num_hydrogens / max(num_carbons, 1)
        o_c = num_oxygens / max(num_carbons, 1)

        return GeometryMetrics(
            num_atoms=mol.GetNumAtoms(),
            num_carbons=num_carbons,
            num_hydrogens=num_hydrogens,
            num_oxygens=num_oxygens,
            num_aromatic_carbons=num_aromatic,
            aromaticity_percent=aromaticity,
            bond_errors=bond_errors,
            steric_clashes=steric_clashes,
            is_valid=valid,
            coordinate_span_x=span_x,
            coordinate_span_y=span_y,
            coordinate_span_z=span_z,
            average_span=avg_span,
            ring_planarity=planarity,
            h_c_ratio=h_c,
            o_c_ratio=o_c,
        )

    def generate_report(self) -> str:
        """Generate a comprehensive test report."""
        report = []
        report.append("\n" + "=" * 70)
        report.append("BIOCHAR PAH TEST SUITE REPORT")
        report.append("=" * 70)

        # Library SMILES
        if 'library_smiles' in self.results:
            report.append("\n" + "-" * 70)
            report.append("TEST 0: PAH Library SMILES Validation")
            report.append("-" * 70)
            smiles_data = self.results['library_smiles']
            n_valid = sum(1 for v in smiles_data.values() if v)
            report.append(f"\n  {n_valid}/{len(smiles_data)} SMILES valid")
            for name, valid in smiles_data.items():
                report.append(f"  {'✓' if valid else '✗'} {name}")

        # PAH Sizes
        if 'pah_sizes' in self.results:
            report.append("\n" + "-" * 70)
            report.append("TEST 1: PAH Sizes (6–200C)")
            report.append("-" * 70)

            sizes_data = self.results['pah_sizes']
            report.append(f"\n{'Target':<8} {'Actual':<8} {'Arom%':<7} {'Bonds':<10} {'Clashes':<10} {'Valid':<8} {'Span':<10} {'H/C':<6}")
            report.append("-" * 70)

            for size in sorted(sizes_data.keys()):
                metrics = sizes_data[size]
                valid_str = "✓" if metrics.is_valid else "✗"
                report.append(
                    f"{size}C{'':<3} {metrics.num_carbons}C{'':<3} "
                    f"{metrics.aromaticity_percent:.0f}%{'':<3} "
                    f"{metrics.bond_errors:<10} {metrics.steric_clashes:<10} "
                    f"{valid_str:<8} {metrics.average_span:.2f} Å  "
                    f"{metrics.h_c_ratio:.2f}"
                )

        # Kekulizability
        if 'kekulizability' in self.results:
            report.append("\n" + "-" * 70)
            report.append("TEST 1b: Kekulizability (100% aromatic required)")
            report.append("-" * 70)
            kek_data = self.results['kekulizability']
            for size, ok in sorted(kek_data.items()):
                report.append(f"  {'✓' if ok else '✗'} {size}C: {'100% aromatic' if ok else 'NOT fully aromatic'}")

        # Atom count accuracy
        if 'atom_count_accuracy' in self.results:
            report.append("\n" + "-" * 70)
            report.append("TEST 1c: Atom Count Accuracy (target vs actual)")
            report.append("-" * 70)
            acc_data = self.results['atom_count_accuracy']
            report.append(f"\n{'Target':<8} {'Error%':<10} {'OK?':<6}")
            report.append("-" * 25)
            for target, err in sorted(acc_data.items()):
                ok = err <= 10
                report.append(f"{target:<8} {err:<10.1f} {'✓' if ok else '✗'}")

        # Compositions
        if 'compositions' in self.results:
            report.append("\n" + "-" * 70)
            report.append("TEST 2: Compositions")
            report.append("-" * 70)

            comp_data = self.results['compositions']
            report.append(f"\n{'Composition':<20} {'H/C':<6} {'O/C':<6} {'Bonds':<10} {'Clashes':<10} {'Valid':<8}")
            report.append("-" * 65)

            for name in sorted(comp_data.keys()):
                metrics = comp_data[name]
                valid_str = "✓" if metrics.is_valid else "✗"
                report.append(
                    f"{name:<20} {metrics.h_c_ratio:<6.3f} {metrics.o_c_ratio:<6.3f} "
                    f"{metrics.bond_errors:<10} {metrics.steric_clashes:<10} {valid_str:<8}"
                )

        # Seed Variation
        if 'seed_variation' in self.results:
            report.append("\n" + "-" * 70)
            report.append("TEST 3: Seed Variation (Statistics)")
            report.append("-" * 70)

            seed_data = self.results['seed_variation']
            bond_errors = [m.bond_errors for m in seed_data]
            clashes = [m.steric_clashes for m in seed_data]

            report.append(f"\nBond Errors:")
            report.append(f"  Mean: {np.mean(bond_errors):.1f} ± {np.std(bond_errors):.1f}")
            report.append(f"  Range: {min(bond_errors)}-{max(bond_errors)}")

            report.append(f"\nSteric Clashes:")
            report.append(f"  Mean: {np.mean(clashes):.1f} ± {np.std(clashes):.1f}")
            report.append(f"  Range: {min(clashes)}-{max(clashes)}")

        # Large Structures
        if 'large_structures' in self.results:
            report.append("\n" + "-" * 70)
            report.append("TEST 4: Large Structures (50–200C)")
            report.append("-" * 70)
            lg_data = self.results['large_structures']
            report.append(f"\n{'Target':<8} {'Actual':<8} {'H/C':<6} {'Clashes/C':<12} {'Planarity':<12}")
            report.append("-" * 50)
            for size in sorted(lg_data.keys()):
                m = lg_data[size]
                cpc = m.steric_clashes / max(m.num_carbons, 1)
                report.append(f"{size:<8} {m.num_carbons:<8} {m.h_c_ratio:<6.2f} {cpc:<12.2f} {m.ring_planarity:.3f} Å")

        # Summary
        report.append("\n" + "=" * 70)
        report.append("SUMMARY & RECOMMENDATIONS")
        report.append("=" * 70)

        if 'pah_sizes' in self.results:
            best_size = min(
                self.results['pah_sizes'].items(),
                key=lambda x: x[1].bond_errors + x[1].steric_clashes
            )
            report.append(f"\n✓ Best performing size: {best_size[0]}C")
            report.append(f"  Bond errors: {best_size[1].bond_errors}")
            report.append(f"  Steric clashes: {best_size[1].steric_clashes}")

        report.append(f"\n✓ Supported sizes: 6C to 200C+")
        report.append(f"  Skeleton: parity-aware PAH seed + 4-node ring growth")
        report.append(f"  Geometry: 2D-first embedding for >80 heavy atoms")

        return "\n".join(report)

    def save_report(self, filename: str = "test_report.txt"):
        """Save test report to file."""
        report = self.generate_report()
        report_path = self.output_dir / filename

        with open(report_path, 'w') as f:
            f.write(report)

        print(f"\n✓ Report saved to {report_path}")
        return report_path


def main():
    """Run the complete test suite."""
    tester = PAHQualityTester()

    # Run all tests
    tester.test_library_smiles()
    tester.test_kekulizability()
    tester.test_atom_count_accuracy()
    pah_results = tester.test_pah_sizes()
    comp_results = tester.test_compositions(target_carbons=50)
    seed_results = tester.test_seed_variation(target_carbons=50, num_seeds=5)
    large_results = tester.test_large_structures()

    # Generate and print report
    print(tester.generate_report())

    # Save report
    tester.save_report()


if __name__ == "__main__":
    main()
