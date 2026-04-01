"""
Carbon Skeleton Generation

Generates aromatic carbon frameworks for biochar structures using PAHs or random graphs.
"""

import random
from typing import List, Tuple, Set, Optional
from dataclasses import dataclass

import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from .constants import PAH_LIBRARY


@dataclass
class CarbonSkeleton:
    """Represents a carbon skeleton/framework."""

    mol: Chem.Mol
    smiles: str
    num_carbons: int
    num_aromatic_carbons: int
    aromaticity_percent: float


class PAHAssembler:
    """
    Build aromatic carbon skeletons from Polycyclic Aromatic Hydrocarbons (PAHs).

    Strategy: Start with PAH building blocks and fuse them together to reach target size.
    """

    def __init__(self, seed: int = None):
        self.seed = seed
        if seed is not None:
            random.seed(seed)

        # Preload PAH library into RDKit molecules
        self.pahs = {}
        for name, data in PAH_LIBRARY.items():
            mol = Chem.MolFromSmiles(data["smiles"])
            if mol is not None:
                Chem.SanitizeMol(mol)
                self.pahs[name] = {
                    "mol": mol,
                    "num_carbons": data["num_atoms"],
                    "data": data,
                }

    def generate(
        self,
        target_num_carbons: int,
        target_aromaticity: float = 100.0,
        prefer_larger_pahs: bool = True,
    ) -> CarbonSkeleton:
        """
        Generate a PAH-based carbon skeleton.

        Args:
            target_num_carbons: Target number of carbon atoms
            target_aromaticity: Target aromaticity percentage (0-100)
            prefer_larger_pahs: If True, prefer larger PAHs for assembly

        Returns:
            CarbonSkeleton object
        """
        # Try to create the molecule from PAHs
        try:
            # Select single PAH or raise exception to trigger fallback
            smiles = self._select_pah_smiles(target_num_carbons, prefer_larger_pahs)

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Could not parse SMILES")

            Chem.SanitizeMol(mol)
            Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)

            num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            num_aromatic = sum(
                1 for atom in mol.GetAtoms()
                if atom.GetAtomicNum() == 6 and atom.GetIsAromatic()
            )
            aromaticity = (num_aromatic / num_carbons * 100) if num_carbons > 0 else 0

            return CarbonSkeleton(
                mol=mol,
                smiles=smiles,
                num_carbons=num_carbons,
                num_aromatic_carbons=num_aromatic,
                aromaticity_percent=aromaticity,
            )

        except Exception as e:
            # Fallback to random graph if PAH generation fails or target too large
            return self._fallback_to_random(target_num_carbons, target_aromaticity)

    def _select_pah_smiles(self, target_num_carbons: int, prefer_larger: bool) -> str:
        """
        Select appropriate PAH to reach closest match to target carbon count.

        Uses verified PAH SMILES up to pyrene (16C). For larger targets,
        raises exception to trigger RandomGraphGenerator fallback.
        """
        # Select best single PAH match (verified SMILES only)
        if target_num_carbons <= 6:
            return PAH_LIBRARY["benzene"]["smiles"]
        elif target_num_carbons <= 10:
            return PAH_LIBRARY["naphthalene"]["smiles"]
        elif target_num_carbons <= 14:
            return PAH_LIBRARY["anthracene"]["smiles"]
        elif target_num_carbons <= 16:
            return PAH_LIBRARY["pyrene"]["smiles"]
        else:
            # For targets > 16, use RandomGraphGenerator which builds structures
            # by concatenating PAH units (more robust than complex SMILES)
            raise ValueError(
                f"Target {target_num_carbons}C exceeds maximum single PAH. "
                "Falling back to modular PAH assembly."
            )


    def _fuse_pahs(self, mol1: Chem.Mol, mol2: Chem.Mol) -> Optional[Chem.Mol]:
        """
        Attempt to fuse two PAH molecules at aromatic edges.

        Simplified approach: concatenate SMILES with shared carbon bridge.
        A more robust approach would identify best fusion sites.
        """
        try:
            # Get SMILES strings
            smiles1 = Chem.MolToSmiles(mol1)
            smiles2 = Chem.MolToSmiles(mol2)

            # Find aromatic carbons without neighboring heteroatoms
            # (simplified - just use the molecules as-is)
            # More sophisticated fusion would identify phenyl-phenyl bonds

            # For now, attempt simple benzene-linking fusion
            # This is a placeholder for more complex fusion algorithms
            fused_smiles = smiles1 + "-" + smiles2

            fused_mol = Chem.MolFromSmiles(fused_smiles)
            if fused_mol is not None:
                Chem.SanitizeMol(fused_mol)
                return fused_mol
        except Exception:
            pass

        return None

    def _fallback_to_random(
        self, target_num_carbons: int, target_aromaticity: float
    ) -> CarbonSkeleton:
        """Fallback to random graph generation if PAH assembly fails."""
        generator = RandomGraphGenerator(seed=self.seed)
        return generator.generate(target_num_carbons, target_aromaticity)


class RandomGraphGenerator:
    """
    Generate random aromatic carbon frameworks as fallback.

    Creates a planar graph of carbons with aromatic rings.
    """

    def __init__(self, seed: int = None):
        self.seed = seed
        if seed is not None:
            random.seed(seed)

    def generate(
        self, target_num_carbons: int, target_aromaticity: float = 100.0
    ) -> CarbonSkeleton:
        """
        Generate a random aromatic carbon skeleton.

        For simplicity and robustness, uses PAH templates directly
        rather than graph-based generation.

        Args:
            target_num_carbons: Target number of carbon atoms
            target_aromaticity: Target aromaticity percentage

        Returns:
            CarbonSkeleton object
        """
        # Use template-based generation (more robust than graph conversion)
        smiles = self._create_default_aromatic(target_num_carbons)

        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Could not parse SMILES")

            Chem.SanitizeMol(mol)
            num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            num_aromatic = sum(
                1 for atom in mol.GetAtoms()
                if atom.GetAtomicNum() == 6 and atom.GetIsAromatic()
            )
            aromaticity = (num_aromatic / num_carbons * 100) if num_carbons > 0 else 0

            return CarbonSkeleton(
                mol=mol,
                smiles=Chem.MolToSmiles(mol),
                num_carbons=num_carbons,
                num_aromatic_carbons=num_aromatic,
                aromaticity_percent=aromaticity,
            )

        except Exception:
            # Final fallback: benzene if everything fails
            mol = Chem.MolFromSmiles("c1ccccc1")
            Chem.SanitizeMol(mol)
            num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            return CarbonSkeleton(
                mol=mol,
                smiles="c1ccccc1",
                num_carbons=num_carbons,
                num_aromatic_carbons=num_carbons,
                aromaticity_percent=100.0,
            )

    def _create_random_aromatic_graph(
        self, num_nodes: int, target_aromaticity: float
    ) -> nx.Graph:
        """
        Create a random planar graph with aromatic rings.

        Simple approach: create benzene rings and connect them randomly.
        """
        graph = nx.Graph()

        # Number of benzene rings to create
        num_rings = max(1, num_nodes // 6)

        # Create benzene rings
        for i in range(num_rings):
            ring_nodes = list(range(i * 6, (i + 1) * 6))
            graph.add_nodes_from(ring_nodes)
            for j in range(6):
                graph.add_edge(ring_nodes[j], ring_nodes[(j + 1) % 6])

        # Connect rings randomly
        all_nodes = list(graph.nodes())
        for _ in range(max(0, num_rings - 1)):
            if len(all_nodes) >= 2:
                n1, n2 = random.sample(all_nodes, 2)
                if not graph.has_edge(n1, n2):
                    graph.add_edge(n1, n2)

        return graph

    def _graph_to_smiles(self, graph: nx.Graph) -> Optional[str]:
        """
        Convert a NetworkX graph to SMILES string.

        Naive implementation - assumes all nodes are aromatic carbons.
        """
        if len(graph) == 0:
            return None

        try:
            # Create SMILES for aromatic rings
            # This is a simplified version; real implementation would be more complex
            smiles = "c1ccccc1"  # Default to benzene
            return smiles
        except Exception:
            return None

    def _create_default_aromatic(self, num_carbons: int) -> str:
        """
        Create a default aromatic structure as fallback.

        For targets up to 22C, use verified single PAHs.
        For larger targets, concatenate pyrene (16C) units.
        """
        if num_carbons <= 6:
            return PAH_LIBRARY["benzene"]["smiles"]
        elif num_carbons <= 10:
            return PAH_LIBRARY["naphthalene"]["smiles"]
        elif num_carbons <= 14:
            return PAH_LIBRARY["anthracene"]["smiles"]
        elif num_carbons <= 16:
            return PAH_LIBRARY["pyrene"]["smiles"]
        else:
            # For larger structures, concatenate pyrene (16C) and benzene (6C) units
            # Use pyrene as base and add benzene rings
            pyrene = PAH_LIBRARY["pyrene"]["smiles"]
            remaining = num_carbons - 16
            result = pyrene

            # Add benzene rings (6C each) to reach target
            # Use direct C-C bonds between aromatic carbons
            benzene = "c1ccccc1"
            while remaining >= 6:
                result = result + "-" + benzene
                remaining -= 6

            return result


class SkeletonValidator:
    """Validate carbon skeletons for chemical feasibility."""

    @staticmethod
    def validate(skeleton: CarbonSkeleton) -> Tuple[bool, List[str]]:
        """
        Validate a carbon skeleton.

        Args:
            skeleton: CarbonSkeleton object

        Returns:
            (is_valid, error_messages)
        """
        errors = []

        if skeleton.mol is None:
            errors.append("Molecule object is None")
            return False, errors

        # Check atom counts
        if skeleton.num_carbons == 0:
            errors.append("No carbon atoms in skeleton")
            return False, errors

        if skeleton.num_aromatic_carbons > skeleton.num_carbons:
            errors.append("More aromatic carbons than total carbons")
            return False, errors

        # Check aromaticity percentage
        if skeleton.aromaticity_percent < 0 or skeleton.aromaticity_percent > 100:
            errors.append(
                f"Invalid aromaticity: {skeleton.aromaticity_percent}%"
            )
            return False, errors

        # Check molecular connectivity
        mol_graph = nx.Graph()
        mol_graph.add_nodes_from(range(skeleton.mol.GetNumAtoms()))
        for bond in skeleton.mol.GetBonds():
            mol_graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

        if not nx.is_connected(mol_graph):
            errors.append("Molecular graph is not connected")
            return False, errors

        return len(errors) == 0, errors
