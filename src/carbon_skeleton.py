"""
Carbon Skeleton Generation

Generates aromatic carbon frameworks for biochar structures using PAHs or random graphs.
"""

import random
from typing import List, Tuple, Set, Optional
from dataclasses import dataclass

import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, RWMol

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
        try:
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

        except Exception:
            return self._fallback_to_random(target_num_carbons, target_aromaticity)

    def _select_pah_smiles(self, target_num_carbons: int, prefer_larger: bool) -> str:
        """
        Select appropriate PAH to reach closest match to target carbon count.

        Uses verified PAH SMILES up to pyrene (16C). For larger targets,
        raises exception to trigger RandomGraphGenerator fallback.
        """
        if target_num_carbons <= 6:
            return PAH_LIBRARY["benzene"]["smiles"]
        elif target_num_carbons <= 10:
            return PAH_LIBRARY["naphthalene"]["smiles"]
        elif target_num_carbons <= 14:
            return PAH_LIBRARY["anthracene"]["smiles"]
        elif target_num_carbons <= 16:
            return PAH_LIBRARY["pyrene"]["smiles"]
        else:
            raise ValueError(
                f"Target {target_num_carbons}C exceeds maximum single PAH. "
                "Falling back to modular PAH assembly."
            )

    def _fuse_pahs(self, mol1: Chem.Mol, mol2: Chem.Mol) -> Optional[Chem.Mol]:
        """Attempt to fuse two PAH molecules at aromatic edges."""
        try:
            smiles1 = Chem.MolToSmiles(mol1)
            smiles2 = Chem.MolToSmiles(mol2)
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

    For targets ≤16C: uses verified PAH templates.
    For larger targets: builds 2D graphene nanoflakes (fused aromatic sheets).
    """

    def __init__(self, seed: int = None):
        self.seed = seed
        if seed is not None:
            random.seed(seed)

    def generate(
        self, target_num_carbons: int, target_aromaticity: float = 100.0
    ) -> CarbonSkeleton:
        """
        Generate an aromatic carbon skeleton.

        For ≤16C: uses PAH templates.
        For >16C: builds 2D graphene nanoflakes via ring fusion.

        Args:
            target_num_carbons: Target number of carbon atoms
            target_aromaticity: Target aromaticity percentage

        Returns:
            CarbonSkeleton object
        """
        mol = None
        smiles = None

        if target_num_carbons <= 16:
            smiles = self._create_pah_template(target_num_carbons)
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                Chem.SanitizeMol(mol)
        else:
            # Build 2D graphene nanoflake
            try:
                mol = self._build_graphene_flake_mol(target_num_carbons)
            except Exception:
                mol = None

            if mol is not None:
                try:
                    smiles = Chem.MolToSmiles(mol)
                except Exception:
                    # If SMILES generation fails, it's okay - we'll use the mol directly
                    smiles = None

        # Fallback to polyphenyl chain if nanoflake builder failed
        if mol is None:
            smiles = self._build_polyphenyl(max(3, round(target_num_carbons / 6)))
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                Chem.SanitizeMol(mol)

        # Final fallback: benzene
        if mol is None:
            smiles = "c1ccccc1"
            mol = Chem.MolFromSmiles(smiles)
            Chem.SanitizeMol(mol)

        num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        num_aromatic = sum(
            1 for atom in mol.GetAtoms()
            if atom.GetAtomicNum() == 6 and atom.GetIsAromatic()
        )
        aromaticity = (num_aromatic / num_carbons * 100) if num_carbons > 0 else 0

        return CarbonSkeleton(
            mol=mol,
            smiles=smiles or Chem.MolToSmiles(mol),
            num_carbons=num_carbons,
            num_aromatic_carbons=num_aromatic,
            aromaticity_percent=aromaticity,
        )

    # ------------------------------------------------------------------
    # 2D graphene nanoflake builder
    # ------------------------------------------------------------------

    def _build_graphene_flake_mol(self, target_carbons: int) -> Optional[Chem.Mol]:
        """
        Build a 2D graphene nanoflake by iteratively fusing benzene rings.

        Strategy:
        1. Start from benzene.
        2. Repeatedly fuse a new 6-membered ring at the most "interior"
           available edge (highest sum of neighbor degrees) to grow a compact
           2D sheet rather than a linear chain.
        3. Convert to RDKit molecule using a Kekulé bond assignment found via
           networkx maximum matching, then let SanitizeMol perceive aromaticity.

        Each new fused ring adds 4 carbons (2 shared with existing ring).

        Args:
            target_carbons: Approximate target carbon count

        Returns:
            Sanitized RDKit Mol, or None if sanitization fails.
        """
        # Build the connectivity graph
        G = nx.Graph()

        # Initial benzene ring: nodes 0–5
        for i in range(6):
            G.add_node(i)
        for i in range(6):
            G.add_edge(i, (i + 1) % 6)

        next_node = 6

        while G.number_of_nodes() + 4 <= target_carbons:
            fusable = self._get_fusable_edges(G)
            if not fusable:
                break

            # Prefer edges whose non-shared neighbors have higher degree
            # → promotes compact 2D growth over linear acene chains
            fusable.sort(key=lambda e: -self._edge_neighbor_degree(G, e[0], e[1]))
            u, v = fusable[0]

            # Add 4 new nodes and wire them into a new 6-membered ring
            new_nodes = list(range(next_node, next_node + 4))
            next_node += 4
            for n in new_nodes:
                G.add_node(n)
            # Ring: u - v - n0 - n1 - n2 - n3 - u
            G.add_edge(v, new_nodes[0])
            G.add_edge(new_nodes[0], new_nodes[1])
            G.add_edge(new_nodes[1], new_nodes[2])
            G.add_edge(new_nodes[2], new_nodes[3])
            G.add_edge(new_nodes[3], u)

        # Assign Kekulé bonds via maximum matching
        # In a valid PAH, every carbon in the π system participates in one
        # double bond (perfect matching exists for even-carbon systems).
        matching = nx.max_weight_matching(G, maxcardinality=True)
        double_bonds: Set[Tuple[int, int]] = {
            (min(u, v), max(u, v)) for u, v in matching
        }

        # Build RDKit molecule with explicit SINGLE/DOUBLE bonds
        node_list = sorted(G.nodes())
        node_to_idx = {n: i for i, n in enumerate(node_list)}

        rwmol = RWMol()
        for _ in node_list:
            rwmol.AddAtom(Chem.Atom(6))

        for u, v in G.edges():
            iu, iv = node_to_idx[u], node_to_idx[v]
            key = (min(u, v), max(u, v))
            btype = (
                Chem.BondType.DOUBLE if key in double_bonds else Chem.BondType.SINGLE
            )
            rwmol.AddBond(iu, iv, btype)

        mol = rwmol.GetMol()
        try:
            Chem.SanitizeMol(mol)
            # Re-perceive aromaticity to ensure all aromatic systems are recognized
            Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
            return mol
        except Exception:
            # If Kekulé assignment failed (e.g. imperfect matching for odd-carbon
            # intermediate), try with aromatic bonds as a last resort
            return self._build_graphene_flake_aromatic(G)

    def _build_graphene_flake_aromatic(self, G: nx.Graph) -> Optional[Chem.Mol]:
        """Fallback: build with AROMATIC bonds and let RDKit sanitize."""
        node_list = sorted(G.nodes())
        node_to_idx = {n: i for i, n in enumerate(node_list)}

        rwmol = RWMol()
        for _ in node_list:
            atom = Chem.Atom(6)
            atom.SetNoImplicit(False)
            rwmol.AddAtom(atom)

        for u, v in G.edges():
            iu, iv = node_to_idx[u], node_to_idx[v]
            rwmol.AddBond(iu, iv, Chem.BondType.AROMATIC)

        mol = rwmol.GetMol()
        try:
            Chem.SanitizeMol(mol)
            # Re-perceive aromaticity to ensure all aromatic systems are recognized
            Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
            return mol
        except Exception:
            return None

    @staticmethod
    def _get_fusable_edges(G: nx.Graph) -> List[Tuple[int, int]]:
        """Return edges where both endpoints have degree ≤ 2 (open for ring fusion)."""
        return [
            (u, v)
            for u, v in G.edges()
            if G.degree(u) == 2 and G.degree(v) == 2
        ]

    @staticmethod
    def _edge_neighbor_degree(G: nx.Graph, a: int, b: int) -> int:
        """Sum of degrees of neighbors of a and b (excluding each other)."""
        neighbors_a = {n for n in G.neighbors(a) if n != b}
        neighbors_b = {n for n in G.neighbors(b) if n != a}
        return sum(G.degree(n) for n in neighbors_a | neighbors_b)

    # ------------------------------------------------------------------
    # Template and SMILES helpers
    # ------------------------------------------------------------------

    def _create_pah_template(self, num_carbons: int) -> str:
        """Return a verified PAH SMILES for the given carbon count."""
        if num_carbons <= 6:
            return PAH_LIBRARY["benzene"]["smiles"]
        elif num_carbons <= 10:
            return PAH_LIBRARY["naphthalene"]["smiles"]
        elif num_carbons <= 14:
            return PAH_LIBRARY["anthracene"]["smiles"]
        else:
            return PAH_LIBRARY["pyrene"]["smiles"]

    def _build_polyphenyl(self, n_rings: int) -> str:
        """
        Build a para-polyphenyl chain SMILES with n_rings benzene rings.
        Kept as ultimate fallback for cases where the nanoflake builder fails.
        """
        smiles = "c1ccccc1"
        for _ in range(n_rings - 1):
            smiles = f"c1ccc(-{smiles})cc1"
        return smiles

    # Legacy stubs kept for API compatibility
    def _create_random_aromatic_graph(self, num_nodes: int, target_aromaticity: float) -> nx.Graph:
        graph = nx.Graph()
        num_rings = max(1, num_nodes // 6)
        for i in range(num_rings):
            ring_nodes = list(range(i * 6, (i + 1) * 6))
            graph.add_nodes_from(ring_nodes)
            for j in range(6):
                graph.add_edge(ring_nodes[j], ring_nodes[(j + 1) % 6])
        all_nodes = list(graph.nodes())
        for _ in range(max(0, num_rings - 1)):
            if len(all_nodes) >= 2:
                n1, n2 = random.sample(all_nodes, 2)
                if not graph.has_edge(n1, n2):
                    graph.add_edge(n1, n2)
        return graph

    def _graph_to_smiles(self, graph: nx.Graph) -> Optional[str]:
        return None


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

        if skeleton.num_carbons == 0:
            errors.append("No carbon atoms in skeleton")
            return False, errors

        if skeleton.num_aromatic_carbons > skeleton.num_carbons:
            errors.append("More aromatic carbons than total carbons")
            return False, errors

        if skeleton.aromaticity_percent < 0 or skeleton.aromaticity_percent > 100:
            errors.append(f"Invalid aromaticity: {skeleton.aromaticity_percent}%")
            return False, errors

        mol_graph = nx.Graph()
        mol_graph.add_nodes_from(range(skeleton.mol.GetNumAtoms()))
        for bond in skeleton.mol.GetBonds():
            mol_graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

        if not nx.is_connected(mol_graph):
            errors.append("Molecular graph is not connected")
            return False, errors

        return len(errors) == 0, errors
