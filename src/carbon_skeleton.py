"""
Carbon Skeleton Generation

Generates aromatic carbon frameworks for biochar structures using PAHs or random graphs.
"""

import math
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


# ---------------------------------------------------------------------------
# Utility: convert RDKit Mol (carbon-only) to NetworkX graph
# ---------------------------------------------------------------------------

def _mol_to_graph(mol: Chem.Mol) -> nx.Graph:
    """Convert an RDKit carbon-only Mol to a NetworkX graph."""
    G = nx.Graph()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            G.add_node(atom.GetIdx())
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i in G.nodes and j in G.nodes:
            G.add_edge(i, j)
    return G


def _graph_to_mol(G: nx.Graph) -> Optional[Chem.Mol]:
    """
    Convert a NetworkX carbon graph to an RDKit Mol.

    Uses aromatic bond type so RDKit handles kekulization internally,
    which correctly perceives all fused-ring aromatic carbons.
    Falls back to max-matching Kekule assignment if sanitization fails.
    """
    if G.number_of_nodes() == 0:
        return None

    node_list = sorted(G.nodes())
    node_to_idx = {n: i for i, n in enumerate(node_list)}

    # Primary: all-aromatic bonds, let SanitizeMol kekulize
    rwmol = RWMol()
    for _ in node_list:
        a = Chem.Atom(6)
        a.SetNoImplicit(False)
        rwmol.AddAtom(a)
    for u, v in G.edges():
        rwmol.AddBond(node_to_idx[u], node_to_idx[v], Chem.BondType.AROMATIC)
    mol = rwmol.GetMol()
    try:
        Chem.SanitizeMol(mol)
        Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
        return mol
    except Exception:
        pass

    # Fallback: explicit Kekule via max matching (requires even node count)
    if G.number_of_nodes() % 2 != 0:
        return None

    matching = nx.max_weight_matching(G, maxcardinality=True)
    double_bonds: Set[Tuple[int, int]] = {
        (min(u, v), max(u, v)) for u, v in matching
    }
    rwmol2 = RWMol()
    for _ in node_list:
        rwmol2.AddAtom(Chem.Atom(6))
    for u, v in G.edges():
        iu, iv = node_to_idx[u], node_to_idx[v]
        key = (min(u, v), max(u, v))
        btype = Chem.BondType.DOUBLE if key in double_bonds else Chem.BondType.SINGLE
        rwmol2.AddBond(iu, iv, btype)
    mol2 = rwmol2.GetMol()
    try:
        Chem.SanitizeMol(mol2)
        Chem.SetAromaticity(mol2, Chem.AromaticityModel.AROMATICITY_MDL)
        return mol2
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Hex-lattice seed builder
# ---------------------------------------------------------------------------

def _build_hex_lattice_graph(num_rings: int) -> nx.Graph:
    """
    Build a compact graphene-like PAH graph from a hexagonal lattice.

    Uses axial hex coordinates. Starts with the centre cell and spirals
    outward.  Each hex cell is a 6-membered ring; shared edges between
    adjacent cells are automatically deduplicated.

    The resulting graph is always bipartite (A/B sublattice of the
    honeycomb), which guarantees a perfect matching exists (by König's
    theorem) when the graph has an even number of nodes.

    Args:
        num_rings: Number of hexagonal cells to place (1 = benzene, 7 = coronene)

    Returns:
        NetworkX graph of the carbon backbone.
    """
    # Generate hex cell centres in axial coordinates using a spiral
    cells = _hex_spiral(num_rings)

    # For each cell, compute its 6 corner positions (using "pointy-top" hex).
    # Corners are shared between adjacent cells.
    # We identify corners by rounding their (x, y) to avoid float duplicates.
    corner_id_map: dict = {}  # (rounded_x, rounded_y) -> node_id
    cell_corners: list = []   # list of lists of node_ids per cell
    next_id = 0

    for q, r in cells:
        # Centre of hex cell in Cartesian (pointy-top orientation)
        cx = math.sqrt(3) * (q + r / 2.0)
        cy = 1.5 * r

        corners = []
        for k in range(6):
            angle = math.pi / 3.0 * k + math.pi / 6.0
            px = cx + math.cos(angle)
            py = cy + math.sin(angle)
            key = (round(px, 4), round(py, 4))
            if key not in corner_id_map:
                corner_id_map[key] = next_id
                next_id += 1
            corners.append(corner_id_map[key])
        cell_corners.append(corners)

    # Build the graph: nodes are corner IDs, edges are ring edges
    G = nx.Graph()
    G.add_nodes_from(range(next_id))
    for corners in cell_corners:
        for i in range(6):
            G.add_edge(corners[i], corners[(i + 1) % 6])

    return G


def _hex_spiral(n: int) -> List[Tuple[int, int]]:
    """Return *n* hex-cell axial coordinates in a compact spiral.

    Cells are emitted from the centre outwards, ring by ring.  Each ring
    is traversed using the correct ring-edge direction vectors so that
    consecutive cells on the same ring are always adjacent (share an
    edge), producing a compact 2-D cluster rather than a strip.
    """
    if n <= 0:
        return []
    cells = [(0, 0)]
    if n == 1:
        return cells

    # Direction vectors for traversing the *perimeter* of ring r.
    # Starting at (r, 0) and going counter-clockwise in axial coords.
    ring_dirs = [(-1, 1), (-1, 0), (0, -1), (1, -1), (1, 0), (0, 1)]
    radius = 1
    while len(cells) < n:
        q, r = radius, 0  # Start of ring `radius`
        for d in range(6):
            for _ in range(radius):
                if len(cells) >= n:
                    return cells
                cells.append((q, r))
                q += ring_dirs[d][0]
                r += ring_dirs[d][1]
        radius += 1
    return cells[:n]


# ---------------------------------------------------------------------------
# Ring-growth engine (shared by PAHAssembler and RandomGraphGenerator)
# ---------------------------------------------------------------------------

def _grow_graph(G: nx.Graph, target_nodes: int, seed: Optional[int] = None) -> nx.Graph:
    """
    Grow a carbon graph by iteratively fusing 6-membered rings.

    Each new ring shares one edge with the existing graph and adds 4
    new carbon nodes.  The algorithm prefers edges whose neighbours have
    higher total degree, promoting compact 2D growth.

    A parity guard ensures the graph always has an even number of nodes
    (required for Kekule / perfect-matching assignment).

    Args:
        G: Seed graph (will NOT be mutated; a copy is made).
        target_nodes: Desired number of nodes.
        seed: Optional RNG seed for tie-breaking.

    Returns:
        Grown graph (new object).
    """
    G = G.copy()
    rng = random.Random(seed)
    next_node = max(G.nodes()) + 1 if G.nodes() else 0

    while G.number_of_nodes() + 4 <= target_nodes:
        fusable = _get_fusable_edges(G)
        if not fusable:
            break

        # Sort by compactness heuristic; break ties randomly
        fusable.sort(key=lambda e: (-_edge_neighbor_degree(G, e[0], e[1]), rng.random()))
        u, v = fusable[0]

        # Add 4 new nodes forming a new 6-membered ring sharing edge u-v
        new_nodes = list(range(next_node, next_node + 4))
        next_node += 4
        for n in new_nodes:
            G.add_node(n)
        G.add_edge(v, new_nodes[0])
        G.add_edge(new_nodes[0], new_nodes[1])
        G.add_edge(new_nodes[1], new_nodes[2])
        G.add_edge(new_nodes[2], new_nodes[3])
        G.add_edge(new_nodes[3], u)

    # Parity guard: if odd node count, try adding one more ring
    if G.number_of_nodes() % 2 != 0:
        fusable = _get_fusable_edges(G)
        if fusable:
            fusable.sort(key=lambda e: (-_edge_neighbor_degree(G, e[0], e[1]), rng.random()))
            u, v = fusable[0]
            new_nodes = list(range(next_node, next_node + 4))
            for n in new_nodes:
                G.add_node(n)
            G.add_edge(v, new_nodes[0])
            G.add_edge(new_nodes[0], new_nodes[1])
            G.add_edge(new_nodes[1], new_nodes[2])
            G.add_edge(new_nodes[2], new_nodes[3])
            G.add_edge(new_nodes[3], u)

    return G


def _get_fusable_edges(G: nx.Graph) -> List[Tuple[int, int]]:
    """Return edges where both endpoints have degree <= 2 (open for ring fusion)."""
    return [(u, v) for u, v in G.edges() if G.degree(u) == 2 and G.degree(v) == 2]


def _edge_neighbor_degree(G: nx.Graph, a: int, b: int) -> int:
    """Sum of degrees of neighbours of a and b (excluding each other)."""
    na = {n for n in G.neighbors(a) if n != b}
    nb = {n for n in G.neighbors(b) if n != a}
    return sum(G.degree(n) for n in na | nb)


# ---------------------------------------------------------------------------
# PAH Assembler (main public API)
# ---------------------------------------------------------------------------

class PAHAssembler:
    """
    Build aromatic carbon skeletons from validated PAH building blocks.

    Strategy:
      <= 24C  : use a pre-validated PAH from the library (exact or nearest).
      > 24C   : build a hex-lattice seed, then grow via ring fusion.
    """

    # Sorted list of (num_carbons, pah_name) for quick lookup
    _SIZE_INDEX: Optional[List[Tuple[int, str]]] = None

    def __init__(self, seed: int = None):
        self.seed = seed
        if seed is not None:
            random.seed(seed)

        # Preload PAH library
        self.pahs = {}
        for name, data in PAH_LIBRARY.items():
            mol = Chem.MolFromSmiles(data["smiles"])
            if mol is not None:
                try:
                    Chem.SanitizeMol(mol)
                    nC = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
                    self.pahs[name] = {"mol": mol, "num_carbons": nC, "data": data}
                except Exception:
                    pass

        # Build size index (sorted ascending by carbon count)
        if PAHAssembler._SIZE_INDEX is None:
            PAHAssembler._SIZE_INDEX = sorted(
                [(info["num_carbons"], name) for name, info in self.pahs.items()]
            )

    def generate(
        self,
        target_num_carbons: int,
        target_aromaticity: float = 100.0,
        prefer_larger_pahs: bool = True,
    ) -> CarbonSkeleton:
        """Generate a PAH carbon skeleton of approximately *target_num_carbons*."""
        mol = None

        # --- Route 1: exact library match ---
        for name, info in self.pahs.items():
            if info["num_carbons"] == target_num_carbons:
                mol = Chem.Mol(info["mol"])
                break

        # --- Route 2: library seed + graph growth ---
        if mol is None:
            mol = self._build_from_seed(target_num_carbons)

        # --- Fallback: pyrene ---
        if mol is None:
            mol = Chem.Mol(self.pahs["pyrene"]["mol"])

        Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
        return self._make_skeleton(mol)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_from_seed(self, target: int) -> Optional[Chem.Mol]:
        """
        Build a PAH of *target* carbons using a library seed plus graph growth.

        Seed selection uses parity: each ring addition in _grow_graph adds
        exactly 4 nodes, so we need (target - seed_size) % 4 == 0.  Starting
        from an even seed and adding 4 at a time always produces an even node
        count, which is required for a valid Kekule structure.
        """
        # Pick the largest library PAH where (target - nC) % 4 == 0 and nC <= target
        seed_mol = None
        seed_carbons = 0
        for nC, name in reversed(self._SIZE_INDEX):
            if nC <= target and (target - nC) % 4 == 0 and name in self.pahs:
                seed_mol = self.pahs[name]["mol"]
                seed_carbons = nC
                break

        if seed_mol is None:
            # No library match with correct parity; round target up to nearest valid size
            # from benzene (6C): valid targets are 6 + 4k
            rem = (target - 6) % 4
            if rem != 0:
                target = target + (4 - rem)
            seed_mol = Chem.MolFromSmiles(PAH_LIBRARY["benzene"]["smiles"])
            seed_carbons = 6

        if seed_carbons >= target:
            return Chem.Mol(seed_mol)

        G = _mol_to_graph(seed_mol)
        G = _grow_graph(G, target, seed=self.seed)
        return _graph_to_mol(G)

    @staticmethod
    def _make_skeleton(mol: Chem.Mol) -> CarbonSkeleton:
        num_carbons = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
        num_aromatic = sum(
            1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetIsAromatic()
        )
        aromaticity = (num_aromatic / num_carbons * 100) if num_carbons > 0 else 0
        try:
            smiles = Chem.MolToSmiles(mol)
        except Exception:
            smiles = f"C{num_carbons}_PAH"
        return CarbonSkeleton(
            mol=mol,
            smiles=smiles,
            num_carbons=num_carbons,
            num_aromatic_carbons=num_aromatic,
            aromaticity_percent=aromaticity,
        )


# ---------------------------------------------------------------------------
# RandomGraphGenerator (kept for backward compatibility)
# ---------------------------------------------------------------------------

class RandomGraphGenerator:
    """
    Generate random aromatic carbon frameworks.

    Delegates to hex-lattice builder + ring growth for all sizes.
    """

    def __init__(self, seed: int = None):
        self.seed = seed
        if seed is not None:
            random.seed(seed)

    def generate(
        self, target_num_carbons: int, target_aromaticity: float = 100.0
    ) -> CarbonSkeleton:
        # Delegate to PAHAssembler for consistent parity-aware seed selection
        asm = PAHAssembler(seed=self.seed)
        return asm.generate(target_num_carbons, target_aromaticity)

    def _create_pah_template(self, num_carbons: int) -> str:
        if num_carbons <= 6:
            return PAH_LIBRARY["benzene"]["smiles"]
        elif num_carbons <= 10:
            return PAH_LIBRARY["naphthalene"]["smiles"]
        elif num_carbons <= 14:
            return PAH_LIBRARY["anthracene"]["smiles"]
        else:
            return PAH_LIBRARY["pyrene"]["smiles"]


# ---------------------------------------------------------------------------
# Skeleton Validator
# ---------------------------------------------------------------------------

class SkeletonValidator:
    """Validate carbon skeletons for chemical feasibility."""

    @staticmethod
    def validate(skeleton: CarbonSkeleton) -> Tuple[bool, List[str]]:
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
