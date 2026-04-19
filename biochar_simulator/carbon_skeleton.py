"""
Carbon Skeleton Generation

Generates aromatic carbon frameworks for biochar structures using PAHs or random graphs.
"""

import math
import random
from typing import List, Tuple, Set, Optional, Dict
from dataclasses import dataclass, field

import networkx as nx
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, RWMol

from .constants import PAH_LIBRARY


# Target aromatic C–C bond length (Å) for hex-lattice layout.  Used to
# rescale seed 2D coords from RDKit (which returns ~1.5 Å bonds) and as the
# edge length for all ring-fusion geometry.
_AROMATIC_CC_BOND = 1.42

# When position-tracking is active, two computed ring-vertex positions that
# lie within this distance (Å) of each other are treated as the same atom.
# This handles shared corners between adjacent hexagons in a fused lattice:
# both growth fronts independently predict the same Cartesian coordinate for
# what is physically one carbon atom.
_POSITION_MERGE_TOLERANCE = 0.05  # Å


@dataclass
class CarbonSkeleton:
    """Represents a carbon skeleton/framework."""

    mol: Chem.Mol
    smiles: str
    num_carbons: int
    num_aromatic_carbons: int
    aromaticity_percent: float


# ---------------------------------------------------------------------------
# Ring-fusion geometry helpers (propagate hex-lattice positions during growth)
# ---------------------------------------------------------------------------

def _compute_fused_ring_positions(
    pu: Tuple[float, float],
    pv: Tuple[float, float],
    outward: Tuple[float, float],
    n_sides: int,
    num_new: int,
) -> List[Tuple[float, float]]:
    """
    Compute Cartesian positions for *num_new* new vertices of a regular
    *n_sides*-gon sharing edge (pu, pv), placed on the *outward* side.

    Returns positions in the order matching ``_fuse_hexagon`` /
    ``_fuse_pentagon``:  v → new[0] → new[1] → … → new[num_new-1] → u.
    """
    pu_a = np.array(pu, dtype=float)
    pv_a = np.array(pv, dtype=float)
    edge = pv_a - pu_a
    L = float(np.linalg.norm(edge))
    if L < 1e-8:
        # Degenerate; fall back to a small lattice step along +y
        return [(pu_a[0], pu_a[1] + (i + 1) * _AROMATIC_CC_BOND) for i in range(num_new)]

    e_hat = edge / L
    # Perpendicular (right-handed in 2D)
    n_right = np.array([-e_hat[1], e_hat[0]])
    out = np.array(outward, dtype=float)
    if np.linalg.norm(out) < 1e-8:
        out = n_right
    n_hat = n_right if float(np.dot(n_right, out)) >= 0.0 else -n_right

    apothem = L / (2.0 * math.tan(math.pi / n_sides))
    R = L / (2.0 * math.sin(math.pi / n_sides))
    midpoint = (pu_a + pv_a) / 2.0
    center = midpoint + n_hat * apothem

    theta_v = math.atan2(pv_a[1] - center[1], pv_a[0] - center[0])
    # Choose angular direction (CCW or CW) so traversal wraps AWAY from the
    # shared edge (through the new vertices) rather than across it to u.
    cross = e_hat[0] * n_hat[1] - e_hat[1] * n_hat[0]
    direction = 1.0 if cross > 0 else -1.0

    step = 2.0 * math.pi / n_sides
    positions = []
    for k in range(num_new):
        theta = theta_v + direction * (k + 1) * step
        px = center[0] + R * math.cos(theta)
        py = center[1] + R * math.sin(theta)
        positions.append((float(px), float(py)))
    return positions


def _outward_direction(
    G: nx.Graph,
    positions: Dict[int, Tuple[float, float]],
    u: int,
    v: int,
) -> Tuple[float, float]:
    """
    Estimate the outward normal direction at fusable edge (u, v).

    A fusable edge has both endpoints at degree 2.  Each endpoint has
    exactly one "other neighbour" in the existing graph; those neighbours
    sit on the interior side.  The outward direction points from their
    centroid toward the edge midpoint.
    """
    mid = ((positions[u][0] + positions[v][0]) / 2.0,
           (positions[u][1] + positions[v][1]) / 2.0)
    interior = []
    for nbr in G.neighbors(u):
        if nbr != v and nbr in positions:
            interior.append(positions[nbr])
    for nbr in G.neighbors(v):
        if nbr != u and nbr in positions:
            interior.append(positions[nbr])
    if not interior:
        return (0.0, 1.0)
    cx = sum(p[0] for p in interior) / len(interior)
    cy = sum(p[1] for p in interior) / len(interior)
    return (mid[0] - cx, mid[1] - cy)


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


def _graph_to_mol(
    G: nx.Graph,
    positions: Optional[Dict[int, Tuple[float, float]]] = None,
) -> Optional[Chem.Mol]:
    """
    Convert a NetworkX carbon graph to an RDKit Mol.

    Uses aromatic bond type so RDKit handles kekulization internally,
    which correctly perceives all fused-ring aromatic carbons.
    Falls back to max-matching Kekule assignment if sanitization fails.

    If *positions* is provided, each carbon atom is tagged with double
    properties ``init_x`` / ``init_y`` so downstream geometry generation
    can use the hex-lattice layout directly (bypassing Compute2DCoords /
    kamada-kawai on large fused graphs).
    """
    if G.number_of_nodes() == 0:
        return None

    node_list = sorted(G.nodes())
    node_to_idx = {n: i for i, n in enumerate(node_list)}

    def _attach_positions(m: Chem.Mol) -> None:
        if positions is None:
            return
        for orig_node, idx in node_to_idx.items():
            if orig_node in positions:
                x, y = positions[orig_node]
                atom = m.GetAtomWithIdx(idx)
                atom.SetDoubleProp("init_x", float(x))
                atom.SetDoubleProp("init_y", float(y))

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
        _attach_positions(mol)
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
        _attach_positions(mol2)
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

def _find_node_at_position(
    positions: Dict[int, Tuple[float, float]],
    target: Tuple[float, float],
    tolerance: float = _POSITION_MERGE_TOLERANCE,
) -> Optional[int]:
    """
    Return the id of a node already at *target* (within *tolerance* Å),
    or ``None`` if no such node exists.

    Called during ring fusion to detect shared corners in the hex lattice:
    when a newly-computed vertex position coincides with an existing node the
    two should be the same atom, not independent duplicates.
    """
    tx, ty = target
    for node_id, (px, py) in positions.items():
        if math.hypot(px - tx, py - ty) < tolerance:
            return node_id
    return None


def _fuse_hexagon(
    G: nx.Graph,
    u: int,
    v: int,
    next_node: int,
    positions: Optional[Dict[int, Tuple[float, float]]] = None,
) -> int:
    """
    Fuse a 6-membered ring onto edge (u, v).

    When *positions* is provided the four ring vertices are located
    analytically **before** any graph modification.  If a computed position
    coincides with an existing node (within ``_POSITION_MERGE_TOLERANCE``)
    that node is **reused** rather than duplicated — correctly modelling the
    shared-corner topology of a graphene lattice when two growth fronts meet.

    Returns the next free node id.
    """
    if positions is not None and u in positions and v in positions:
        # 1. Compute the four vertex positions of the new ring FIRST.
        outward = _outward_direction(G, positions, u, v)
        new_pos_list = _compute_fused_ring_positions(
            positions[u], positions[v], outward, n_sides=6, num_new=4
        )

        # 2. For each vertex: reuse an existing node if the position is
        #    already occupied, otherwise allocate a fresh node.
        new_nodes: List[int] = []
        for pos in new_pos_list:
            existing = _find_node_at_position(positions, pos)
            if existing is not None:
                new_nodes.append(existing)
            else:
                G.add_node(next_node)
                positions[next_node] = pos
                new_nodes.append(next_node)
                next_node += 1

        # 3. Add the ring edges (nx.Graph silently ignores duplicate edges).
        ring_seq = [v, new_nodes[0], new_nodes[1], new_nodes[2], new_nodes[3], u]
        for a, b in zip(ring_seq, ring_seq[1:]):
            G.add_edge(a, b)

        return next_node

    else:
        # No position tracking: original behaviour (always allocate 4 new nodes).
        new_nodes = list(range(next_node, next_node + 4))
        for n in new_nodes:
            G.add_node(n)
        G.add_edge(v, new_nodes[0])
        G.add_edge(new_nodes[0], new_nodes[1])
        G.add_edge(new_nodes[1], new_nodes[2])
        G.add_edge(new_nodes[2], new_nodes[3])
        G.add_edge(new_nodes[3], u)
        return next_node + 4


def _fuse_pentagon(
    G: nx.Graph,
    u: int,
    v: int,
    next_node: int,
    positions: Optional[Dict[int, Tuple[float, float]]] = None,
) -> int:
    """
    Fuse a 5-membered ring onto edge (u, v).

    Applies the same position-merge logic as ``_fuse_hexagon``: computed
    vertices that coincide with existing nodes are reused rather than
    duplicated.

    Returns the next free node id.
    """
    if positions is not None and u in positions and v in positions:
        outward = _outward_direction(G, positions, u, v)
        new_pos_list = _compute_fused_ring_positions(
            positions[u], positions[v], outward, n_sides=5, num_new=3
        )

        new_nodes: List[int] = []
        for pos in new_pos_list:
            existing = _find_node_at_position(positions, pos)
            if existing is not None:
                new_nodes.append(existing)
            else:
                G.add_node(next_node)
                positions[next_node] = pos
                new_nodes.append(next_node)
                next_node += 1

        ring_seq = [v, new_nodes[0], new_nodes[1], new_nodes[2], u]
        for a, b in zip(ring_seq, ring_seq[1:]):
            G.add_edge(a, b)

        return next_node

    else:
        new_nodes = list(range(next_node, next_node + 3))
        for n in new_nodes:
            G.add_node(n)
        G.add_edge(v, new_nodes[0])
        G.add_edge(new_nodes[0], new_nodes[1])
        G.add_edge(new_nodes[1], new_nodes[2])
        G.add_edge(new_nodes[2], u)
        return next_node + 3


def _grow_graph(
    G: nx.Graph,
    target_nodes: int,
    seed: Optional[int] = None,
    defect_fraction: float = 0.0,
    positions: Optional[Dict[int, Tuple[float, float]]] = None,
) -> nx.Graph:
    """
    Grow a carbon graph by iteratively fusing rings.

    By default only 6-membered rings are added (+4 nodes each).  When
    *defect_fraction* > 0 the algorithm occasionally fuses a 5-membered
    ring (+3 nodes) instead, introducing pentagonal defects analogous to
    those observed in disordered graphitic carbon.

    Each new ring shares one edge with the existing graph.  The algorithm
    prefers edges whose neighbours have higher total degree, promoting
    compact 2D growth.

    A parity guard ensures the final graph has an even number of nodes
    (required for a valid Kekulé / perfect-matching assignment):
      - Pure hexagon mode (defect_fraction == 0): adds one extra hexagon.
      - Defect mode (defect_fraction >  0): adds one extra pentagon
        because odd + 3 = even, whereas odd + 4 remains odd.

    Args:
        G: Seed graph (will NOT be mutated; a copy is made).
        target_nodes: Desired number of nodes.
        seed: Optional RNG seed for tie-breaking and pentagon selection.
        defect_fraction: Probability [0, 1) that any given ring addition
            is a 5-membered (pentagon) ring rather than a 6-membered ring.

    Returns:
        Grown graph (new object).
    """
    G = G.copy()
    rng = random.Random(seed)
    next_node = max(G.nodes()) + 1 if G.nodes() else 0

    # When position-tracking is active some ring fusions add zero new nodes
    # (shared corners are merged into existing nodes).  We must therefore
    # continue the loop until the *unique* node count reaches target_nodes,
    # not merely until we've attempted target/min_step iterations.
    #
    # For the common no-merging case the two conditions are equivalent:
    # the seed is chosen so (target - seed) % 4 == 0, each step adds
    # exactly 4 nodes, and we reach target_nodes on the last iteration
    # regardless of whether the guard reads  "< target"  or
    # "+ 4 <= target".

    while G.number_of_nodes() < target_nodes:
        fusable = _get_fusable_edges(G)
        if not fusable:
            break

        # Sort by compactness heuristic; break ties randomly
        fusable.sort(key=lambda e: (-_edge_neighbor_degree(G, e[0], e[1]), rng.random()))
        u, v = fusable[0]

        nodes_remaining = target_nodes - G.number_of_nodes()

        # Decide ring size for this step.
        # Pentagon adds ≤ 3 new nodes; hexagon adds ≤ 4.
        # In defect mode: if exactly 3 slots remain use a pentagon so we
        # land on the target rather than overshooting by 1.
        if defect_fraction > 0 and nodes_remaining == 3:
            use_pentagon = True
        elif defect_fraction > 0 and nodes_remaining >= 4:
            use_pentagon = rng.random() < defect_fraction
        else:
            use_pentagon = False

        if use_pentagon:
            next_node = _fuse_pentagon(G, u, v, next_node, positions=positions)
        else:
            next_node = _fuse_hexagon(G, u, v, next_node, positions=positions)

    # When position-tracking was active, some ring fusions may have created
    # non-adjacent pairs at bond distance (1.42 Å) that were never directly
    # bonded (see ``_add_missing_lattice_bonds`` docstring).  Fix these before
    # building the mol so the graph has the correct topology.  The 1.45 Å
    # tolerance safely distinguishes genuine bonds from all non-bonded C-C
    # pairs (min non-bonded = 2.30 Å in pentagons, 2.46 Å in hexagons).
    if positions is not None:
        _add_missing_lattice_bonds(G, positions)

    # Parity guard: ensure even node count for Kekulé assignment.
    # Run as a bounded loop: when position-merging is active a single ring
    # fusion may add zero new nodes (all corners merged), leaving the count
    # still odd — we must try additional edges until parity is restored or
    # no more fusable edges remain.
    _parity_limit = max(10, G.number_of_nodes())
    for _ in range(_parity_limit):
        if G.number_of_nodes() % 2 == 0:
            break
        fusable = _get_fusable_edges(G)
        if not fusable:
            break
        fusable.sort(key=lambda e: (-_edge_neighbor_degree(G, e[0], e[1]), rng.random()))
        u, v = fusable[0]
        if defect_fraction > 0:
            # Pentagon adds odd → odd + odd = even  ✓
            next_node = _fuse_pentagon(G, u, v, next_node, positions=positions)
        else:
            # Hexagon adds even → pure-hexagon safety net
            next_node = _fuse_hexagon(G, u, v, next_node, positions=positions)

    return G


def _add_missing_lattice_bonds(
    G: nx.Graph,
    positions: Dict[int, Tuple[float, float]],
    tolerance: float = 1.45,
) -> None:
    """
    Add any bond-distance C-C bonds that are missing from the graph.

    This can occur when position-merging places two existing nodes at adjacent
    hex-lattice sites (1.42 Å) as non-adjacent members of the same fused ring
    (e.g. positions new[0] and new[3]).  The ring-sequence edges connect them
    through intermediates but not directly; the direct edge must be added here.

    The tolerance (1.45 Å < 2.46 Å = min non-bonded distance in graphene)
    guarantees that only genuine lattice bonds are added, never spurious ones.

    Operates **in place** on G.
    """
    pos_list = [(n, positions[n]) for n in G.nodes() if n in positions]
    for i, (na, (ax, ay)) in enumerate(pos_list):
        for nb, (bx, by) in pos_list[i + 1:]:
            if G.has_edge(na, nb):
                continue
            if math.hypot(ax - bx, ay - by) < tolerance:
                G.add_edge(na, nb)


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

    @staticmethod
    def _compute_seed_positions(
        seed_mol: Chem.Mol,
    ) -> Optional[Dict[int, Tuple[float, float]]]:
        """
        Generate 2D positions for a library PAH seed and rescale so the
        mean C–C bond length equals ``_AROMATIC_CC_BOND``.

        Works reliably on small library PAHs (benzene, naphthalene,
        pyrene, coronene, …) where RDKit's Compute2DCoords produces a
        clean hexagonal layout.  Returns ``None`` if the layout can't be
        computed.
        """
        try:
            mol = Chem.Mol(seed_mol)
            AllChem.Compute2DCoords(mol)
            conf = mol.GetConformer()
            positions: Dict[int, Tuple[float, float]] = {}
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() != 6:
                    continue
                p = conf.GetAtomPosition(atom.GetIdx())
                positions[atom.GetIdx()] = (float(p.x), float(p.y))
            if not positions:
                return None

            # Rescale so mean bond length ≈ _AROMATIC_CC_BOND (RDKit's 2D
            # coords come out with ~1.5 Å bonds; we want 1.42 Å).
            bond_lengths = []
            for bond in mol.GetBonds():
                a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                if a in positions and b in positions:
                    dx = positions[a][0] - positions[b][0]
                    dy = positions[a][1] - positions[b][1]
                    bond_lengths.append(math.hypot(dx, dy))
            if bond_lengths:
                mean_bl = sum(bond_lengths) / len(bond_lengths)
                if mean_bl > 1e-6:
                    scale = _AROMATIC_CC_BOND / mean_bl
                    positions = {
                        k: (v[0] * scale, v[1] * scale) for k, v in positions.items()
                    }
            return positions
        except Exception:
            return None

    def generate(
        self,
        target_num_carbons: int,
        target_aromaticity: float = 100.0,
        prefer_larger_pahs: bool = True,
        defect_fraction: float = 0.0,
    ) -> CarbonSkeleton:
        """
        Generate a carbon skeleton of approximately *target_num_carbons*.

        Args:
            target_num_carbons: Desired carbon count.
            target_aromaticity: Unused (kept for backward compatibility).
                Aromaticity is an output determined by the ring topology.
            defect_fraction: Probability [0, 1) that any ring addition
                during graph growth is a 5-membered (pentagon) ring.
                0.0 = pure hexagonal PAH (default).
                0.1 = roughly 10% pentagons.
        """
        mol = None

        if defect_fraction <= 0.0:
            # --- Route 1: exact library match (only for pure hexagon mode) ---
            for name, info in self.pahs.items():
                if info["num_carbons"] == target_num_carbons:
                    mol = Chem.Mol(info["mol"])
                    # Attach hex-lattice positions from Compute2DCoords so
                    # downstream geometry generation preserves planarity and
                    # uniform bond lengths.
                    seed_positions = self._compute_seed_positions(mol)
                    if seed_positions:
                        for atom in mol.GetAtoms():
                            idx = atom.GetIdx()
                            if idx in seed_positions:
                                x, y = seed_positions[idx]
                                atom.SetDoubleProp("init_x", x)
                                atom.SetDoubleProp("init_y", y)
                    break

        # --- Route 2: library seed + graph growth (with optional pentagons) ---
        if mol is None:
            mol = self._build_from_seed(target_num_carbons, defect_fraction)

        # --- Fallback: pyrene (no pentagons) ---
        if mol is None:
            mol = Chem.Mol(self.pahs["pyrene"]["mol"])
            seed_positions = self._compute_seed_positions(mol)
            if seed_positions:
                for atom in mol.GetAtoms():
                    idx = atom.GetIdx()
                    if idx in seed_positions:
                        x, y = seed_positions[idx]
                        atom.SetDoubleProp("init_x", x)
                        atom.SetDoubleProp("init_y", y)

        Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
        return self._make_skeleton(mol)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_from_seed(
        self, target: int, defect_fraction: float = 0.0
    ) -> Optional[Chem.Mol]:
        """
        Build a carbon skeleton of *target* carbons using a library seed
        plus ring-growth.

        Pure hexagon mode (defect_fraction == 0):
            Each ring addition adds exactly 4 nodes, so seed selection uses
            the parity constraint (target - seed_size) % 4 == 0.

        Defect mode (defect_fraction > 0):
            Pentagon additions add 3 nodes, so the parity constraint is
            relaxed — any even-count seed works because the growth loop
            handles parity itself.  Up to *_MAX_DEFECT_RETRIES* attempts
            are made with different sub-seeds if kekulization fails.
        """
        _MAX_DEFECT_RETRIES = 5

        if defect_fraction <= 0.0:
            # --- Pure hexagon: strict parity constraint ---
            seed_mol = None
            seed_carbons = 0
            for nC, name in reversed(self._SIZE_INDEX):
                if nC <= target and (target - nC) % 4 == 0 and name in self.pahs:
                    seed_mol = self.pahs[name]["mol"]
                    seed_carbons = nC
                    break

            if seed_mol is None:
                rem = (target - 6) % 4
                if rem != 0:
                    target = target + (4 - rem)
                seed_mol = Chem.MolFromSmiles(PAH_LIBRARY["benzene"]["smiles"])
                seed_carbons = 6

            seed_positions = self._compute_seed_positions(seed_mol)

            if seed_carbons >= target:
                out = Chem.Mol(seed_mol)
                if seed_positions:
                    for atom in out.GetAtoms():
                        idx = atom.GetIdx()
                        if idx in seed_positions:
                            x, y = seed_positions[idx]
                            atom.SetDoubleProp("init_x", x)
                            atom.SetDoubleProp("init_y", y)
                return out

            G = _mol_to_graph(seed_mol)
            # Copy positions so _grow_graph / _fuse_* can mutate in place
            grown_positions = dict(seed_positions) if seed_positions else None
            G = _grow_graph(
                G, target, seed=self.seed, defect_fraction=0.0,
                positions=grown_positions,
            )
            return _graph_to_mol(G, positions=grown_positions)

        else:
            # --- Defect mode: relax parity, retry on kekulization failure ---
            # Use the largest even-count library seed that fits
            seed_mol = None
            seed_carbons = 0
            for nC, name in reversed(self._SIZE_INDEX):
                if nC <= target and name in self.pahs:
                    seed_mol = self.pahs[name]["mol"]
                    seed_carbons = nC
                    break

            if seed_mol is None:
                seed_mol = Chem.MolFromSmiles(PAH_LIBRARY["benzene"]["smiles"])
                seed_carbons = 6

            seed_positions = self._compute_seed_positions(seed_mol)

            if seed_carbons >= target:
                out = Chem.Mol(seed_mol)
                if seed_positions:
                    for atom in out.GetAtoms():
                        idx = atom.GetIdx()
                        if idx in seed_positions:
                            x, y = seed_positions[idx]
                            atom.SetDoubleProp("init_x", x)
                            atom.SetDoubleProp("init_y", y)
                return out

            base_seed = self.seed if self.seed is not None else 0
            for attempt in range(_MAX_DEFECT_RETRIES):
                G = _mol_to_graph(seed_mol)
                grown_positions = dict(seed_positions) if seed_positions else None
                G = _grow_graph(
                    G, target,
                    seed=base_seed + attempt,
                    defect_fraction=defect_fraction,
                    positions=grown_positions,
                )
                mol = _graph_to_mol(G, positions=grown_positions)
                if mol is not None:
                    return mol
                # kekulization failed — try again with a different sub-seed

            return None  # all retries exhausted

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
