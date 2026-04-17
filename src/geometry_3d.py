"""
3D Geometry Generation and Validation

Generate 3D coordinates for biochar molecules and validate geometry.
"""

from typing import Tuple, List, Optional
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.distance import pdist, squareform

from .constants import VDW_RADII, COVALENT_RADII


def _get_excluded_pairs(mol: Chem.Mol) -> set:
    """
    Return the set of (min_idx, max_idx) atom-index pairs that should be
    **excluded** from non-bonded clash detection.

    Standard molecular-mechanics exclusion rules:
    * 1-2 pairs (directly bonded)
    * 1-3 pairs (angle neighbours, sharing one atom)

    In a flat graphene-like sheet the 1-3 distance is
    √3 × L ≈ 2.46 Å (L = 1.42 Å aromatic bond), which is below the
    0.75 × vdW-sum threshold for C-C (2.55 Å).  Without this exclusion
    every second-neighbour carbon pair is flagged as a clash, producing
    hundreds of false positives that the resolver can never fix.
    """
    excluded: set = set()
    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        for nbr in atom.GetNeighbors():
            j = nbr.GetIdx()
            # 1-2 pair
            excluded.add((min(i, j), max(i, j)))
            # 1-3 pairs: all neighbours of j except i
            for nbr2 in nbr.GetNeighbors():
                k = nbr2.GetIdx()
                if k != i:
                    excluded.add((min(i, k), max(i, k)))
    return excluded


def _read_skeleton_positions(mol: Chem.Mol) -> dict:
    """
    Collect hex-lattice positions tagged on the molecule by
    :mod:`carbon_skeleton` (atom properties ``init_x`` / ``init_y``).

    Returns ``{atom_idx: np.array([x, y, 0])}`` for heavy atoms that have
    both props.  Empty dict if no atom is tagged.
    """
    positions: dict = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        if atom.HasProp("init_x") and atom.HasProp("init_y"):
            try:
                x = atom.GetDoubleProp("init_x")
                y = atom.GetDoubleProp("init_y")
            except Exception:
                continue
            positions[atom.GetIdx()] = np.array([x, y, 0.0])
    return positions


def _place_untagged_heavy_atom(
    mol: Chem.Mol, idx: int, placed: dict, fallback_pos: np.ndarray
) -> np.ndarray:
    """
    Compute a 2D position for a heavy atom (typically an oxygen) that
    lacks the hex-lattice tag, by radiating outward from its tagged
    heavy-atom neighbour.
    """
    atom = mol.GetAtomWithIdx(idx)
    heavy_nbrs = [
        n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1
    ]
    parents = [n for n in heavy_nbrs if n.GetIdx() in placed]
    if not parents:
        return fallback_pos.copy()

    parent = parents[0]
    par_pos = placed[parent.GetIdx()]

    # Direction: away from parent's other heavy neighbours
    other = [
        n for n in parent.GetNeighbors()
        if n.GetAtomicNum() != 1 and n.GetIdx() != idx and n.GetIdx() in placed
    ]
    if other:
        nbr_avg = np.mean([placed[n.GetIdx()] for n in other], axis=0)
        direction = par_pos - nbr_avg
    else:
        direction = np.array([1.0, 0.0, 0.0])
    norm = float(np.linalg.norm(direction))
    if norm < 1e-8:
        direction = np.array([1.0, 0.0, 0.0])
    else:
        direction = direction / norm

    # Bond length by element pair
    sym_parent = parent.GetSymbol()
    sym_this = atom.GetSymbol()
    if sym_parent == "C" and sym_this == "O":
        bl = 1.36
    elif sym_parent == "O" and sym_this == "C":
        bl = 1.36
    elif sym_this == "O" or sym_parent == "O":
        bl = 1.36
    else:
        bl = 1.42
    return par_pos + direction * bl


def _heavy_atom_flat_layout(mol: Chem.Mol, seed: Optional[int] = None) -> np.ndarray:
    """
    Produce a flat (z=0) 2D layout for the heavy-atom scaffold of a PAH.

    Uses networkx.kamada_kawai_layout on the heavy-atom bond graph, then
    scales so the mean bond length is ~1.4 Å.  This is dramatically more
    robust on large fused polycyclic biochar graphs than either
    AllChem.Compute2DCoords (which collapses ~25 pairs to distance 0) or
    ETKDGv3 3D embedding + flatten (which folds then loses separation when
    projected).  Biochar graphs are planar by construction, so a planar
    spring-like layout converges cleanly.

    Args:
        mol: molecule with heavy atoms (no explicit hydrogens).
        seed: optional RNG seed for reproducibility of the layout.

    Returns:
        (N, 3) ndarray of 3D coords in Å, with z = 0 for every atom.
    """
    import networkx as nx

    G = nx.Graph()
    for a in mol.GetAtoms():
        G.add_node(a.GetIdx())
    for b in mol.GetBonds():
        G.add_edge(b.GetBeginAtomIdx(), b.GetEndAtomIdx())

    # Kamada-Kawai needs a connected graph; handle components independently.
    n = mol.GetNumAtoms()
    coords = np.zeros((n, 3))

    components = list(nx.connected_components(G))
    x_offset = 0.0
    for comp in components:
        subG = G.subgraph(comp).copy()
        if subG.number_of_nodes() == 1:
            only_node = next(iter(subG.nodes()))
            coords[only_node] = [x_offset, 0.0, 0.0]
            x_offset += 2.0
            continue
        pos = nx.kamada_kawai_layout(subG)
        sub_coords = np.array([pos[i] for i in sorted(pos.keys())])

        # Scale so the mean bond length is ~1.4 Å (aromatic C-C baseline).
        bls = [
            np.linalg.norm(np.array(pos[u]) - np.array(pos[v]))
            for u, v in subG.edges()
        ]
        if bls:
            scale = 1.4 / np.mean(bls)
            sub_coords = sub_coords * scale
        # Shift so component's min-x sits past previous component's max-x.
        min_x = sub_coords[:, 0].min()
        sub_coords[:, 0] += x_offset - min_x
        x_offset = sub_coords[:, 0].max() + 2.0

        for sub_idx, orig_idx in enumerate(sorted(pos.keys())):
            coords[orig_idx] = [sub_coords[sub_idx, 0], sub_coords[sub_idx, 1], 0.0]

    return coords


def _kekulize_or_dearomatize(mol: Chem.Mol) -> Chem.Mol:
    """
    Ensure the returned mol is safe for FF/embedding by either kekulizing it
    or, if that fails (non-kekulizable graph — common for large biochar PAHs
    with odd-parity subgraphs), converting every aromatic bond to SINGLE and
    clearing all aromatic flags.

    After this call, `AllChem.Compute2DCoords`, `AllChem.EmbedMolecule`,
    `AllChem.UFFGetMoleculeForceField`, and `AllChem.MMFFGetMoleculeProperties`
    all work deterministically.
    """
    m = Chem.RWMol(mol)
    try:
        Chem.Kekulize(m, clearAromaticFlags=False)
        return m
    except Exception:
        pass

    # Fallback: strip aromaticity entirely. Downstream FF treats the ring as
    # sp3/sp2 single-bonded; MMFF may still refuse without conjugation info,
    # but UFF works and produces reasonable bond lengths.
    for bond in m.GetBonds():
        if bond.GetBondType() == Chem.BondType.AROMATIC:
            bond.SetBondType(Chem.BondType.SINGLE)
            bond.SetIsAromatic(False)
    for atom in m.GetAtoms():
        atom.SetIsAromatic(False)
    try:
        Chem.SanitizeMol(
            m,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
            ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
        )
    except Exception:
        pass
    return m


def _perpendicular_unit(v: np.ndarray) -> np.ndarray:
    """Return an arbitrary unit vector perpendicular to *v* (v must be non-zero)."""
    v = np.asarray(v, dtype=float)
    # Pick the component with the smallest absolute value to avoid a near-parallel
    # cross product, which would give a numerically unstable result.
    if abs(v[0]) <= abs(v[1]) and abs(v[0]) <= abs(v[2]):
        perp = np.array([0.0, -v[2], v[1]])
    elif abs(v[1]) <= abs(v[2]):
        perp = np.array([-v[2], 0.0, v[0]])
    else:
        perp = np.array([-v[1], v[0], 0.0])
    return perp / np.linalg.norm(perp)


def _optimize_h_positions(
    mol: Chem.Mol,
    all_coords: np.ndarray,
    threshold: float = 1.8,
    num_angles: int = 12,
    allow_out_of_plane: bool = True,
    num_passes: int = 2,
) -> np.ndarray:
    """
    Reduce steric clashes for hydrogen atoms by rotating them around the bond
    from their grandparent to their parent atom.

    For each H atom whose nearest non-bonded neighbour is closer than
    *threshold* Å, the H is swept through *num_angles* equally-spaced dihedral
    angles around the grandparent→parent bond axis (e.g. the C–O axis for a
    phenolic OH, or the ring-C bond for a peripheral C-H).  The position that
    maximises the minimum non-bonded distance is kept.

    When *allow_out_of_plane* is True, off-plane positions (z ≠ 0, valid
    because hex-lattice molecules lie flat in the xy plane) are sampled at
    ±0.3, ±0.6, and ±0.9 Å above/below the sheet.  This is the primary fix
    for OH hydrogens that are pushed into the ring plane by geometric placement.

    *num_passes* sweeps are performed (default 2).  A second pass is necessary
    because optimising H_i changes the non-bonded environment seen by H_j
    (which may be H_i's nearest neighbour), so H_j may need to be re-evaluated
    after H_i has moved.

    Clashing H atoms are processed sequentially within each pass; each update
    is immediately visible to subsequent H optimisations in the same pass.

    Args:
        mol:               RDKit molecule (used for connectivity only).
        all_coords:        (N, 3) coordinate array in Å — **modified in place**.
        threshold:         Minimum non-bonded distance (Å) that triggers
                           optimisation.  Default 1.8 Å.
        num_angles:        Number of in-plane dihedral samples (default 12 →
                           every 30°).
        allow_out_of_plane: If True, also sample positions above/below the
                           molecular plane.
        num_passes:        Number of full sweeps over all H atoms (default 2).

    Returns:
        The (possibly modified) coordinate array (same object as *all_coords*).
    """
    excluded = _get_excluded_pairs(mol)

    def _min_nb_dist(h_idx: int, h_pos: np.ndarray) -> float:
        """Minimum distance from *h_pos* to any non-excluded, non-self atom."""
        min_d = float("inf")
        for j in range(mol.GetNumAtoms()):
            if j == h_idx:
                continue
            pair = (min(h_idx, j), max(h_idx, j))
            if pair in excluded:
                continue
            d = float(np.linalg.norm(h_pos - all_coords[j]))
            if d < min_d:
                min_d = d
        return min_d

    # Out-of-plane z-offsets (Å): 0 = flat (in-plane); ±Δz = above/below sheet.
    z_offsets: List[float] = [0.0]
    if allow_out_of_plane:
        z_offsets += [0.3, -0.3, 0.6, -0.6, 0.9, -0.9]

    angles = np.linspace(0.0, 2.0 * np.pi, num_angles, endpoint=False)

    for _pass in range(num_passes):
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 1:
                continue
            h_idx = atom.GetIdx()
            h_pos = all_coords[h_idx].copy()

            # Only optimise H atoms that currently have a short contact
            if _min_nb_dist(h_idx, h_pos) >= threshold:
                continue

            # Parent (O or C that this H is bonded to)
            nbrs = atom.GetNeighbors()
            if not nbrs:
                continue
            parent = nbrs[0]
            par_idx = parent.GetIdx()
            par_pos = all_coords[par_idx].copy()

            # Grandparent: first heavy-atom neighbour of parent (excludes self)
            grand_nbrs = [
                n for n in parent.GetNeighbors()
                if n.GetAtomicNum() != 1 and n.GetIdx() != h_idx
            ]
            if not grand_nbrs:
                continue
            grand_pos = all_coords[grand_nbrs[0].GetIdx()].copy()

            # Rotation axis = grandparent → parent (normalised)
            axis = par_pos - grand_pos
            a_norm = float(np.linalg.norm(axis))
            if a_norm < 1e-8:
                continue
            axis = axis / a_norm

            # Decompose h_vec into components along / perpendicular to the axis
            h_vec = h_pos - par_pos
            bl = float(np.linalg.norm(h_vec))
            if bl < 1e-8:
                continue

            h_par = np.dot(h_vec, axis) * axis    # component along rotation axis
            h_perp = h_vec - h_par                 # component in the sweep plane
            perp_norm = float(np.linalg.norm(h_perp))

            # Orthonormal basis in the sweep plane
            u_hat = (h_perp / perp_norm if perp_norm > 1e-8
                     else _perpendicular_unit(axis))
            v_hat = np.cross(axis, u_hat)

            best_pos = h_pos.copy()
            best_d = _min_nb_dist(h_idx, h_pos)

            for z_off in z_offsets:
                z_vec = np.array([0.0, 0.0, z_off])
                for theta in angles:
                    new_perp = perp_norm * (
                        np.cos(theta) * u_hat + np.sin(theta) * v_hat
                    )
                    candidate = par_pos + h_par + new_perp + z_vec
                    # Rescale to the original O-H / C-H bond length
                    cv = candidate - par_pos
                    cv_norm = float(np.linalg.norm(cv))
                    if cv_norm > 1e-8:
                        candidate = par_pos + cv * (bl / cv_norm)
                    d = _min_nb_dist(h_idx, candidate)
                    if d > best_d:
                        best_d = d
                        best_pos = candidate.copy()

            all_coords[h_idx] = best_pos

    return all_coords


class CoordinateGenerator:
    """Generate valid 3D coordinates for molecules."""

    def __init__(self, seed: int = None):
        self.seed = seed
        # Set to True after generate_3d_coordinates() if hex-lattice positions
        # (from carbon_skeleton init_x / init_y props) were used.  When True,
        # callers should skip any subsequent force-field minimization because
        # running MMFF/UFF on the dearomatized mol (all-SINGLE bonds) would
        # incorrectly treat sp2 carbons as sp3 and buckle the flat sheet.
        self.used_hex_lattice: bool = False
        if seed is not None:
            np.random.seed(seed)

    def generate_3d_coordinates(
        self,
        mol: Chem.Mol,
        force_aromatic_planarity: bool = True,
        num_conformers: int = 1,
        use_distance_geometry: bool = True,
        max_embedding_attempts: int = 3,
        max_ff_iterations: int = None,
    ) -> Tuple[Chem.Mol, np.ndarray]:
        """
        Generate 3D coordinates for molecule with enhanced relaxation.

        Args:
            mol: RDKit molecule
            force_aromatic_planarity: If True, enforce planarity on aromatic rings
            num_conformers: Number of conformers to generate
            use_distance_geometry: Use distance geometry for coordinate generation
            max_embedding_attempts: Number of independent embedding attempts (default 3)
            max_ff_iterations: Max force field iterations (default: 500 for large molecules, 300 for small)

        Returns:
            (mol_with_coords, coords_array)
        """
        # Ensure aromaticity is properly perceived before embedding.
        # Do NOT call Kekulize here: it fails for large graphene-like systems and
        # corrupts the aromatic bond types.  ETKDGv3 can embed aromatic bond types
        # natively; UFF/MMFF operate on the sanitized mol after embedding.
        try:
            Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
        except Exception:
            pass

        # Add hydrogens if not already present
        if mol.GetNumHeavyAtoms() == mol.GetNumAtoms():
            mol = Chem.AddHs(mol)

        # Adaptive force field iterations based on molecule size
        if max_ff_iterations is None:
            num_heavy_atoms = mol.GetNumHeavyAtoms()
            if num_heavy_atoms > 300:
                max_ff_iterations = 500
            elif num_heavy_atoms > 100:
                max_ff_iterations = 400
            else:
                max_ff_iterations = 300

        # For large flat aromatic systems (> 80 heavy atoms), ETKDGv3 frequently
        # folds or collapses the structure.  Use 2D-first embedding instead:
        #   1. Compute 2D layout (flat aromatic ring system)
        #   2. Promote to 3D with z = 0
        #   3. Perturb non-ring atoms slightly in z for FF convergence
        #   4. Run force field minimization
        use_2d_first = mol.GetNumHeavyAtoms() > 80

        best_mol = None
        best_energy = float('inf')
        best_converged = False
        self.used_hex_lattice = False  # reset each call

        if use_2d_first:
            # Produce a kekulized (or de-aromatized) working copy so that
            # Compute2DCoords, FF setup, etc. don't raise KekulizeException
            # on non-kekulizable biochar graphs (all 100+ atoms flagged).
            mol_copy = Chem.RWMol(_kekulize_or_dearomatize(mol))
            try:
                # 2D-first embedding for large flat aromatics:
                #  1. Compute 2D layout on the heavy-atom-only scaffold
                #  2. Place heavy atoms at z=0 (flat PAH)
                #  3. Compute H positions geometrically (not from 2D conventions)
                #  4. FF minimization to relax H and functional groups
                heavy_atom_indices = [
                    a.GetIdx() for a in mol_copy.GetAtoms() if a.GetAtomicNum() != 1
                ]
                mol_no_h = Chem.RemoveAllHs(mol_copy.GetMol(), sanitize=False)

                # --- Layout strategy (priority order) ---
                # 1. Hex-lattice positions threaded from carbon_skeleton.py
                #    via atom props init_x / init_y.  Gives perfect 1.42 Å
                #    aromatic bonds and planar hex geometry for all sizes.
                # 2. Kamada-Kawai spring layout (fallback when props absent).
                #    Works well up to ~100 atoms; degrades at 200-500 atoms.
                skel_pos = _read_skeleton_positions(mol_no_h)
                n_carbons_nh = sum(
                    1 for a in mol_no_h.GetAtoms() if a.GetAtomicNum() == 6
                )
                use_skel = bool(skel_pos) and len(skel_pos) >= n_carbons_nh
                if use_skel:
                    self.used_hex_lattice = True

                flat_heavy = np.zeros((mol_no_h.GetNumAtoms(), 3))
                if use_skel:
                    # Carbon positions come directly from hex-lattice geometry.
                    # Oxygens / other heteroatoms are placed by radiating from
                    # their tagged parent carbon.
                    placed: dict = {}
                    # First pass: place carbons (all tagged)
                    for atom in mol_no_h.GetAtoms():
                        idx = atom.GetIdx()
                        if idx in skel_pos:
                            flat_heavy[idx] = skel_pos[idx]
                            placed[idx] = skel_pos[idx]
                    # Second pass: place untagged heteroatoms (oxygens)
                    for atom in mol_no_h.GetAtoms():
                        idx = atom.GetIdx()
                        if idx not in placed:
                            pos = _place_untagged_heavy_atom(
                                mol_no_h, idx, placed,
                                fallback_pos=flat_heavy[max(placed, key=lambda k: k)] + np.array([1.36, 0.0, 0.0])
                                if placed else np.zeros(3),
                            )
                            flat_heavy[idx] = pos
                            placed[idx] = pos
                else:
                    # Fallback: kamada-kawai spring layout
                    flat_heavy = _heavy_atom_flat_layout(mol_no_h, seed=self.seed)

                # Map mol_no_h atom indices (0..Nheavy-1) back to mol_copy
                # heavy-atom indices (which can interleave with H indices).
                heavy_pos = {}  # mol_copy_idx -> np.array([x, y, z])
                for nh_idx, orig_idx in enumerate(heavy_atom_indices):
                    heavy_pos[orig_idx] = flat_heavy[nh_idx].copy()

                # Build full coordinate array; start with heavy atoms
                all_coords = np.zeros((mol_copy.GetNumAtoms(), 3))
                for idx, coord in heavy_pos.items():
                    all_coords[idx] = coord

                # Place H atoms: radiate from parent at correct bond length
                for atom in mol_copy.GetAtoms():
                    if atom.GetAtomicNum() != 1:
                        continue
                    idx = atom.GetIdx()
                    parent = atom.GetNeighbors()[0]
                    par_idx = parent.GetIdx()
                    par_pos = heavy_pos.get(par_idx, all_coords[par_idx])

                    # Direction: opposite to average of other neighbor directions
                    other_nbrs = [
                        n.GetIdx() for n in parent.GetNeighbors() if n.GetIdx() != idx
                    ]
                    if other_nbrs:
                        nbr_avg = np.mean(
                            [heavy_pos.get(n, all_coords[n]) for n in other_nbrs], axis=0
                        )
                        h_dir = par_pos - nbr_avg
                    else:
                        h_dir = np.array([1.0, 0.0, 0.0])
                    norm = np.linalg.norm(h_dir)
                    if norm < 1e-8:
                        h_dir = np.array([1.0, 0.0, 0.0])
                    else:
                        h_dir = h_dir / norm

                    # Choose bond length by parent-H element
                    par_elem = parent.GetAtomicNum()
                    bl = 1.09 if par_elem == 6 else 0.96  # C-H or O-H
                    all_coords[idx] = par_pos + h_dir * bl

                # Rotate clashing H atoms to minimise non-bonded contacts.
                # Called for both the hex-lattice (use_skel=True) and the
                # Kamada-Kawai fallback paths: both produce a flat z=0 layout,
                # so out-of-plane sampling (z ≠ 0) is safe and effective.
                _optimize_h_positions(mol_copy, all_coords)

                # Create 3D conformer
                conf3d = Chem.Conformer(mol_copy.GetNumAtoms())
                for i, xyz in enumerate(all_coords):
                    conf3d.SetAtomPosition(i, xyz.tolist())
                mol_copy.RemoveAllConformers()
                mol_copy.AddConformer(conf3d, assignId=True)
                best_mol = mol_copy

                # FF minimization:
                # When hex-lattice positions are used (use_skel=True) we SKIP the
                # force field entirely.  The hex lattice gives perfect 1.42 Å CC
                # bonds; running MMFF/UFF on the dearomatized mol (all-SINGLE bonds)
                # would incorrectly treat carbons as sp3 and buckle the sheet into a
                # 3D blob.  H and O positions from geometric placement are good enough
                # for GROMACS input — the production MD run will relax them further.
                #
                # For the KK fallback (use_skel=False) the FF is still useful to
                # refine bond lengths from the spring layout.
                if not use_skel:
                    try:
                        mol_for_ff = mol_copy.GetMol()
                        mmff_props = AllChem.MMFFGetMoleculeProperties(mol_for_ff)
                        ff = (
                            AllChem.MMFFGetMoleculeForceField(mol_for_ff, mmff_props)
                            if mmff_props is not None
                            else None
                        )
                        if ff is None:
                            ff = AllChem.UFFGetMoleculeForceField(mol_for_ff)
                        if ff is not None:
                            minimize_result = ff.Minimize(maxIts=max_ff_iterations, forceTol=1e-4)
                            best_converged = minimize_result == 0
                            best_energy = ff.CalcEnergy()
                            # Pull minimized coordinates back into mol_copy
                            conf_min = mol_for_ff.GetConformer(0)
                            conf_target = mol_copy.GetConformer(0)
                            for i in range(mol_for_ff.GetNumAtoms()):
                                conf_target.SetAtomPosition(i, conf_min.GetAtomPosition(i))
                    except Exception:
                        pass
            except Exception:
                best_mol = None
        else:
            # Kekulize (or de-aromatize on failure) once outside the retry loop
            # so EmbedMolecule / MMFF don't throw on the same issue repeatedly.
            mol_safe = _kekulize_or_dearomatize(mol)
            for attempt in range(max_embedding_attempts):
                seed = (self.seed + attempt) if self.seed is not None else attempt
                mol_copy = Chem.RWMol(mol_safe)

                result = -1
                for params_fn in [
                    lambda: AllChem.ETKDGv3(),
                    lambda: AllChem.ETKDGv2(),
                    lambda: AllChem.EmbedParameters(),
                ]:
                    params = params_fn()
                    params.randomSeed = seed
                    params.maxIterations = 200
                    result = AllChem.EmbedMolecule(mol_copy, params)
                    if result == 0:
                        break

                if result != 0:
                    continue

                ff_converged = False
                ff_energy = float('inf')
                # MMFF first — requires MMFFMolProperties, else call throws.
                try:
                    mmff_props = AllChem.MMFFGetMoleculeProperties(mol_copy)
                    if mmff_props is not None:
                        ff = AllChem.MMFFGetMoleculeForceField(mol_copy, mmff_props)
                        if ff is not None:
                            minimize_result = ff.Minimize(maxIts=max_ff_iterations, forceTol=1e-6)
                            ff_converged = minimize_result == 0
                            ff_energy = ff.CalcEnergy()
                except Exception:
                    pass

                if ff_energy == float('inf'):
                    try:
                        ff = AllChem.UFFGetMoleculeForceField(mol_copy)
                        if ff is not None:
                            minimize_result = ff.Minimize(maxIts=max_ff_iterations, forceTol=1e-6)
                            ff_converged = minimize_result == 0
                            ff_energy = ff.CalcEnergy()
                    except Exception:
                        pass

                if ff_energy < best_energy:
                    best_energy = ff_energy
                    best_mol = mol_copy
                    best_converged = ff_converged

        # Use best embedding, or fall back to 2D if all attempts failed
        if best_mol is not None:
            mol = best_mol
        else:
            # Final fallback: use 2D coordinates with z=0
            AllChem.Compute2DCoords(mol)

        # Extract coordinates
        try:
            conf = mol.GetConformer(0)
            coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
        except Exception:
            coords = np.zeros((mol.GetNumAtoms(), 3))

        # Apply planarity constraints to aromatic rings if requested
        if force_aromatic_planarity:
            coords = self._enforce_aromatic_planarity(mol, coords)

        return mol, coords

    def _enforce_aromatic_planarity(
        self, mol: Chem.Mol, coords: np.ndarray, tolerance: float = 0.5
    ) -> np.ndarray:
        """
        Enforce planarity on aromatic ring systems.

        Args:
            mol: RDKit molecule
            coords: Current 3D coordinates
            tolerance: Maximum deviation from plane (Angstrom)

        Returns:
            Adjusted coordinates
        """
        # Find aromatic rings
        ring_info = mol.GetRingInfo()
        aromatic_rings = [
            ring for ring in ring_info.AtomRings()
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        ]

        adjusted_coords = coords.copy()

        for ring_atoms in aromatic_rings:
            if len(ring_atoms) < 3:
                continue

            # Get ring atom coordinates
            ring_coords = coords[list(ring_atoms)]

            # Fit plane to ring atoms
            plane_normal, plane_point = self._fit_plane(ring_coords)

            # Project ring atoms onto plane
            for i, atom_idx in enumerate(ring_atoms):
                point = adjusted_coords[atom_idx]
                projected = point - np.dot(point - plane_point, plane_normal) * plane_normal
                adjusted_coords[atom_idx] = projected

        return adjusted_coords

    @staticmethod
    def _fit_plane(points: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Fit a plane to 3D points using PCA.

        Returns:
            (plane_normal, plane_point)
        """
        center = points.mean(axis=0)
        centered_points = points - center

        # SVD
        u, s, vt = np.linalg.svd(centered_points)
        normal = vt[2]  # Normal is the smallest singular vector

        return normal, center

    def validate_and_relax(
        self, mol: Chem.Mol, coords: np.ndarray, max_iterations: int = 100
    ) -> Tuple[np.ndarray, bool]:
        """
        Validate and relax coordinates using force field.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates
            max_iterations: Maximum relaxation iterations

        Returns:
            (relaxed_coords, converged)
        """
        # Set coordinates in molecule
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, coord in enumerate(coords):
            conf.SetAtomPosition(i, tuple(coord))
        mol.RemoveAllConformers()
        mol.AddConformer(conf)

        # Try MMFF94 force field.
        # NOTE: MMFFGetMoleculeForceField requires MMFFMolProperties.
        try:
            mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
            if mmff_props is not None:
                ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props)
                if ff is not None:
                    result = ff.Minimize(maxIts=max_iterations)
                    converged = result == 0
                    conf = mol.GetConformer(0)
                    relaxed_coords = np.array([
                        list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())
                    ])
                    return relaxed_coords, converged
        except Exception:
            pass

        # Fallback to UFF
        try:
            ff = AllChem.UFFGetMoleculeForceField(mol)
            if ff is not None:
                result = ff.Minimize(maxIts=max_iterations)
                converged = result == 0
                conf = mol.GetConformer(0)
                relaxed_coords = np.array([
                    list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())
                ])
                return relaxed_coords, converged
        except Exception:
            pass

        return coords, False


class ClashResolver:
    """Resolve steric clashes in molecular geometries."""

    @staticmethod
    def resolve_clashes(
        mol: Chem.Mol,
        coords: np.ndarray,
        max_iterations: int = 10,
        displacement_step: float = 0.1,
        use_vdw_radii: bool = True,
    ) -> np.ndarray:
        """
        Iteratively resolve steric clashes by displacing atoms away from each other.

        Algorithm:
        - Detect pairs of atoms that are too close (violating Van der Waals radii)
        - Displace the lighter atom (usually H) away from its neighbor
        - Repeat until no clashes remain or max iterations reached

        Args:
            mol: RDKit molecule
            coords: 3D coordinates
            max_iterations: Maximum number of clash resolution iterations
            displacement_step: Step size for atom displacement (Angstrom)
            use_vdw_radii: If True, use Van der Waals radii; else use fixed threshold

        Returns:
            Adjusted coordinates with reduced clashes
        """
        adjusted_coords = coords.copy()

        for iteration in range(max_iterations):
            clashes = ClashResolver._detect_clashes(mol, adjusted_coords, use_vdw_radii)

            if not clashes:
                # No more clashes
                break

            # Displace clashing atoms
            for atom_i, atom_j in clashes:
                # Prefer to move the lighter atom (usually H)
                mass_i = mol.GetAtomWithIdx(atom_i).GetMass()
                mass_j = mol.GetAtomWithIdx(atom_j).GetMass()

                if mass_i < mass_j:
                    moving_atom = atom_i
                else:
                    moving_atom = atom_j

                # Calculate displacement direction (away from the other atom)
                if moving_atom == atom_i:
                    other_atom = atom_j
                else:
                    other_atom = atom_i

                direction = adjusted_coords[moving_atom] - adjusted_coords[other_atom]
                distance = np.linalg.norm(direction)

                if distance > 1e-6:  # Avoid division by zero
                    direction_normalized = direction / distance
                    adjusted_coords[moving_atom] += displacement_step * direction_normalized

        return adjusted_coords

    @staticmethod
    def _detect_clashes(
        mol: Chem.Mol,
        coords: np.ndarray,
        use_vdw_radii: bool = True,
    ) -> List[Tuple[int, int]]:
        """
        Detect steric clashes between atoms.

        1-2 (bonded) and 1-3 (angle) atom pairs are excluded because in flat
        graphene-like sheets the 1-3 C-C distance is ≈2.46 Å — below the
        0.75 × vdW threshold (2.55 Å) — so they would produce hundreds of
        false positives.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates
            use_vdw_radii: If True, use Van der Waals radii; else use fixed 2.0 Å threshold

        Returns:
            List of (atom_i, atom_j) pairs with clashes
        """
        # Pre-compute 1-2 and 1-3 exclusions (see _get_excluded_pairs).
        excluded = _get_excluded_pairs(mol)
        clashes = []
        distances = squareform(pdist(coords))

        for i in range(mol.GetNumAtoms()):
            for j in range(i + 1, mol.GetNumAtoms()):
                if (i, j) in excluded:
                    continue

                if use_vdw_radii:
                    # A "clash" is real overlap — use 0.75 × vdW-sum as the hard
                    # threshold. Adding +0.2 to the vdW sum (as the old code did)
                    # flagged every peri-H in a PAH as a clash, and the resolver
                    # then destroyed bond geometry trying to fix "clashes" that
                    # were just normal aromatic-ring proximity.
                    atom_i = mol.GetAtomWithIdx(i)
                    atom_j = mol.GetAtomWithIdx(j)
                    symbol_i = atom_i.GetSymbol()
                    symbol_j = atom_j.GetSymbol()
                    r_vdw_i = VDW_RADII.get(symbol_i, 1.70)
                    r_vdw_j = VDW_RADII.get(symbol_j, 1.70)
                    min_distance = 0.75 * (r_vdw_i + r_vdw_j)
                else:
                    # Use fixed threshold
                    min_distance = 1.5

                if distances[i, j] < min_distance:
                    clashes.append((i, j))

        return clashes


class GeometryValidator:
    """Validate 3D molecular geometry."""

    @staticmethod
    def validate_geometry(
        mol: Chem.Mol, coords: np.ndarray
    ) -> Tuple[bool, List[str]]:
        """
        Validate 3D geometry for issues.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates

        Returns:
            (is_valid, error_messages)
        """
        errors = []

        # Check for NaN coordinates
        if np.any(np.isnan(coords)):
            errors.append("NaN values found in coordinates")
            return False, errors

        # Check for steric clashes
        clash_errors = GeometryValidator._check_steric_clashes(mol, coords)
        errors.extend(clash_errors)

        # Check bond lengths
        bond_errors = GeometryValidator._validate_bond_lengths(mol, coords)
        errors.extend(bond_errors)

        # Check for unphysical geometries
        if np.max(coords) > 1000 or np.min(coords) < -1000:
            errors.append("Coordinates are unphysically large")

        return len(errors) == 0, errors

    @staticmethod
    def _check_steric_clashes(
        mol: Chem.Mol, coords: np.ndarray, clash_threshold: float = 2.0
    ) -> List[str]:
        """
        Check for steric clashes between atoms with enhanced reporting.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates
            clash_threshold: Minimum allowed distance (Angstrom)

        Returns:
            List of error messages with clash severity
        """
        # Pre-compute 1-2 and 1-3 exclusions: in a flat graphene sheet the
        # 1-3 C-C distance (≈2.46 Å) falls below the 0.75×vdW threshold
        # (2.55 Å), so without this exclusion every angle neighbour pair is
        # reported as a clash.
        excluded = _get_excluded_pairs(mol)

        errors = []
        clashes = []

        # Calculate pairwise distances
        distances = squareform(pdist(coords))

        for i in range(mol.GetNumAtoms()):
            for j in range(i + 1, mol.GetNumAtoms()):
                if (i, j) in excluded:
                    continue

                atom_i = mol.GetAtomWithIdx(i)
                atom_j = mol.GetAtomWithIdx(j)

                # Get Van der Waals radii. A real clash is atom overlap — use
                # 0.75 × vdW-sum so that normal non-bonded proximity (e.g. peri
                # H-H ≈ 1.95 Å in PAHs) is not flagged as a clash.
                symbol_i = atom_i.GetSymbol()
                symbol_j = atom_j.GetSymbol()
                r_vdw_i = VDW_RADII.get(symbol_i, 1.70)
                r_vdw_j = VDW_RADII.get(symbol_j, 1.70)
                min_distance = 0.75 * (r_vdw_i + r_vdw_j)

                distance = distances[i, j]
                if distance < min_distance:
                    severity = min_distance - distance  # How far below minimum
                    clash_type = "H-H" if symbol_i == "H" and symbol_j == "H" else \
                                 "H-C" if symbol_i in ["H", "C"] and symbol_j in ["H", "C"] else \
                                 "Other"

                    clashes.append((i, j, distance, min_distance, severity, clash_type))
                    errors.append(
                        f"Steric clash: atoms {i} and {j} "
                        f"(distance: {distance:.2f}, min: {min_distance:.2f}, "
                        f"severity: {severity:.2f} Å, type: {clash_type})"
                    )

        return errors

    @staticmethod
    def _validate_bond_lengths(
        mol: Chem.Mol, coords: np.ndarray
    ) -> List[str]:
        """
        Validate bond lengths are reasonable.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates

        Returns:
            List of error messages
        """
        errors = []

        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()

            atom_i = mol.GetAtomWithIdx(i)
            atom_j = mol.GetAtomWithIdx(j)

            symbol_i = atom_i.GetSymbol()
            symbol_j = atom_j.GetSymbol()

            # Calculate distance
            distance = np.linalg.norm(coords[i] - coords[j])

            # Expected distance: sum of covalent radii
            r_cov_i = COVALENT_RADII.get(symbol_i, 0.76)
            r_cov_j = COVALENT_RADII.get(symbol_j, 0.76)
            expected_distance = r_cov_i + r_cov_j

            # Allow 20% deviation
            min_dist = expected_distance * 0.8
            max_dist = expected_distance * 1.5

            if distance < min_dist or distance > max_dist:
                errors.append(
                    f"Unusual bond length: atoms {i}-{j} "
                    f"(distance: {distance:.2f}, expected: {expected_distance:.2f})"
                )

        return errors

    @staticmethod
    def check_aromaticity(mol: Chem.Mol) -> Tuple[bool, int]:
        """
        Check if aromaticity is preserved.

        Returns:
            (is_valid, num_aromatic_atoms)
        """
        num_aromatic = sum(
            1 for atom in mol.GetAtoms() if atom.GetIsAromatic()
        )
        return num_aromatic > 0, num_aromatic

    @staticmethod
    def measure_ring_planarity(mol: Chem.Mol, coords: np.ndarray) -> Tuple[float, str]:
        """
        Measure planarity of aromatic rings.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates

        Returns:
            (avg_deviation, assessment)
        """
        ring_info = mol.GetRingInfo()
        aromatic_rings = [
            ring for ring in ring_info.AtomRings()
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        ]

        if not aromatic_rings:
            return 0.0, "No aromatic rings"

        deviations = []

        for ring_atoms in aromatic_rings:
            if len(ring_atoms) < 3:
                continue

            ring_coords = coords[list(ring_atoms)]
            normal, point = CoordinateGenerator._fit_plane(ring_coords)

            # Calculate deviations from plane
            for coord in ring_coords:
                deviation = abs(np.dot(coord - point, normal))
                deviations.append(deviation)

        avg_deviation = np.mean(deviations) if deviations else 0.0
        if avg_deviation < 0.1:
            assessment = "Good planarity"
        elif avg_deviation < 0.3:
            assessment = "Acceptable planarity"
        else:
            assessment = "Poor planarity"

        return avg_deviation, assessment
