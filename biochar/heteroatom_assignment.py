"""
Heteroatom Assignment

Assign hydrogen, oxygen, and nitrogen atoms to carbon skeletons to satisfy
compositional ratios. Respects valence constraints for all atoms.
"""

import logging
import random
from typing import List, Dict, Tuple, Optional, Set
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

from rdkit import Chem

from .constants import FUNCTIONAL_GROUPS
from .valence import ValenceValidator


@dataclass
class CompositionResult:
    """
    Mutable composition record threaded through the generation pipeline.

    Created by :class:`OxygenAssigner` and updated in-place by
    :class:`HydrogenAssigner`, so callers always hold one authoritative object.
    """

    num_carbons: int = 0
    num_hydrogens: int = 0
    num_oxygens: int = 0
    num_nitrogens: int = 0
    num_sulfurs: int = 0
    H_C_ratio: float = 0.0
    O_C_ratio: float = 0.0
    N_C_ratio: float = 0.0
    S_C_ratio: float = 0.0
    # Ring-substituting nitrogen census.  These are a subset of num_nitrogens
    # (which counts every N atom, including pendant amino N).
    num_pyridinic: int = 0
    num_pyrrolic: int = 0
    num_graphitic: int = 0
    functional_groups: Dict[str, int] = field(default_factory=dict)
    placed_counts: Dict[str, int] = field(default_factory=dict)   # groups actually placed
    requested_counts: Dict[str, int] = field(default_factory=dict)  # groups requested
    # pH-dependent protonation census, set by ``ProtonationAssigner``.
    # ``net_charge`` is the molecule's total formal charge -- 0 when no pH was
    # requested, and what ``genion -neutral`` must balance at solvation time.
    # ``titratable_counts`` is every site that *could* ionize;
    # ``ionized_counts`` is the subset that actually did.  Both are reported so
    # the per-site sampling stays visible: near a pKa a single structure is one
    # draw, not the ensemble mean.
    net_charge: int = 0
    ionized_counts: Dict[str, int] = field(default_factory=dict)
    titratable_counts: Dict[str, int] = field(default_factory=dict)
    # Set by HydrogenAssigner when the requested H/C exceeds what the built
    # skeleton can physically carry.  ``h_c_ceiling`` is the maximum H/C
    # achievable for this molecule (all free edge valences saturated with H);
    # ``h_c_target_unreachable`` flags that the target was above that ceiling and
    # the structure was capped rather than reaching the request.
    h_c_ceiling: Optional[float] = None
    h_c_target_unreachable: bool = False

    @property
    def molecular_formula(self) -> str:
        """Hill-order molecular formula string, e.g. 'C24H12O2N1'."""
        parts = [f"C{self.num_carbons}"]
        if self.num_hydrogens:
            parts.append(f"H{self.num_hydrogens}")
        if self.num_oxygens:
            parts.append(f"O{self.num_oxygens}")
        if self.num_nitrogens:
            parts.append(f"N{self.num_nitrogens}")
        if self.num_sulfurs:
            parts.append(f"S{self.num_sulfurs}")
        return "".join(parts)

    @property
    def molecular_weight(self) -> float:
        """Molecular weight in g/mol."""
        return (
            self.num_carbons * 12.011
            + self.num_hydrogens * 1.008
            + self.num_oxygens * 15.999
            + self.num_nitrogens * 14.007
            + self.num_sulfurs * 32.06
        )


# Backward-compatible alias
CompositionInfo = CompositionResult


# Flags for sanitisation that skip SANITIZE_KEKULIZE and SANITIZE_SETAROMATICITY.
#
# SANITIZE_KEKULIZE:      fails on large graphene-like systems and partially
#                         corrupts aromatic bond types before raising.
#
# SANITIZE_SETAROMATICITY: uses RDKit's own aromaticity model (RDKIT, not MDL),
#                          which marks ether C-O bonds as AROMATIC when the bridge
#                          closes a ring that satisfies Hückel's rule (furan-like).
#                          We re-perceive aromaticity separately with
#                          SetAromaticity(mol, AROMATICITY_MDL) which correctly
#                          leaves C-O bonds as SINGLE while still marking ring C-C
#                          bonds as AROMATIC.
_SAFE_FLAGS = (
    Chem.SanitizeFlags.SANITIZE_FINDRADICALS
    | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
    | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
    | Chem.SanitizeFlags.SANITIZE_SYMMRINGS
    | Chem.SanitizeFlags.SANITIZE_ADJUSTHS
)


def _safe_sanitize(mol: Chem.Mol) -> None:
    """Sanitize mol in-place without Kekulization, then re-perceive aromaticity."""
    try:
        Chem.SanitizeMol(mol, _SAFE_FLAGS)
    except Exception:
        pass
    try:
        Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
    except Exception:
        pass


def _fix_heteroatom_bond_types(mol: Chem.Mol) -> Chem.Mol:
    """
    Correct bonds that involve heteroatoms (non-C) and were spuriously marked
    AROMATIC by aromaticity perception.

    This can happen when an ether oxygen bridges two aromatic edge carbons and
    the resulting large ring satisfies Hückel's rule (e.g., furan-like 5-ring).
    Oxygen, hydrogen, and carboxyl carbons outside the aromatic ring system should
    never have AROMATIC bonds.

    Exception — ring-substituted nitrogen (pyridinic / pyrrolic / graphitic) is a
    *genuine* aromatic ring member.  Its in-ring C-N bonds are correctly AROMATIC
    and must NOT be reset to SINGLE.  Such an N is identified by being aromatic
    and a ring member with at least two heavy (non-H) neighbours.

    Returns a new mol with corrected bond types.
    """
    ring_info = mol.GetRingInfo()

    def _is_aromatic_ring_nitrogen(atom: Chem.Atom) -> bool:
        if atom.GetAtomicNum() != 7:
            return False
        if not atom.GetIsAromatic():
            return False
        if ring_info.NumAtomRings(atom.GetIdx()) == 0:
            return False
        heavy = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() != 1)
        return heavy >= 2

    rwmol = Chem.RWMol(mol)
    changed = False
    for bond in rwmol.GetBonds():
        if bond.GetBondTypeAsDouble() == 1.5:  # currently AROMATIC
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Only fix bonds where at least one atom is a heteroatom (not C)
            if a1.GetAtomicNum() != 6 or a2.GetAtomicNum() != 6:
                # Leave aromatic ring nitrogen alone — it belongs in the ring.
                if _is_aromatic_ring_nitrogen(a1) or _is_aromatic_ring_nitrogen(a2):
                    continue
                bond.SetBondType(Chem.BondType.SINGLE)
                if a1.GetAtomicNum() != 6:
                    a1.SetIsAromatic(False)
                if a2.GetAtomicNum() != 6:
                    a2.SetIsAromatic(False)
                changed = True
    return rwmol.GetMol() if changed else mol


def attach_aliphatic_carbons(
    mol: Chem.Mol, n_aliphatic: int, rng: random.Random
) -> Chem.Mol:
    """
    Attach *n_aliphatic* pendant sp3 carbons (methyl groups) to distinct
    aromatic edge carbons and return the new molecule.

    A fully aromatic (sp2) flake caps hydrogen only on its perimeter, so its
    H/C is bounded by the perimeter/area ratio.  Reaching the higher H/C of
    low-temperature biochar requires genuine aliphatic carbon; each pendant
    methyl converts one edge C-H into C-CH3 (net +2 H, +1 C).  The pendant
    carbons are left bare -- their three hydrogens are added later by
    :class:`HydrogenAssigner`, OPLS typing maps them to CT with HC hydrogens,
    and the geometry stage places them like any other off-lattice atom.

    Attaches at most as many methyls as there are free aromatic edge sites.
    """
    if n_aliphatic <= 0:
        return mol
    sites = [
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIsAromatic() and a.GetDegree() < 3
    ]
    if not sites:
        return mol
    rng.shuffle(sites)
    sites = sites[:n_aliphatic]

    emol = Chem.RWMol(mol)
    for c_ring in sites:
        c_new = emol.AddAtom(Chem.Atom(6))
        emol.AddBond(c_ring, c_new, Chem.BondType.SINGLE)
    out = emol.GetMol()
    _safe_sanitize(out)
    out = _fix_heteroatom_bond_types(out)
    return out


# Groups that are not implementable on pure all-aromatic PAH edge carbons.
# They require sp2 C with ≥2 free valence (absent in our skeletons).
_FALLBACK_GROUPS: Dict[str, str] = {
    "carbonyl": "phenolic",
    "quinone":  "phenolic",
    "lactone":  "phenolic",
}


class OxygenAssigner:
    """
    Assign oxygen atoms and functional groups to carbon skeleton.

    Supported groups and oxygen count per placement:
        phenolic  — Ar-OH           (1 O)
        hydroxyl  — Ar-OH (= phenolic for pure PAH)  (1 O)
        carboxyl  — Ar-C(=O)(OH)    (2 O)
        ether     — Ar-O-Ar bridge  (1 O, needs 2 edge sites)
        amino     — Ar-NH2           (0 O, 1 N)
        thiol     — Ar-SH            (0 O, 1 S)
        thioether — Ar-S-Ar bridge   (0 O, 1 S, needs 2 edge sites)
        carbonyl  — not supported for pure aromatic PAH; falls back to phenolic
        quinone   — not supported for pure aromatic PAH; falls back to phenolic
        lactone   — not supported for pure aromatic PAH; falls back to phenolic

    Usage — exact-count dict:
        assigner.assign_oxygens(mol, ...,
            functional_group_preference={"phenolic": 3, "carboxyl": 1})

    Usage — O_C_ratio-driven (default, phenolic only):
        assigner.assign_oxygens(mol, target_O_C_ratio=0.1)
    """

    def __init__(self, seed: int = None, max_ether_span: int = 3):
        if max_ether_span < 3:
            raise ValueError(
                f"max_ether_span must be ≥ 3 (minimum for a 5-membered ring), "
                f"got {max_ether_span}"
            )
        self.seed = seed
        # Maximum C-skeleton shortest-path distance (bonds) between the two
        # ring carbons that an ether oxygen may bridge.  Enforced minimum = 3
        # (5-membered ring); the ring formed has (max_ether_span + 2) members.
        #   3 → 5-membered ring (furan/benzofuran-like) ← safe default
        #   4 → 6-membered ring (pyran/chromene-like)
        #   5 → 7-membered ring
        self._max_ether_span = max_ether_span
        # Instance-level RNG — does not modify the global random state.
        self._rng = random.Random(seed)

    # ------------------------------------------------------------------ #
    #  Public interface                                                    #
    # ------------------------------------------------------------------ #

    # Groups that attach to sp3 (aliphatic) carbons; everything else draws on
    # the aromatic edge pool.
    _ALIPHATIC_GROUPS = frozenset({"aliphatic_hydroxyl"})

    def assign_oxygens(
        self,
        mol: Chem.Mol,
        target_O_C_ratio: float,
        O_C_tolerance: float = 0.10,
        functional_group_preference: Optional[Dict[str, int]] = None,
        allow_aliphatic_oxygen: bool = True,
    ) -> Tuple[Chem.Mol, CompositionResult]:
        """
        Assign oxygen-containing functional groups to the molecule.

        Args:
            mol: RDKit molecule (carbon skeleton, all-aromatic)
            target_O_C_ratio: Target O/C ratio (used only when
                functional_group_preference is None).
            O_C_tolerance: Tolerance for ratio check (informational only).
            functional_group_preference: Dict mapping group name → exact count,
                e.g. {"phenolic": 3, "carboxyl": 1}.
                If None, places phenolic groups to reach target_O_C_ratio.
            allow_aliphatic_oxygen: In O_C_ratio mode, when the aromatic edge
                sites cannot hold the target oxygen, place the remainder as
                aliphatic hydroxyls (-CH2-OH) on the sp3 carbons.  Set False to
                reproduce the aromatic-only behaviour exactly.  Has no effect in
                dict mode, where the requested groups are always honoured, and
                no effect when the skeleton carries no aliphatic carbons.

        Returns:
            (modified_mol, composition_result) where composition_result carries
            ``placed_counts`` (groups actually placed) and ``requested_counts``
            (groups requested after normalisation).
        """
        num_carbons = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)

        # ── Build the two carbon pools ─────────────────────────────────────
        # Aromatic edge carbons take phenolic/carboxyl/ether/etc; sp3 carbons
        # take aliphatic hydroxyls.  Both are shuffled so placement is spread.
        aromatic_pool = self._get_edge_sites(mol)
        aliphatic_pool = self._get_aliphatic_sites(mol)
        self._rng.shuffle(aromatic_pool)
        self._rng.shuffle(aliphatic_pool)

        # ── Mode 1: exact-count dict ───────────────────────────────────────
        if isinstance(functional_group_preference, dict) and functional_group_preference:
            groups_spec = self._validate_and_normalise_spec(functional_group_preference)
        # ── Mode 2: O_C_ratio-driven (legacy / default) ────────────────────
        else:
            target_num_oxygens = int(round(num_carbons * target_O_C_ratio))
            if target_num_oxygens <= 0:
                return mol, self._calculate_composition(mol, {}, {}, {})
            # Fill aromatic edges with phenolic first; spill the shortfall onto
            # sp3 carbons as aliphatic hydroxyls.  Without the spill, a
            # low-aromaticity char whose edges are mostly consumed by aliphatic
            # decoration + H-saturation cannot reach its oxygen target at all.
            n_phenolic = min(target_num_oxygens, len(aromatic_pool))
            groups_spec = {}
            if n_phenolic > 0:
                groups_spec["phenolic"] = n_phenolic
            if allow_aliphatic_oxygen:
                n_aliphatic = min(target_num_oxygens - n_phenolic, len(aliphatic_pool))
                if n_aliphatic > 0:
                    groups_spec["aliphatic_hydroxyl"] = n_aliphatic

        if not groups_spec:
            return mol, self._calculate_composition(mol, {}, {}, {})

        requested_counts: Dict[str, int] = dict(groups_spec)
        used_sites: Set[int] = set()

        # ── Warn if requested counts exceed a feasible estimate ───────────────
        if isinstance(functional_group_preference, dict) and functional_group_preference:
            n_arom = len(aromatic_pool)
            n_aliph = len(aliphatic_pool)
            max_estimates = {"ether": n_arom // 2, "thioether": n_arom // 2,
                             "aliphatic_hydroxyl": n_aliph}
            for _g in ("phenolic", "hydroxyl", "carboxyl", "amino", "thiol"):
                max_estimates[_g] = n_arom
            for _g, _cnt in groups_spec.items():
                _max = max_estimates.get(_g, n_arom)
                if _cnt > _max * 1.5:
                    logger.warning(
                        "Requested %d '%s' groups but estimated max for this "
                        "%d-carbon molecule is ~%d. "
                        "Actual placement may be much lower.",
                        _cnt, _g, num_carbons, _max,
                    )

        emol = Chem.EditableMol(mol)
        new_mol = mol
        functional_groups_added: Dict[str, int] = {}

        # ── Place each group type the requested number of times ────────────
        for group, count in groups_spec.items():
            pool = aliphatic_pool if group in self._ALIPHATIC_GROUPS else aromatic_pool
            placed = 0
            while placed < count:
                free = [s for s in pool if s not in used_sites]
                ok, new_mol, emol, n_O, sites_used = self._place_group(
                    group, emol, new_mol, free
                )
                if not ok:
                    logger.warning(
                        "Could not place all '%s' groups (%d/%d placed — not enough sites)",
                        group, placed, count,
                    )
                    break
                placed += 1
                used_sites.update(sites_used)
                functional_groups_added[group] = (
                    functional_groups_added.get(group, 0) + 1
                )

        _safe_sanitize(new_mol)
        # Fix any heteroatom bonds spuriously marked AROMATIC (e.g., ether O
        # closing a ring that satisfies Hückel's rule through the PAH system).
        new_mol = _fix_heteroatom_bond_types(new_mol)
        return new_mol, self._calculate_composition(
            new_mol, functional_groups_added, functional_groups_added, requested_counts
        )

    # ------------------------------------------------------------------ #
    #  Private helpers                                                    #
    # ------------------------------------------------------------------ #

    def _validate_and_normalise_spec(
        self, spec: Dict[str, int]
    ) -> Dict[str, int]:
        """
        Validate group names and substitute unsupported ones with their fallback.
        Returns the cleaned spec dict.
        """
        valid_keys = set(FUNCTIONAL_GROUPS.keys())

        # Remove unknown keys
        unknown = [k for k in spec if k not in valid_keys]
        if unknown:
            logger.warning("Unknown functional groups ignored: %s", unknown)
        cleaned: Dict[str, int] = {k: v for k, v in spec.items() if k in valid_keys}

        # Substitute unsupported groups (warn once each)
        warned: Set[str] = set()
        for g in list(cleaned):
            if g in _FALLBACK_GROUPS:
                fallback = _FALLBACK_GROUPS[g]
                if g not in warned:
                    logger.warning(
                        "'%s' is not supported for pure aromatic PAH "
                        "(needs ≥2 free valence on one carbon); "
                        "substituting with '%s'",
                        g, fallback,
                    )
                    warned.add(g)
                cleaned[fallback] = cleaned.get(fallback, 0) + cleaned.pop(g)

        return cleaned

    def _get_edge_sites(self, mol: Chem.Mol) -> List[int]:
        """Return indices of edge aromatic carbons with ≥1 free valence."""
        return [
            a.GetIdx()
            for a in mol.GetAtoms()
            if a.GetAtomicNum() == 6
            and a.GetIsAromatic()
            and ValenceValidator.get_valence_info(mol, a.GetIdx()).available_valence >= 1
        ]

    def _get_aliphatic_sites(self, mol: Chem.Mol) -> List[int]:
        """Return indices of sp3 (non-aromatic) carbons with ≥1 free valence.

        These are the pendant methyls/methylenes the H/C-shaping stage adds.
        Each can host one hydroxyl (a -CH3 becomes -CH2-OH); a second -OH on the
        same carbon would be a geminal diol, so the caller consumes each site
        once, exactly as it does for aromatic edge sites.
        """
        return [
            a.GetIdx()
            for a in mol.GetAtoms()
            if a.GetAtomicNum() == 6
            and not a.GetIsAromatic()
            and ValenceValidator.get_valence_info(mol, a.GetIdx()).available_valence >= 1
        ]

    def _place_group(
        self,
        group: str,
        emol: Chem.EditableMol,
        new_mol: Chem.Mol,
        free_sites: List[int],
    ) -> Tuple[bool, Chem.Mol, Chem.EditableMol, int, Set[int]]:
        """
        Attempt to attach one instance of *group* to the molecule.

        Returns:
            (success, updated_mol, updated_emol, oxygens_added, sites_consumed)
        """
        BT = Chem.BondType

        if group in ("phenolic", "hydroxyl", "aliphatic_hydroxyl"):
            # C-OH: single-bond O then H.  Identical construction whether the
            # carbon is an aromatic edge (phenolic) or an sp3 pendant carbon
            # (aliphatic_hydroxyl); the caller chooses the pool of carbons.
            if not free_sites:
                return False, new_mol, emol, 0, set()
            c = free_sites[0]
            o = emol.AddAtom(Chem.Atom(8))
            h = emol.AddAtom(Chem.Atom(1))
            emol.AddBond(c, o, BT.SINGLE)
            emol.AddBond(o, h, BT.SINGLE)
            return True, emol.GetMol(), emol, 1, {c}

        elif group == "carboxyl":
            # Ar-C(=O)(OH): one ring C → new sp2 C → =O and -OH
            if not free_sites:
                return False, new_mol, emol, 0, set()
            c_ring = free_sites[0]
            c_new  = emol.AddAtom(Chem.Atom(6))  # carboxyl carbon
            o1     = emol.AddAtom(Chem.Atom(8))   # C=O
            o2     = emol.AddAtom(Chem.Atom(8))   # C-OH
            h      = emol.AddAtom(Chem.Atom(1))
            emol.AddBond(c_ring, c_new, BT.SINGLE)
            emol.AddBond(c_new,  o1,   BT.DOUBLE)
            emol.AddBond(c_new,  o2,   BT.SINGLE)
            emol.AddBond(o2,     h,    BT.SINGLE)
            return True, emol.GetMol(), emol, 2, {c_ring}

        elif group == "ether":
            # Ar-O-Ar bridge: two edge C atoms connected through one O.
            #
            # Ring size formed = span + 2  (span intermediate C atoms + the
            # two bridged C atoms + the O).
            #   span 1 → 3-membered ring (epoxide-like, forbidden)
            #   span 2 → 4-membered ring (too strained, forbidden)
            #   span 3 → 5-membered ring (furan/benzofuran-like — minimum)
            #   span 4 → 6-membered ring (pyran/chromene-like)
            #   span 5 → 7-membered ring (default maximum)
            #
            # Large spans fold the flat sheet into a tube; enforce max_ether_span.
            max_span = getattr(self, "_max_ether_span", 5)
            c1, c2 = None, None
            for i in range(len(free_sites)):
                for j in range(i + 1, len(free_sites)):
                    s1, s2 = free_sites[i], free_sites[j]
                    # Shortest C-skeleton path must be within the allowed span
                    path = Chem.GetShortestPath(new_mol, s1, s2)
                    path_len = len(path) - 1  # number of bonds
                    # Minimum 3: span < 3 gives strained 3- or 4-membered rings
                    # that are geometrically incompatible with the flat hex lattice.
                    if path_len < 3 or path_len > max_span:
                        continue
                    c1, c2 = s1, s2
                    break
                if c1 is not None:
                    break
            if c1 is None:
                return False, new_mol, emol, 0, set()
            o = emol.AddAtom(Chem.Atom(8))
            emol.AddBond(c1, o, BT.SINGLE)
            emol.AddBond(c2, o, BT.SINGLE)
            return True, emol.GetMol(), emol, 1, {c1, c2}

        elif group == "amino":
            # Ar-NH2: single-bond N then two H atoms
            if not free_sites:
                return False, new_mol, emol, 0, set()
            c = free_sites[0]
            n = emol.AddAtom(Chem.Atom(7))
            h1 = emol.AddAtom(Chem.Atom(1))
            h2 = emol.AddAtom(Chem.Atom(1))
            emol.AddBond(c, n, BT.SINGLE)
            emol.AddBond(n, h1, BT.SINGLE)
            emol.AddBond(n, h2, BT.SINGLE)
            return True, emol.GetMol(), emol, 0, {c}

        elif group == "thiol":
            # Ar-SH: single-bond S then one H atom
            if not free_sites:
                return False, new_mol, emol, 0, set()
            c = free_sites[0]
            s = emol.AddAtom(Chem.Atom(16))
            h = emol.AddAtom(Chem.Atom(1))
            emol.AddBond(c, s, BT.SINGLE)
            emol.AddBond(s, h, BT.SINGLE)
            return True, emol.GetMol(), emol, 0, {c}

        elif group == "thioether":
            # Ar-S-Ar bridge: two edge C atoms connected through one S.
            # Mirrors the ether logic; span constraints keep the bridged ring
            # within a geometrically reasonable 5- to 7-membered size.
            max_span = getattr(self, "_max_ether_span", 5)
            c1, c2 = None, None
            for i in range(len(free_sites)):
                for j in range(i + 1, len(free_sites)):
                    s1, s2 = free_sites[i], free_sites[j]
                    path = Chem.GetShortestPath(new_mol, s1, s2)
                    path_len = len(path) - 1  # number of bonds
                    if path_len < 3 or path_len > max_span:
                        continue
                    c1, c2 = s1, s2
                    break
                if c1 is not None:
                    break
            if c1 is None:
                return False, new_mol, emol, 0, set()
            s = emol.AddAtom(Chem.Atom(16))
            emol.AddBond(c1, s, BT.SINGLE)
            emol.AddBond(c2, s, BT.SINGLE)
            return True, emol.GetMol(), emol, 0, {c1, c2}

        # Unknown / unsupported (should have been filtered before reaching here)
        return False, new_mol, emol, 0, set()

    def _find_attachment_sites(self, mol: Chem.Mol) -> List[int]:
        """Find carbon atoms suitable for functional group attachment."""
        sites = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                if atom.GetTotalValence() < 4:
                    sites.append(atom.GetIdx())
        return sites

    def _calculate_composition(
        self,
        mol: Chem.Mol,
        functional_groups: Dict[str, int],
        placed_counts: Dict[str, int],
        requested_counts: Dict[str, int],
    ) -> CompositionResult:
        """Calculate composition information from molecule."""
        num_C = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        num_H = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
        num_O = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
        num_N = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
        num_S = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)

        H_C_ratio = num_H / max(num_C, 1)
        O_C_ratio = num_O / max(num_C, 1)
        N_C_ratio = num_N / max(num_C, 1)
        S_C_ratio = num_S / max(num_C, 1)

        return CompositionResult(
            num_carbons=num_C,
            num_hydrogens=num_H,
            num_oxygens=num_O,
            num_nitrogens=num_N,
            num_sulfurs=num_S,
            H_C_ratio=H_C_ratio,
            O_C_ratio=O_C_ratio,
            N_C_ratio=N_C_ratio,
            S_C_ratio=S_C_ratio,
            functional_groups=functional_groups,
            placed_counts=placed_counts,
            requested_counts=requested_counts,
        )


class NitrogenSubstitutor:
    """
    Substitute ring carbon atoms with nitrogen (biochar N-doping).

    Unlike :class:`OxygenAssigner`, which *adds* pendant heteroatoms to the
    skeleton, this class *replaces* selected skeleton carbon atoms with
    nitrogen.  It is run after oxygen assignment but before hydrogen
    saturation, so that :class:`HydrogenAssigner` can satisfy the new valence
    requirements correctly.

    Three substitution modes (matching the chemistry of N-doped carbons):

        pyridinic  — edge 6-ring C with exactly one H is swapped for N; the H
                     is removed.  Pyridine-like N contributes a lone pair to the
                     sp2 plane and carries no H (2 ring bonds + lone pair).
        pyrrolic   — a 5-ring C (only present when defect_fraction > 0 created
                     pentagons) is swapped for N and keeps/gains an H, like the
                     N-H of pyrrole.
        graphitic  — an interior ring C (degree 3, no H) is swapped for N
                     (quaternary / "graphitic" N).

    The substitutor is robust: if fewer suitable sites exist than requested,
    it places as many as possible, logs a warning, and records the actual
    counts placed.
    """

    def __init__(self, seed: Optional[int] = None):
        self.seed = seed
        self._rng = random.Random(seed)
        # Number of each N type actually placed in the most recent call.
        self.placed_pyridinic = 0
        self.placed_pyrrolic = 0
        self.placed_graphitic = 0

    def substitute(
        self,
        mol: Chem.Mol,
        n_pyridinic: int = 0,
        n_pyrrolic: int = 0,
        n_graphitic: int = 0,
    ) -> Chem.Mol:
        """
        Replace ring carbons with nitrogen.

        Args:
            mol: RDKit molecule (carbon skeleton, possibly with pendant O groups,
                before hydrogen saturation).
            n_pyridinic: Number of pyridinic (edge 6-ring, no-H) N to place.
            n_pyrrolic: Number of pyrrolic (5-ring, N-H) N to place.
            n_graphitic: Number of graphitic (interior, no-H) N to place.

        Returns:
            Modified RDKit molecule.  The actual numbers placed are stored on
            ``self.placed_pyridinic`` / ``placed_pyrrolic`` / ``placed_graphitic``.
        """
        self.placed_pyridinic = 0
        self.placed_pyrrolic = 0
        self.placed_graphitic = 0

        if n_pyridinic <= 0 and n_pyrrolic <= 0 and n_graphitic <= 0:
            return mol

        used: Set[int] = set()

        # Graphitic first: interior atoms are scarcest and most constrained.
        if n_graphitic > 0:
            mol, placed = self._substitute_graphitic(mol, n_graphitic, used)
            self.placed_graphitic = placed
            if placed < n_graphitic:
                logger.warning(
                    "Could only place %d/%d graphitic N (not enough interior "
                    "ring carbons)", placed, n_graphitic,
                )

        if n_pyrrolic > 0:
            mol, placed = self._substitute_pyrrolic(mol, n_pyrrolic, used)
            self.placed_pyrrolic = placed
            if placed < n_pyrrolic:
                logger.warning(
                    "Could only place %d/%d pyrrolic N (not enough 5-ring "
                    "carbons — set defect_fraction > 0 to create pentagons)",
                    placed, n_pyrrolic,
                )

        if n_pyridinic > 0:
            mol, placed = self._substitute_pyridinic(mol, n_pyridinic, used)
            self.placed_pyridinic = placed
            if placed < n_pyridinic:
                logger.warning(
                    "Could only place %d/%d pyridinic N (not enough edge 6-ring "
                    "carbons with a free H)", placed, n_pyridinic,
                )

        # NOTE: ring-substituted N (pyridinic / graphitic) is a genuine
        # aromatic ring member — its in-ring bonds must stay AROMATIC.  We must
        # NOT call _fix_heteroatom_bond_types here (that helper is for pendant
        # ether O, which de-aromatises C-O bonds).  Re-perceive aromaticity so
        # the new N atoms are correctly flagged aromatic and their ring bonds
        # remain order 1.5.
        _safe_sanitize(mol)
        return mol

    # ------------------------------------------------------------------ #
    #  Private helpers                                                    #
    # ------------------------------------------------------------------ #

    def _ring_sizes(self, mol: Chem.Mol, idx: int) -> List[int]:
        """Sizes of all rings containing atom *idx*."""
        ring_info = mol.GetRingInfo()
        return [len(r) for r in ring_info.AtomRings() if idx in r]

    def _candidate_carbons(self, mol: Chem.Mol, used: Set[int]) -> List[int]:
        """Aromatic ring carbons not yet substituted, shuffled."""
        cands = [
            a.GetIdx()
            for a in mol.GetAtoms()
            if a.GetAtomicNum() == 6
            and mol.GetRingInfo().NumAtomRings(a.GetIdx()) > 0
            and a.GetIdx() not in used
        ]
        self._rng.shuffle(cands)
        return cands

    def _swap_carbon_to_nitrogen(
        self, mol: Chem.Mol, c_idx: int, remove_h: bool, formal_charge: int = 0
    ) -> Chem.Mol:
        """
        Replace the carbon at *c_idx* with nitrogen, optionally removing one
        bonded hydrogen.  Returns a new mol (atom index of the N is preserved).

        *formal_charge* is applied to the new nitrogen.  Graphitic (quaternary)
        N has three aromatic ring bonds and must carry a +1 formal charge so the
        aromatic π system remains closed-shell and RDKit can kekulize it (this
        is the standard pyridinium-like representation of graphitic N).
        """
        rw = Chem.RWMol(mol)
        # Remove a bonded explicit H first (indices shift only for atoms after
        # the removed H, and c_idx < any later-added H in our skeletons; we
        # remove after recording to be safe).
        h_to_remove = None
        if remove_h:
            for nbr in rw.GetAtomWithIdx(c_idx).GetNeighbors():
                if nbr.GetAtomicNum() == 1:
                    h_to_remove = nbr.GetIdx()
                    break
        n_atom = rw.GetAtomWithIdx(c_idx)
        n_atom.SetAtomicNum(7)
        n_atom.SetNoImplicit(False)
        n_atom.SetFormalCharge(formal_charge)
        if h_to_remove is not None:
            rw.RemoveAtom(h_to_remove)
        return rw.GetMol()

    def _substitute_pyridinic(
        self, mol: Chem.Mol, count: int, used: Set[int]
    ) -> Tuple[Chem.Mol, int]:
        """Edge 6-ring carbon with exactly one H → N, drop the H."""
        placed = 0
        for c_idx in self._candidate_carbons(mol, used):
            atom = mol.GetAtomWithIdx(c_idx)
            sizes = self._ring_sizes(mol, c_idx)
            if 6 not in sizes:
                continue
            # Edge carbon: degree 2 ring + (implicit or explicit) H available.
            h_count = atom.GetTotalNumHs() + sum(
                1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 1
            )
            heavy_deg = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() != 1)
            if heavy_deg != 2:
                continue
            mol = self._swap_carbon_to_nitrogen(mol, c_idx, remove_h=(h_count >= 1))
            used.add(c_idx)
            placed += 1
            if placed >= count:
                break
        return mol, placed

    def _substitute_pyrrolic(
        self, mol: Chem.Mol, count: int, used: Set[int]
    ) -> Tuple[Chem.Mol, int]:
        """
        5-ring carbon → N with an explicit N-H (pyrrole-like).

        A pyrrolic N donates its lone pair to the aromatic π system, so it forms
        single σ-bonds to its two ring neighbours plus an N-H (3 σ-bonds total).
        After aromaticity perception its ring bonds become order 1.5 and the
        valence checker would otherwise consider it saturated, so we attach the
        N-H explicitly here (and pin it with SetNoImplicit) rather than relying
        on HydrogenAssigner.
        """
        placed = 0
        for c_idx in self._candidate_carbons(mol, used):
            sizes = self._ring_sizes(mol, c_idx)
            if 5 not in sizes:
                continue
            mol = self._swap_carbon_to_nitrogen(mol, c_idx, remove_h=False)
            # Attach an explicit N-H and forbid implicit Hs on this nitrogen.
            rw = Chem.RWMol(mol)
            n_atom = rw.GetAtomWithIdx(c_idx)
            h_idx = rw.AddAtom(Chem.Atom(1))
            rw.AddBond(c_idx, h_idx, Chem.BondType.SINGLE)
            n_atom.SetNumExplicitHs(0)
            n_atom.SetNoImplicit(True)
            mol = rw.GetMol()
            used.add(c_idx)
            placed += 1
            if placed >= count:
                break
        return mol, placed

    def _substitute_graphitic(
        self, mol: Chem.Mol, count: int, used: Set[int]
    ) -> Tuple[Chem.Mol, int]:
        """Interior ring carbon (3 heavy neighbours, no H) → N."""
        placed = 0
        for c_idx in self._candidate_carbons(mol, used):
            atom = mol.GetAtomWithIdx(c_idx)
            sizes = self._ring_sizes(mol, c_idx)
            if 6 not in sizes:
                continue
            heavy_deg = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() != 1)
            h_count = atom.GetTotalNumHs() + sum(
                1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 1
            )
            # Interior junction: 3 heavy ring neighbours, no hydrogen.
            if heavy_deg != 3 or h_count != 0:
                continue
            # Graphitic N is pyridinium-like → formal +1 keeps the ring
            # closed-shell and kekulizable.
            mol = self._swap_carbon_to_nitrogen(
                mol, c_idx, remove_h=False, formal_charge=1
            )
            used.add(c_idx)
            placed += 1
            if placed >= count:
                break
        return mol, placed


class HydrogenAssigner:
    """
    Assign hydrogen atoms to satisfy valence and H/C ratio.

    Strategy:
    1. First, ensure all atoms have minimum valence satisfied
    2. Then, adjust hydrogen count to match target H/C ratio
    """

    def __init__(self, seed: int = None):
        self.seed = seed
        self._rng = random.Random(seed)

    def assign_hydrogens(
        self,
        mol: Chem.Mol,
        target_H_C_ratio: float,
        H_C_tolerance: float = 0.10,
        result: Optional[CompositionResult] = None,
    ) -> Tuple[Chem.Mol, CompositionResult]:
        """
        Assign hydrogens to molecule to satisfy minimum valence and target H/C ratio.

        Strategy:
        1. Add hydrogens to satisfy minimum valences for all atoms
        2. Trim excess hydrogens to match target H/C ratio (if needed)
        3. Ensure no violations of valence constraints

        Args:
            mol: RDKit molecule (may already have oxygens)
            target_H_C_ratio: Target hydrogen-to-carbon ratio
            H_C_tolerance: Tolerance for ratio
            result: Existing :class:`CompositionResult` from the oxygen step.
                If provided, its hydrogen-related fields are updated in-place and
                the same object is returned, preserving ``functional_groups``,
                ``placed_counts``, and ``requested_counts`` set earlier.
                If ``None``, a new :class:`CompositionResult` is created.

        Returns:
            (modified_mol_with_H, composition_result)
        """
        num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        target_num_hydrogens = int(round(num_carbons * target_H_C_ratio))

        # Step 1: Ensure minimum valences are satisfied
        mol_saturated = self._saturate_valences(mol)

        # Step 2: Adjust hydrogen count to match target ratio
        num_hydrogens = sum(1 for atom in mol_saturated.GetAtoms() if atom.GetAtomicNum() == 1)

        if num_hydrogens > target_num_hydrogens:
            mol_trimmed = self._trim_hydrogens(
                mol_saturated,
                num_hydrogens - target_num_hydrogens
            )
            mol_with_H = mol_trimmed
        else:
            # The saturated skeleton already has at most the target number of
            # hydrogens: every free edge valence is capped and there is no site
            # to add more without changing the carbon topology.  When the target
            # is strictly above this ceiling the request is structurally
            # unreachable -- record it and warn honestly instead of silently
            # returning a hydrogen-deficient structure that only fails later in
            # composition validation.
            mol_with_H = mol_saturated
            if num_hydrogens < target_num_hydrogens and num_carbons > 0:
                ceiling = num_hydrogens / num_carbons
                if result is not None:
                    result.h_c_ceiling = ceiling
                    result.h_c_target_unreachable = True
                logger.warning(
                    "Requested H/C %.3f exceeds the structural ceiling %.3f for "
                    "this %d-carbon skeleton (max %d H when all edge valences are "
                    "saturated; %d requested). The structure was capped at the "
                    "ceiling. Raising H/C requires a less-condensed skeleton or "
                    "aliphatic/sp3 carbons, not more edge hydrogens.",
                    target_H_C_ratio, ceiling, num_carbons,
                    num_hydrogens, target_num_hydrogens,
                )

        _safe_sanitize(mol_with_H)
        # After aromaticity re-perception, ether oxygens that bridge two ring
        # carbons can be spuriously marked AROMATIC.  Fix them back to SINGLE.
        mol_with_H = _fix_heteroatom_bond_types(mol_with_H)

        composition = self._update_composition(mol_with_H, result)
        return mol_with_H, composition

    def _saturate_valences(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Add hydrogens to ensure all atoms meet minimum valence requirements.

        Uses valence validation to ensure each atom gets enough bonds.
        """
        emol = Chem.EditableMol(mol)

        atoms_to_saturate = []
        for atom in mol.GetAtoms():
            val_info = ValenceValidator.get_valence_info(mol, atom.GetIdx())
            if val_info.needed_valence > 0:
                atoms_to_saturate.append((atom.GetIdx(), val_info.needed_valence))

        for atom_idx, needed_bonds in atoms_to_saturate:
            for _ in range(needed_bonds):
                h_idx = emol.AddAtom(Chem.Atom(1))
                emol.AddBond(atom_idx, h_idx, Chem.BondType.SINGLE)

        return emol.GetMol()

    def _trim_hydrogens(self, mol: Chem.Mol, num_to_remove: int) -> Chem.Mol:
        """
        Remove excess hydrogens while preserving minimum valence constraints.

        Only removes H atoms from atoms that still have available valence.
        """
        if num_to_remove <= 0:
            return mol

        emol = Chem.EditableMol(mol)

        removable_h = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                neighbors = atom.GetNeighbors()
                if neighbors:
                    parent = neighbors[0]
                    parent_val = ValenceValidator.get_valence_info(mol, parent.GetIdx())
                    if parent_val.available_valence > 0:
                        removable_h.append(atom.GetIdx())

        for h_idx in sorted(removable_h, reverse=True)[:num_to_remove]:
            emol.RemoveAtom(h_idx)

        return emol.GetMol()

    def _update_composition(
        self, mol: Chem.Mol, result: Optional[CompositionResult] = None
    ) -> CompositionResult:
        """
        Update atom-count fields on *result* in-place, or create a new one.

        Preserves ``functional_groups``, ``placed_counts``, and
        ``requested_counts`` from the oxygen step when *result* is provided.
        """
        num_C = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        num_H = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
        num_O = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
        num_N = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
        num_S = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)

        H_C_ratio = num_H / max(num_C, 1)
        O_C_ratio = num_O / max(num_C, 1)
        N_C_ratio = num_N / max(num_C, 1)
        S_C_ratio = num_S / max(num_C, 1)

        if result is not None:
            result.num_carbons = num_C
            result.num_hydrogens = num_H
            result.num_oxygens = num_O
            result.num_nitrogens = num_N
            result.num_sulfurs = num_S
            result.H_C_ratio = H_C_ratio
            result.O_C_ratio = O_C_ratio
            result.N_C_ratio = N_C_ratio
            result.S_C_ratio = S_C_ratio
            return result

        return CompositionResult(
            num_carbons=num_C,
            num_hydrogens=num_H,
            num_oxygens=num_O,
            num_nitrogens=num_N,
            num_sulfurs=num_S,
            H_C_ratio=H_C_ratio,
            O_C_ratio=O_C_ratio,
            N_C_ratio=N_C_ratio,
            S_C_ratio=S_C_ratio,
            functional_groups={},
        )


class HeteroatomValidator:
    """Validate heteroatom assignments."""

    @staticmethod
    def validate_ratios(
        composition: CompositionInfo,
        target_H_C: float,
        target_O_C: float,
        H_C_tolerance: float = 0.10,
        O_C_tolerance: float = 0.10,
    ) -> Tuple[bool, List[str]]:
        """
        Validate that ratios are within tolerance.

        Args:
            composition: CompositionInfo object
            target_H_C: Target H/C ratio
            target_O_C: Target O/C ratio
            H_C_tolerance: Tolerance for H/C
            O_C_tolerance: Tolerance for O/C

        Returns:
            (is_valid, error_messages)
        """
        errors = []

        H_C_error = abs(composition.H_C_ratio - target_H_C) / max(target_H_C, 0.01)
        if H_C_error > H_C_tolerance:
            errors.append(
                f"H/C ratio {composition.H_C_ratio:.3f} outside tolerance "
                f"(target: {target_H_C:.3f}, tolerance: {H_C_tolerance:.1%})"
            )

        if target_O_C > 0:
            O_C_error = abs(composition.O_C_ratio - target_O_C) / max(target_O_C, 0.01)
            if O_C_error > O_C_tolerance:
                errors.append(
                    f"O/C ratio {composition.O_C_ratio:.3f} outside tolerance "
                    f"(target: {target_O_C:.3f}, tolerance: {O_C_tolerance:.1%})"
                )

        return len(errors) == 0, errors

    @staticmethod
    def validate_functional_groups(
        composition: CompositionInfo, expected_groups: Optional[List[str]] = None
    ) -> Tuple[bool, List[str]]:
        """
        Validate functional group assignment.

        Args:
            composition: CompositionInfo object
            expected_groups: List of expected functional groups

        Returns:
            (is_valid, error_messages)
        """
        errors = []

        if expected_groups is not None:
            for group in expected_groups:
                if group not in composition.functional_groups:
                    errors.append(f"Expected functional group '{group}' not found")

        return len(errors) == 0, errors
