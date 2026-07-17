"""
pH-Dependent Protonation

Decide, per ionizable site, whether that site is protonated or ionized at a
given pH, and apply the decision to the molecule.

Scope: discrete protonation states assigned **once, at build time**, sampled
from Henderson-Hasselbalch. This is not constant-pH MD -- there is no
lambda-dynamics and no state change during a run. A generated structure is one
sample from the ensemble at that pH, which is why pH near a pKa needs replicates
rather than a single structure (see :class:`ProtonationAssigner`).
"""

import logging
import random
from typing import Dict, List, Optional, Tuple

from rdkit import Chem

from .constants import PH_MAX, PH_MIN, PROTONATION_STATES
from .heteroatom_assignment import CompositionResult, _safe_sanitize
from .opls_typing import AtomTyper

logger = logging.getLogger(__name__)


def fraction_ionized(pKa: float, pH: float, kind: str) -> float:
    """
    Fraction of sites in the **ionized** form, from Henderson-Hasselbalch.

    ``pKa`` always belongs to the protonated member of the pair, so the fraction
    deprotonated is ``1 / (1 + 10 ** (pKa - pH))``. What counts as *ionized*
    then depends on which side the charge sits:

        acidic — the neutral form is protonated and loses H+ to become an anion,
                 so ionized == deprotonated, and the fraction RISES with pH.
        basic  — the neutral form is deprotonated and gains H+ to become a
                 cation, so ionized == protonated, and the fraction FALLS with pH.

    Collapsing those two cases is the single easiest way to get this feature
    silently wrong: it inverts the charge on every nitrogen while still
    producing plausible-looking output.

    Args:
        pKa: Acid dissociation constant of the protonated form.
        pH: Environmental pH.
        kind: "acidic" or "basic".

    Returns:
        Probability in [0, 1] that a given site is ionized.
    """
    if kind not in ("acidic", "basic"):
        raise ValueError(
            f"kind must be 'acidic' or 'basic', got {kind!r}"
        )
    # Guard the exponent: at |pKa - pH| > ~300 the power overflows, and the
    # answer is saturated long before that anyway.
    exponent = max(-300.0, min(300.0, pKa - pH))
    f_deprotonated = 1.0 / (1.0 + 10.0 ** exponent)
    return f_deprotonated if kind == "acidic" else 1.0 - f_deprotonated


class ProtonationAssigner:
    """
    Set each ionizable site's protonation state from the environmental pH.

    Each site is drawn **independently** against its Henderson-Hasselbalch
    probability, using an instance-local RNG seeded for reproducibility. A
    deterministic ``pH > pKa`` threshold would flip every carboxyl on every
    flake at the same pH, producing a step-function titration curve and no
    site-to-site variation; sampling makes a sweep across seeds a real ensemble.

    The consequence is that a single small structure is one *sample*, not the
    ensemble mean. Near a pKa the outcome is close to a coin flip per site, so
    interpreting one structure as representative is a mistake --
    :attr:`CompositionResult.ionized_counts` reports what was actually placed so
    the sampling stays visible.

    Runs after all heteroatoms exist (so there are sites to titrate) and before
    :class:`HydrogenAssigner` (which owns acidic-H placement and would otherwise
    fight the decision).
    """

    def __init__(self, seed: Optional[int] = None):
        self.seed = seed
        # Instance-level RNG — does not modify the global random state.
        self._rng = random.Random(seed)
        self._typer = AtomTyper()
        # Reverse index: OPLS type of the neutral form -> group name.
        self._by_neutral_type: Dict[str, str] = {
            state.neutral_type: group
            for group, state in PROTONATION_STATES.items()
        }

    # ------------------------------------------------------------------ #
    #  Public interface                                                    #
    # ------------------------------------------------------------------ #

    def assign(
        self,
        mol: Chem.Mol,
        pH: float,
        result: Optional[CompositionResult] = None,
    ) -> Tuple[Chem.Mol, CompositionResult]:
        """
        Apply pH-dependent protonation to *mol*.

        Args:
            mol: RDKit molecule with heteroatoms already placed.
            pH: Environmental pH, within [PH_MIN, PH_MAX].
            result: Existing :class:`CompositionResult` to update in place. A
                new one is created when omitted.

        Returns:
            (modified_mol, composition_result) — the result carries
            ``net_charge``, ``ionized_counts``, and ``titratable_counts``.
        """
        if not PH_MIN <= pH <= PH_MAX:
            raise ValueError(
                f"pH {pH} is outside the aqueous range [{PH_MIN}, {PH_MAX}]. "
                f"Henderson-Hasselbalch saturates well inside these bounds, so "
                f"a value beyond them is far more likely a unit error."
            )

        sites = self._find_sites(mol)
        decisions = self._decide(sites, pH)
        out = self._apply(mol, decisions)

        return out, self._record(out, sites, decisions, result)

    # ------------------------------------------------------------------ #
    #  Private helpers                                                    #
    # ------------------------------------------------------------------ #

    def _find_sites(self, mol: Chem.Mol) -> List[Tuple[int, str]]:
        """
        Locate titratable sites as (atom_idx, group_name).

        Sites are identified through :class:`AtomTyper`, so "what is a carboxyl"
        is answered in exactly one place rather than re-derived here. Anything
        the typer does not recognise as a titratable neutral form -- ether O,
        thioether S, graphitic N -- is simply not a site.
        """
        types = self._typer.assign_atom_types(mol)
        sites: List[Tuple[int, str]] = []
        for idx, opls_type in types.items():
            group = self._by_neutral_type.get(opls_type)
            if group is not None:
                sites.append((idx, group))
        return sites

    def _decide(
        self, sites: List[Tuple[int, str]], pH: float
    ) -> List[Tuple[int, str]]:
        """Draw each site independently; return those that ionize."""
        ionizing: List[Tuple[int, str]] = []
        for idx, group in sites:
            state = PROTONATION_STATES[group]
            f = fraction_ionized(state.pKa, pH, state.kind)
            if self._rng.random() < f:
                ionizing.append((idx, group))
        return ionizing

    def _apply(
        self, mol: Chem.Mol, decisions: List[Tuple[int, str]]
    ) -> Chem.Mol:
        """
        Apply the ionization decisions.

        Ordering matters. Added atoms are appended, so they never disturb
        existing indices; removed atoms shift every later index, so removals are
        deferred and then applied in descending order. Doing the removals inline
        would invalidate the indices of every site not yet visited.
        """
        if not decisions:
            return mol

        rw = Chem.RWMol(mol)
        h_to_remove: List[int] = []

        for idx, group in decisions:
            state = PROTONATION_STATES[group]
            atom = rw.GetAtomWithIdx(idx)

            if state.kind == "acidic":
                h_idx = next(
                    (
                        n.GetIdx()
                        for n in atom.GetNeighbors()
                        if n.GetAtomicNum() == 1
                    ),
                    None,
                )
                if h_idx is None:
                    # No acidic H to give up; leave the site neutral rather than
                    # inventing a charge.
                    logger.debug(
                        "Site %d (%s) has no hydrogen to remove; left neutral",
                        idx, group,
                    )
                    continue
                h_to_remove.append(h_idx)
                atom.SetFormalCharge(-1)
                atom.SetNumExplicitHs(0)
                atom.SetNoImplicit(True)
            else:
                h_idx = rw.AddAtom(Chem.Atom(1))
                rw.AddBond(idx, h_idx, Chem.BondType.SINGLE)
                atom.SetFormalCharge(1)
                atom.SetNumExplicitHs(0)
                atom.SetNoImplicit(True)

        for h_idx in sorted(h_to_remove, reverse=True):
            rw.RemoveAtom(h_idx)

        out = rw.GetMol()
        _safe_sanitize(out)
        return out

    def _record(
        self,
        mol: Chem.Mol,
        sites: List[Tuple[int, str]],
        decisions: List[Tuple[int, str]],
        result: Optional[CompositionResult],
    ) -> CompositionResult:
        """Write the census and net charge onto the composition record."""
        titratable: Dict[str, int] = {}
        for _, group in sites:
            titratable[group] = titratable.get(group, 0) + 1

        ionized: Dict[str, int] = {}
        for _, group in decisions:
            ionized[group] = ionized.get(group, 0) + 1

        net_charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())

        if result is None:
            result = CompositionResult()
        result.net_charge = net_charge
        result.ionized_counts = ionized
        result.titratable_counts = titratable
        return result
