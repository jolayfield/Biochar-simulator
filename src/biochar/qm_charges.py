"""
LigParGen-style QM partial charges (1.14*CM1A) via an external MOPAC binary.

This module reproduces the charge methodology of the LigParGen server
(Dodda et al., *Nucleic Acids Res.* 45, W331 (2017)) as an opt-in alternative
to the static OPLS lookup table (:mod:`biochar.opls_typing`) and the ML refiner
(:mod:`biochar.ml_charges`).

Pipeline
--------
1. **AM1** semiempirical calculation on the generated 3D structure, run with an
   external **MOPAC** binary (OpenMOPAC >= 22).
2. **CM1A mapping** of the AM1 net (Mulliken) charges using the atomic bond
   orders — the "Class IV" charge model of Storer, Giesen, Cramer & Truhlar,
   *J. Comput.-Aided Mol. Des.* 9, 87 (1995).  The exact AM1 parameter set and
   the mapping algorithm are transcribed from the open-source **AMSOL 7.1**
   reference implementation (``chgmp1.f`` / ``include/PARAMS.i``), which is the
   same model BOSS/LigParGen use.
3. **1.14 scaling** — LigParGen's empirical condensed-phase factor for neutral
   molecules — followed by a tiny residual correction so the total charge is
   exactly the molecular formal charge (MOPAC prints charges/bond orders to
   finite precision).

MOPAC is a compiled Fortran program and is **not** pip-installable.  Install it
with ``conda install -c conda-forge mopac`` or from
https://github.com/openmopac/mopac.  The ``mopac`` executable is discovered on
``PATH`` at run time; override with the ``mopac_cmd`` argument or the
``$MOPAC_CMD`` environment variable.

Fidelity notes
--------------
* LigParGen optimises the geometry at the AM1 level before computing charges.
  By default this module instead runs a **single-point** AM1 calculation on the
  geometry produced by the biochar generator (``optimize=False``).  This keeps
  the charges consistent with the exact coordinates written to the GROMACS files,
  is far faster for large graphitic sheets, and avoids AM1 re-optimisation
  distorting a deliberately-built structure.  Pass ``optimize=True`` for strict
  LigParGen-style AM1 geometry optimisation (recommended only for small
  fragments).
* Charges come from MOPAC's AM1, whereas LigParGen uses BOSS's AM1.  The two are
  the same semiempirical Hamiltonian but independent codes, so values agree
  closely rather than bit-for-bit.
* The base model implemented here is **1.14*CM1A** (LigParGen's default).  The
  newer 1.14*CM1A-LBCC bond-charge corrections are not applied; see the module
  ``_LBCC`` note for how they would slot in.
"""

from __future__ import annotations

import logging
import math
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# CM1A (AM1) parameters — transcribed verbatim from AMSOL 7.1 include/PARAMS.i
# Keys are atomic numbers Z.  Absent keys default to 0.0 (no correction).
# ---------------------------------------------------------------------------

# SCALA1: per-element scale factor applied to the Mulliken charge (QC0).
_SCALA1: Dict[int, float] = {
    7: 0.3846,    # N
    9: 0.1468,    # F
    15: 0.0,      # P
    16: -0.1311,  # S
    17: 0.0405,   # Cl
    35: 0.1761,   # Br
    53: 0.2380,   # I
}

# OFSTA1: per-element additive offset for non-hydrogen atoms (QD0).
_OFSTA1: Dict[int, float] = {
    8: -0.0283,   # O
    9: 0.0399,    # F
    16: -0.0956,  # S
    17: -0.0276,  # Cl
    35: -0.0802,  # Br
    53: -0.1819,  # I
}

# HOFFA1: offset applied to a hydrogen atom, indexed by the heavy atom it is
# bonded to (the "partner").  Only N, O, Si are non-zero in base CM1A.
_HOFFA1: Dict[int, float] = {
    7: 0.0850,    # H on N
    8: 0.1447,    # H on O
    14: 0.0640,   # H on Si
}

# OXOFF1: offset applied to an oxygen atom, indexed by its bonded partner.
_OXOFF1: Dict[int, float] = {
    16: -0.0600,  # O on S
}

# ANTOF1: offset applied to a nitrogen atom, indexed by its bonded partner.
# Z=6 (carbon) is used only inside the special N–C switching term below; all
# other partners enter through the ordinary offset path.
_ANTOF1: Dict[int, float] = {
    6: -0.0880,   # N–C  (special tanh-switched term)
    8: -0.0630,   # N–O
}

# LigParGen's empirical scaling factor for neutral molecules.
_CM1A_SCALE = 1.14

# Localized bond-charge corrections (1.14*CM1A-LBCC, Vilseck et al. 2017) are
# intentionally NOT applied here — LigParGen's default export is plain
# 1.14*CM1A.  A future ``lbcc=True`` option would add per-bond-type increments
# after scaling; the 19 published corrections are recorded in the project notes.


class QMChargeError(RuntimeError):
    """Raised when the MOPAC backend is unavailable or fails."""


# ---------------------------------------------------------------------------
# Pure CM1A mapping (no external dependencies — unit-testable in isolation)
# ---------------------------------------------------------------------------

def cm1a_from_am1(
    mulliken_charges: List[float],
    bond_orders: Dict[Tuple[int, int], float],
    atomic_numbers: List[int],
) -> List[float]:
    """Map AM1 net (Mulliken) charges to CM1A charges.

    Implements the CM1 mapping of ``AMSOL/chgmp1.f`` for the AM1 Hamiltonian::

        Q(i) = QM(i) + Σ_j  B(i,j) · DELTA(j)
        DELTA(i) = QC(i)·QM(i) + QD(i)

    where ``B`` is the bond-order matrix in AMSOL's sign convention
    (``B(i,j) = -B_ij`` off-diagonal, ``B(i,i) = Σ_j B_ij``), so the result is
    algebraically charge-conserving: ``Σ_i Q(i) == Σ_i QM(i)``.

    Args:
        mulliken_charges: AM1 net atomic charges (MOPAC "NET ATOMIC CHARGES"),
            one per atom.  These equal AMSOL's ``QM`` (core charge minus the
            diagonal density population).
        bond_orders: Mayer bond orders keyed by ``(i, j)`` with ``i < j``.
            Pairs not present are treated as zero.
        atomic_numbers: Atomic number per atom (same ordering as the charges).

    Returns:
        Unscaled CM1A net charges, one per atom (apply the 1.14 factor
        separately via :func:`scale_and_neutralize`).
    """
    n = len(mulliken_charges)
    if len(atomic_numbers) != n:
        raise ValueError("atomic_numbers and mulliken_charges length mismatch")

    qm = list(mulliken_charges)
    z = list(atomic_numbers)

    # Bond-order matrix in AMSOL convention.
    b = [[0.0] * n for _ in range(n)]
    for (i, j), bij in bond_orders.items():
        if i == j:
            continue
        b[i][j] = -bij
        b[j][i] = -bij
    for i in range(n):
        b[i][i] = -sum(b[i][j] for j in range(n) if j != i)  # = +Σ_j B_ij

    # FD[partner_j][atom_i]: bond-order-weighted offsets (chgmp1.f labels 10–19).
    fd = [[0.0] * n for _ in range(n)]
    for i in range(n):
        zi = z[i]
        for j in range(n):
            if j == i:
                continue
            zj = z[j]
            if zi == 1:                       # hydrogen
                fd[j][i] = -_HOFFA1.get(zj, 0.0)
            elif zi == 7:                     # nitrogen (N–C handled specially)
                if zj != 6:
                    fd[j][i] = -_ANTOF1.get(zj, 0.0)
            elif zi == 8:                     # oxygen
                fd[j][i] = -_OXOFF1.get(zj, 0.0)
            # other elements: no FD offset in base CM1A

    # QC (scale) and QD (offset) vectors.
    qc = [0.0] * n
    qd = [0.0] * n
    for i in range(n):
        qc[i] = _SCALA1.get(z[i], 0.0)
        acc = _OFSTA1.get(z[i], 0.0)
        for j in range(n):
            acc += b[i][j] * fd[j][i]
        qd[i] = acc

    # Special AM1 N–C term: tanh switch that activates for high-order (nitrile)
    # bonds and vanishes for aromatic/amine N–C.  Mirrors chgmp1.f loops 210/220.
    scal_n = _SCALA1.get(7, 0.0)
    antof_c = _ANTOF1.get(6, 0.0)
    for i in range(n):
        if z[i] != 7:
            continue
        for j in range(n):
            if z[j] != 6:
                continue
            # b[j][i] holds -B_NC; tanh(-x) = -tanh(x) handled by the sign here.
            bhat = 0.5 * (1.0 - math.tanh(10.0 * b[j][i] + 23.0))
            qc[i] -= scal_n * bhat
            qd[i] += antof_c * bhat

    # DELTA and final charges Q = QM + B·DELTA.
    delta = [qc[i] * qm[i] + qd[i] for i in range(n)]
    q = [0.0] * n
    for i in range(n):
        acc = qm[i]
        row = b[i]
        for j in range(n):
            acc += row[j] * delta[j]
        q[i] = acc
    return q


def scale_and_neutralize(
    cm1a_charges: List[float],
    total_charge: float = 0.0,
    scale: float = _CM1A_SCALE,
) -> List[float]:
    """Apply the 1.14 scale factor and force the sum to ``total_charge`` exactly.

    CM1A already conserves charge, so after scaling the sum equals
    ``scale * total_charge``; for neutral molecules that is 0.  The residual
    correction only cleans up MOPAC's finite print precision.
    """
    scaled = [scale * q for q in cm1a_charges]
    if not scaled:
        return scaled
    residual = sum(scaled) - total_charge
    if abs(residual) > 1e-12:
        corr = residual / len(scaled)
        scaled = [s - corr for s in scaled]
    return scaled


# ---------------------------------------------------------------------------
# MOPAC output parsing (validated against OpenMOPAC 23.x output)
# ---------------------------------------------------------------------------

def parse_net_atomic_charges(out_text: str) -> List[float]:
    """Extract AM1 net atomic charges from a MOPAC ``.out`` file.

    Parses the ``NET ATOMIC CHARGES`` table (the diagonal-density / ZDO-Mulliken
    charge that matches AMSOL's ``QM``), not the deorthogonalised "MULLIKEN
    POPULATIONS AND CHARGES" section.
    """
    lines = out_text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if "NET ATOMIC CHARGES" in line:
            start = i
            break
    if start is None:
        raise QMChargeError("MOPAC output has no 'NET ATOMIC CHARGES' section")

    # Advance to the column header row, then read the table.
    idx = start
    while idx < len(lines) and "CHARGE" not in lines[idx].upper():
        idx += 1
    idx += 1

    charges: List[float] = []
    for line in lines[idx:]:
        parts = line.split()
        # Expected row: <num> <element> <charge> <electrons> ...
        if len(parts) >= 3 and parts[0].isdigit():
            try:
                charges.append(float(parts[2]))
                continue
            except ValueError:
                break
        if charges:  # table ended
            break
    if not charges:
        raise QMChargeError("Could not parse any net atomic charges from MOPAC output")
    return charges


def parse_bond_orders(out_text: str) -> Dict[Tuple[int, int], float]:
    """Extract Mayer bond orders from a MOPAC ``(VALENCIES) BOND ORDERS`` block.

    Returns a dict keyed by ``(i, j)`` with ``i < j`` (0-based atom indices).
    Handles MOPAC's per-atom blocks and their wrapped continuation lines.
    """
    lines = out_text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if "BOND ORDERS" in line:
            start = i + 1
            break
    if start is None:
        raise QMChargeError("MOPAC output has no 'BOND ORDERS' section")

    bonds: Dict[Tuple[int, int], float] = {}
    current: Optional[int] = None
    blanks = 0
    for line in lines[start:]:
        if not line.strip():
            blanks += 1
            # Two consecutive blank lines (or a blank after content) closes the
            # section; a single blank just separates atom blocks.
            if blanks >= 2:
                break
            continue
        blanks = 0
        upper = line.upper()
        if "MULLIKEN" in upper or "POPULATION" in upper or "DIPOLE" in upper:
            break

        tokens = line.split()
        # Atom-block header: "<i> <El> (<valency>) <j> <El> <order> ..."
        header_match = (
            len(tokens) >= 3
            and tokens[0].isdigit()
            and tokens[2].startswith("(")
        )
        if header_match:
            current = int(tokens[0]) - 1
            rest = tokens[3:]
        elif tokens and tokens[0].isdigit() and current is not None:
            # Continuation line: bond triples for the current atom.
            rest = tokens
        else:
            continue

        # Consume "<j> <El> <order>" triples from `rest`.
        k = 0
        while k + 3 <= len(rest):
            try:
                j = int(rest[k]) - 1
                order = float(rest[k + 2])
            except ValueError:
                break
            if current is not None and j != current:
                key = (min(current, j), max(current, j))
                bonds[key] = order
            k += 3
    if not bonds:
        raise QMChargeError("Could not parse any bond orders from MOPAC output")
    return bonds


# ---------------------------------------------------------------------------
# MOPAC driver
# ---------------------------------------------------------------------------

class QMChargeAssigner:
    """Assign 1.14*CM1A partial charges to a molecule using MOPAC's AM1.

    Args:
        mopac_cmd: MOPAC executable name or path.  Defaults to ``$MOPAC_CMD`` or
            ``"mopac"`` discovered on ``PATH``.
        optimize: Run AM1 geometry optimisation before computing charges
            (strict LigParGen behaviour).  Default ``False`` — single-point AM1
            on the supplied geometry (see module docstring).
        keep_files: Keep the temporary MOPAC input/output files (for debugging).
        scratch_dir: Directory for temporary MOPAC files (default: a fresh
            system temp directory, removed afterwards unless ``keep_files``).
        timeout: Maximum seconds to wait for MOPAC.
    """

    def __init__(
        self,
        mopac_cmd: Optional[str] = None,
        optimize: bool = False,
        keep_files: bool = False,
        scratch_dir: Optional[str] = None,
        timeout: float = 600.0,
    ):
        self.mopac_cmd = mopac_cmd or os.environ.get("MOPAC_CMD") or "mopac"
        self.optimize = optimize
        self.keep_files = keep_files
        self.scratch_dir = scratch_dir
        self.timeout = timeout

    # -- public API ---------------------------------------------------------

    def assign(
        self,
        mol: Chem.Mol,
        coords: np.ndarray,
        atom_types: Optional[Dict[int, str]] = None,
    ) -> Dict[int, float]:
        """Compute 1.14*CM1A charges for *mol* at geometry *coords* (Ångströms).

        Args:
            mol: RDKit molecule (provides elements and connectivity).
            coords: ``(N, 3)`` array of Cartesian coordinates in Ångströms,
                ordered to match ``mol`` atom indices.
            atom_types: Unused; accepted for signature parity with
                :meth:`biochar.ml_charges.MLChargeRefinement.refine` and
                :meth:`biochar.opls_typing.ChargeAssigner.assign_charges`.

        Returns:
            ``{atom_idx: charge}`` for every atom in *mol*.

        Raises:
            QMChargeError: If MOPAC is not found or the run/parse fails.
        """
        exe = shutil.which(self.mopac_cmd)
        if exe is None:
            raise QMChargeError(
                f"MOPAC executable {self.mopac_cmd!r} not found on PATH. "
                "Install it with `conda install -c conda-forge mopac` or from "
                "https://github.com/openmopac/mopac, or set $MOPAC_CMD."
            )

        n = mol.GetNumAtoms()
        coords = np.asarray(coords, dtype=float)
        if coords.shape != (n, 3):
            raise ValueError(
                f"coords shape {coords.shape} does not match {n} atoms"
            )

        atomic_numbers = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in range(n)]
        total_charge = Chem.GetFormalCharge(mol)

        out_text = self._run_mopac(exe, mol, coords, total_charge)

        qm = parse_net_atomic_charges(out_text)
        if len(qm) != n:
            raise QMChargeError(
                f"MOPAC returned {len(qm)} charges for {n} atoms"
            )
        bond_orders = parse_bond_orders(out_text)

        cm1a = cm1a_from_am1(qm, bond_orders, atomic_numbers)
        scaled = scale_and_neutralize(cm1a, total_charge=float(total_charge))

        logger.info(
            "1.14*CM1A charges assigned via MOPAC AM1 (%s): total charge %.4f e",
            "optimised" if self.optimize else "single-point",
            sum(scaled),
        )
        return {i: scaled[i] for i in range(n)}

    # -- internals ----------------------------------------------------------

    def _build_input(
        self, mol: Chem.Mol, coords: np.ndarray, total_charge: int
    ) -> str:
        """Compose a MOPAC AM1 input deck."""
        keywords = ["AM1", "MULLIK", "BONDS", "GEO-OK", f"CHARGE={total_charge}"]
        if not self.optimize:
            keywords.insert(1, "1SCF")
        opt_flag = 1 if self.optimize else 0

        header = " ".join(keywords)
        lines = [header, "biochar 1.14*CM1A charge job", ""]
        for i in range(mol.GetNumAtoms()):
            sym = mol.GetAtomWithIdx(i).GetSymbol()
            x, y, z = coords[i]
            lines.append(
                f"{sym:<2} {x:14.8f} {opt_flag} {y:14.8f} {opt_flag} "
                f"{z:14.8f} {opt_flag}"
            )
        lines.append("")
        return "\n".join(lines)

    def _run_mopac(
        self,
        exe: str,
        mol: Chem.Mol,
        coords: np.ndarray,
        total_charge: int,
    ) -> str:
        """Write input, invoke MOPAC, and return the ``.out`` text."""
        deck = self._build_input(mol, coords, total_charge)
        tmpdir = tempfile.mkdtemp(prefix="biochar_qm_", dir=self.scratch_dir)
        mop_path = Path(tmpdir) / "charge.mop"
        out_path = Path(tmpdir) / "charge.out"
        mop_path.write_text(deck)

        try:
            proc = subprocess.run(
                [exe, str(mop_path)],
                cwd=tmpdir,
                capture_output=True,
                text=True,
                timeout=self.timeout,
            )
            if not out_path.exists():
                raise QMChargeError(
                    f"MOPAC produced no output file.\n"
                    f"stdout: {proc.stdout}\nstderr: {proc.stderr}"
                )
            out_text = out_path.read_text(errors="replace")
            if "JOB ENDED NORMALLY" not in out_text and "== MOPAC DONE ==" not in out_text:
                logger.warning(
                    "MOPAC did not report a normal exit; parsing output anyway."
                )
            return out_text
        except subprocess.TimeoutExpired as exc:
            raise QMChargeError(
                f"MOPAC timed out after {self.timeout}s"
            ) from exc
        finally:
            if not self.keep_files:
                shutil.rmtree(tmpdir, ignore_errors=True)
            else:
                logger.info("Kept MOPAC scratch files in %s", tmpdir)
