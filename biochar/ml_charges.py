"""
ML-based partial charge refinement for biochar structures.

Provides optional environment-aware charge assignment as an alternative to
static OPLS-AA table lookups.  An atomic featuriser maps each atom to a
feature vector; a pre-trained Gaussian Process Regression (GPR) model then
predicts the partial charge for each atom.

The bundled model (``biochar/data/charges_gpr_cm5.pkl``) is trained on
OPLS-AA reference charges for representative PAH/biochar fragments and serves
as a functional baseline.  For higher accuracy, retrain on DFT-derived (CM5,
RESP, or HLY) charges using :meth:`MLChargeRefinement.train_and_save`.

Requires the ``ml`` optional extra::

    pip install "biochar[ml]"
"""

from __future__ import annotations

import logging
import pickle
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Optional

import numpy as np
from rdkit import Chem

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)

_DEFAULT_MODEL_PATH = Path(__file__).parent / "data" / "charges_gpr_cm5.pkl"


class MLChargeRefinement:
    """
    Refine per-atom partial charges using a trained Gaussian Process model.

    Args:
        model_path: Path to a serialised scikit-learn ``Pipeline`` (.pkl).
                    Defaults to the bundled model at
                    ``biochar/data/charges_gpr_cm5.pkl``.

    Examples::

        refiner = MLChargeRefinement()
        charges = refiner.refine(mol, atom_types)
    """

    def __init__(self, model_path: Optional[Path] = None):
        path = Path(model_path) if model_path is not None else _DEFAULT_MODEL_PATH
        self._model = self._load_model(path)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def refine(self, mol: Chem.Mol, atom_types: Dict[int, str]) -> Dict[int, float]:
        """
        Predict per-atom partial charges for *mol*.

        The total charge is constrained to zero by an additive correction
        distributed equally across all atoms after prediction.

        Args:
            mol: RDKit molecule.
            atom_types: Mapping ``{atom_idx: opls_type}`` from
                :class:`~biochar.opls_typing.AtomTyper`.

        Returns:
            Dict ``{atom_idx: charge}``, same structure as
            :meth:`~biochar.opls_typing.ChargeAssigner.assign_charges`.
        """
        X = self._featurize(mol, atom_types)
        raw_q = self._model.predict(X).astype(float)

        # Enforce charge neutrality via uniform correction
        correction = raw_q.sum() / len(raw_q)
        raw_q -= correction

        return {idx: float(raw_q[idx]) for idx in range(mol.GetNumAtoms())}

    # ------------------------------------------------------------------
    # Feature engineering
    # ------------------------------------------------------------------

    def _featurize(self, mol: Chem.Mol, atom_types: Dict[int, str]) -> np.ndarray:
        """
        Build feature matrix of shape ``(n_atoms, 8)``.

        Columns:
            0  atomic_num           int
            1  is_aromatic          0/1
            2  in_ring              0/1
            3  smallest_ring_size   0 if acyclic, else 5/6/…
            4  num_heavy_neighbors  int
            5  num_h_neighbors      int
            6  formal_charge        int
            7  opls_type_group      0=C, 1=H, 2=O, 3=N, 4=S, 5=other
        """
        ring_info = mol.GetRingInfo()
        rows = []

        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            atomic_num = atom.GetAtomicNum()
            is_aromatic = int(atom.GetIsAromatic())
            in_ring = int(ring_info.NumAtomRings(idx) > 0)

            if in_ring:
                sizes = [len(r) for r in ring_info.AtomRings() if idx in r]
                smallest_ring = min(sizes)
            else:
                smallest_ring = 0

            neighbors = list(atom.GetNeighbors())
            n_heavy = sum(1 for n in neighbors if n.GetAtomicNum() != 1)
            n_h = sum(1 for n in neighbors if n.GetAtomicNum() == 1)
            formal_charge = atom.GetFormalCharge()

            opls_type = atom_types.get(idx, "")
            group = self._opls_group(opls_type)

            rows.append([
                atomic_num, is_aromatic, in_ring, smallest_ring,
                n_heavy, n_h, formal_charge, group,
            ])

        return np.array(rows, dtype=float)

    @staticmethod
    def _opls_group(opls_type: str) -> int:
        """Coarse atom-group integer from OPLS type string."""
        if opls_type in ("CA", "CT", "C"):
            return 0
        if opls_type in ("HA", "HO", "HO2", "HC", "HNA", "HSH", "HNPR"):
            return 1
        if opls_type in ("OH", "OS", "OC", "O", "OH2", "OW"):
            return 2
        if opls_type in ("NA", "N", "NT", "NPY", "NPR", "NGR"):
            return 3
        if opls_type in ("SH_", "SS"):
            return 4
        return 5

    # ------------------------------------------------------------------
    # Model I/O
    # ------------------------------------------------------------------

    @staticmethod
    def _load_model(path: Path):
        """Load a serialised sklearn pipeline from *path*."""
        if path.exists():
            with open(path, "rb") as fh:
                model = pickle.load(fh)
            logger.debug("Loaded ML charge model from %s", path)
            return model

        logger.warning(
            "Bundled ML charge model not found at %s; "
            "building fallback GPR from OPLS reference data.",
            path,
        )
        return MLChargeRefinement._build_fallback_model()

    @staticmethod
    def _build_fallback_model():
        """
        Train and return a GPR pipeline on OPLS-AA reference charges.

        Used when the bundled .pkl is absent.  Reproduces OPLS-AA charges
        from atomic features — a functional baseline that can be retrained
        with DFT data via :meth:`train_and_save`.
        """
        _require_sklearn()
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import RBF, WhiteKernel
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler

        X_train, y_train = _generate_training_data()

        kernel = RBF(length_scale=1.0) + WhiteKernel(noise_level=1e-3)
        gpr = GaussianProcessRegressor(
            kernel=kernel, n_restarts_optimizer=2, normalize_y=True
        )
        pipe = Pipeline([("scaler", StandardScaler()), ("gpr", gpr)])
        pipe.fit(X_train, y_train)
        return pipe

    @classmethod
    def train_and_save(
        cls,
        X: np.ndarray,
        y: np.ndarray,
        output_path: Optional[Path] = None,
    ) -> "MLChargeRefinement":
        """
        Train a GPR on features *X* and charges *y*, then save to disk.

        Call this with DFT-derived charges to replace the bundled baseline::

            refiner = MLChargeRefinement.train_and_save(X_dft, y_dft)

        Args:
            X: Feature matrix ``(n_samples, 8)`` from :meth:`_featurize`.
            y: Per-atom target charges ``(n_samples,)``.
            output_path: Where to write the .pkl.  Defaults to
                ``biochar/data/charges_gpr_cm5.pkl``.

        Returns:
            Fitted :class:`MLChargeRefinement` instance.
        """
        _require_sklearn()
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import RBF, WhiteKernel
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler

        kernel = RBF(length_scale=1.0) + WhiteKernel(noise_level=1e-3)
        gpr = GaussianProcessRegressor(
            kernel=kernel, n_restarts_optimizer=3, normalize_y=True
        )
        pipe = Pipeline([("scaler", StandardScaler()), ("gpr", gpr)])
        pipe.fit(X, y)

        path = output_path or _DEFAULT_MODEL_PATH
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "wb") as fh:
            pickle.dump(pipe, fh)
        logger.info("ML charge model saved to %s", path)

        instance = cls.__new__(cls)
        instance._model = pipe
        return instance


# ---------------------------------------------------------------------------
# Training data generation
# ---------------------------------------------------------------------------

def _generate_training_data() -> tuple[np.ndarray, np.ndarray]:
    """
    Build (X, y) from representative PAH/biochar structures using OPLS charges
    as targets.  Returns arrays suitable for ``sklearn`` fit calls.
    """
    from .constants import OPLS_ATOM_TYPES
    from .opls_typing import AtomTyper, ChargeAssigner

    smiles_list = [
        "c1ccccc1",               # benzene
        "c1ccc2ccccc2c1",         # naphthalene
        "c1ccc2cc3ccccc3cc2c1",   # anthracene
        "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61",  # coronene
        "Oc1ccccc1",              # phenol
        "c1ccc(N)cc1",            # aniline
        "OC(=O)c1ccccc1",         # benzoic acid
        "c1ccc(Oc2ccccc2)cc1",    # diphenyl ether
        "c1ccc(S)cc1",            # thiophenol
        "c1ccc(Sc2ccccc2)cc1",    # diphenyl sulfide
        "Nc1cc2ccccc2cc1",        # 2-aminonaphthalene
        "c1ccncc1",               # pyridine (pyridinic N model)
        "c1cc[nH]c1",             # pyrrole (pyrrolic N model)
    ]

    typer = AtomTyper()
    charger = ChargeAssigner()
    X_rows: list = []
    y_rows: list = []

    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        mol = Chem.AddHs(mol)
        atom_types = typer.assign_atom_types(mol)
        charges = charger.assign_charges(mol, atom_types)

        ring_info = mol.GetRingInfo()

        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            atomic_num = atom.GetAtomicNum()
            is_aromatic = int(atom.GetIsAromatic())
            in_ring = int(ring_info.NumAtomRings(idx) > 0)

            if in_ring:
                sizes = [len(r) for r in ring_info.AtomRings() if idx in r]
                smallest_ring = min(sizes)
            else:
                smallest_ring = 0

            neighbors = list(atom.GetNeighbors())
            n_heavy = sum(1 for n in neighbors if n.GetAtomicNum() != 1)
            n_h = sum(1 for n in neighbors if n.GetAtomicNum() == 1)
            formal_charge = atom.GetFormalCharge()
            opls_type = atom_types.get(idx, "")
            group = MLChargeRefinement._opls_group(opls_type)

            X_rows.append([
                atomic_num, is_aromatic, in_ring, smallest_ring,
                n_heavy, n_h, formal_charge, group,
            ])
            y_rows.append(charges[idx])

    return np.array(X_rows, dtype=float), np.array(y_rows, dtype=float)


def _require_sklearn() -> None:
    """Raise a clear ImportError if scikit-learn is not installed."""
    try:
        import sklearn  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "scikit-learn is required for ML-based charge refinement. "
            "Install it with: pip install 'biochar[ml]'"
        ) from exc


def build_and_save_bundled_model(output_path: Optional[Path] = None) -> Path:
    """
    Build the bundled GPR model and save it to *output_path*.

    Intended to be run once during development to regenerate
    ``biochar/data/charges_gpr_cm5.pkl``::

        python -c "from biochar.ml_charges import build_and_save_bundled_model; build_and_save_bundled_model()"

    Returns:
        Path where the model was saved.
    """
    _require_sklearn()
    path = output_path or _DEFAULT_MODEL_PATH
    MLChargeRefinement.train_and_save(
        *_generate_training_data(), output_path=path
    )
    print(f"Bundled model saved to {path}")
    return path
