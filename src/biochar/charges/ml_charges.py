"""
ML-based partial charge refinement for biochar structures.

Provides optional environment-aware charge assignment as an alternative to
static OPLS-AA table lookups.  An atomic featuriser maps each atom to a
feature vector; a pre-trained Gaussian Process Regression (GPR) model then
predicts the partial charge for each atom.

The bundled model (``biochar/data/charges_gpr_cm5.json``) is trained on
OPLS-AA reference charges for representative PAH/biochar fragments and serves
as a functional baseline.  For higher accuracy, retrain on DFT-derived (CM5,
RESP, or HLY) charges using :meth:`MLChargeRefinement.train_and_save`.

Requires the ``ml`` optional extra::

    pip install "biochar[ml]"
"""

from __future__ import annotations

import json
import logging
import pickle
import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Optional

import numpy as np
import sklearn
from rdkit import Chem

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)

_DEFAULT_MODEL_PATH = Path(__file__).parent.parent / "data" / "charges_gpr_cm5.json"

#: The artifact schema this module reads and writes.  A pickle is only as
#: portable as the library that wrote it; this format is the fitted
#: hyperparameters and the training set, from which the pipeline is rebuilt on
#: load.  The rebuild is not an approximation of the pickle -- refitting with the
#: hyperparameters held fixed reproduces its predictions to ~2e-15 e.
_ARTIFACT_FORMAT = "biochar-gpr-1"

#: How far the reconstruction may drift from the reference charges recorded in
#: the artifact before it stops being the same model.  Well below anything a
#: partial charge means, well above the ~1e-15 e of Cholesky round-off.
_REPRODUCTION_TOLERANCE = 1e-6


class MLChargeRefinement:
    """
    Refine per-atom partial charges using a trained Gaussian Process model.

    Args:
        model_path: Path to a saved model.  Defaults to the bundled one at
                    ``biochar/data/charges_gpr_cm5.json``.  Pickles written by
                    older versions of :meth:`train_and_save` are still read.

    Examples::

        refiner = MLChargeRefinement()
        charges = refiner.refine(mol, atom_types)
    """

    def __init__(self, model_path: Optional[Path] = None):
        path = Path(model_path) if model_path is not None else _DEFAULT_MODEL_PATH
        #: Which model answered: ``"bundled"``, ``"custom"`` or ``"fallback"``.
        #: A caller who asked for one and got another has nothing else to go on
        #: -- the charges look the same either way.
        self.model_source: str = (
            "bundled" if path == _DEFAULT_MODEL_PATH else "custom"
        )
        if not path.exists():
            self.model_source = "fallback"
        self.model_path: Path = path
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
                :class:`~biochar.pipeline.opls_typing.AtomTyper`.

        Returns:
            Dict ``{atom_idx: charge}``, same structure as
            :meth:`~biochar.pipeline.opls_typing.ChargeAssigner.assign_charges`.
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
    def _recorded_sklearn_version(source) -> Optional[str]:
        """The scikit-learn behind the artifact, if it was not this one.

        *source* is whatever reading the artifact yielded, and there are two
        kinds. For the current format it is the parsed payload, which records
        the version that fitted the hyperparameters outright. For a legacy
        pickle it is the list of warnings raised while unpickling: scikit-learn
        stamps its version onto every estimator it pickles and raises
        ``InconsistentVersionWarning`` under a different one, and the stamp is
        popped from the object during unpickling, so the warning is the only
        place the recorded version can still be read.
        """
        if isinstance(source, dict):
            recorded = source.get("sklearn_version")
            return str(recorded) if recorded is not None else None

        for warning in source:
            recorded = getattr(warning.message, "original_sklearn_version", None)
            if recorded is not None:
                return str(recorded)
        return None

    @staticmethod
    def _load_model(path: Path):
        """Load a charge model from *path*.

        Two formats are read. The current one is JSON: the fitted kernel
        hyperparameters, the training set, and the charges the fitting library
        predicted from them. The pipeline is rebuilt from those with the
        hyperparameters held fixed, which is a deterministic Cholesky solve
        rather than a deserialisation, and the rebuilt predictions are checked
        against the recorded ones -- so what the artifact means is verified
        rather than assumed. The other is a pickle, written by
        ``train_and_save`` before this module stopped producing them.

        The two are reported differently, and the difference is the point. A
        pickle is only as portable as the library that wrote it, it records
        nothing to check itself against, and the constraint on scikit-learn is a
        floor with no ceiling -- so for a pickle a version difference is all
        there is to go on, and the strong warning stands. For the current format
        the version is provenance and the reproduction check is the evidence,
        so silence means the charges were verified rather than assumed.
        """
        if path.exists():
            model, source, deviation = MLChargeRefinement._read_model(path)

            if deviation is not None and deviation > _REPRODUCTION_TOLERANCE:
                warnings.warn(
                    f"The ML charge model at {path} does not reproduce under "
                    f"scikit-learn {sklearn.__version__}: rebuilding it from "
                    f"its recorded hyperparameters gives charges differing "
                    f"from the ones recorded alongside them by up to "
                    f"{deviation:.2e} e, against a tolerance of "
                    f"{_REPRODUCTION_TOLERANCE:.0e} e. The charges it predicts "
                    f"may be invalid. Rebuild it with "
                    f"biochar.charges.ml_charges.build_and_save_bundled_model(), "
                    f"or use charge_method='opls' or 'qm'.",
                    UserWarning,
                    stacklevel=3,
                )

            recorded = MLChargeRefinement._recorded_sklearn_version(source)
            if recorded is not None and recorded != sklearn.__version__:
                warnings.warn(
                    MLChargeRefinement._version_mismatch_message(
                        path, recorded, deviation
                    ),
                    UserWarning,
                    stacklevel=3,
                )
            if not isinstance(source, dict):
                for warning in source:
                    if getattr(
                        warning.message, "original_sklearn_version", None
                    ) is None:
                        warnings.warn_explicit(
                            warning.message, warning.category,
                            warning.filename, warning.lineno,
                        )

            logger.debug("Loaded ML charge model from %s", path)
            return model

        warnings.warn(
            f"No ML charge model at {path}; falling back to a Gaussian process "
            f"fitted at run time from OPLS reference data. It is the same recipe "
            f"as the bundled model but not the same model: it follows the "
            f"current reference table, where the bundled one is pinned to the "
            f"table as it stood when it was built. Check "
            f"MLChargeRefinement.model_source to see which answered.",
            UserWarning,
            stacklevel=3,
        )
        logger.warning(
            "ML charge model not found at %s; building fallback GPR from "
            "OPLS reference data.",
            path,
        )
        return MLChargeRefinement._build_fallback_model()

    @staticmethod
    def _version_mismatch_message(
        path: Path, recorded: str, deviation: Optional[float]
    ) -> str:
        """What a version difference means for the charges, per format.

        For a pickle it means the unpickling itself is unverified, which is the
        strong claim. For the current format the numbers were only *fitted* by
        the other library -- the model in hand was rebuilt by this one, and the
        rebuild was checked -- so the message says what the check found rather
        than repeating a warning the evidence does not support.
        """
        if deviation is None:
            return (
                f"The ML charge model at {path} was written by scikit-learn "
                f"{recorded} and is being loaded under {sklearn.__version__}. "
                f"It is a pickle, so the charges it predicts may be invalid. "
                f"Rebuild it with "
                f"biochar.charges.ml_charges.build_and_save_bundled_model() to "
                f"convert it to the version-independent format, or use "
                f"charge_method='opls' or 'qm'."
            )
        return (
            f"The ML charge model at {path} records hyperparameters fitted by "
            f"scikit-learn {recorded}, and scikit-learn {sklearn.__version__} "
            f"rebuilt the model from them. The rebuilt model reproduces the "
            f"charges recorded alongside those hyperparameters to "
            f"{deviation:.1e} e, so the predictions are the artifact's; this "
            f"is a note on provenance, not a defect."
        )

    @staticmethod
    def _read_model(path: Path):
        """Read *path*, returning ``(model, source, deviation)``.

        *source* is what :meth:`_recorded_sklearn_version` reads the recording
        library from. *deviation* is how far the rebuilt model's predictions sit
        from the ones recorded in the artifact, or ``None`` for a pickle, which
        records nothing to check against.

        The format is decided by content rather than by suffix: ``train_and_save``
        writes wherever the caller points it, and callers point it at ``.pkl``.
        """
        raw = path.read_bytes()
        if raw.lstrip()[:1] == b"{":
            payload = json.loads(raw)
            model = MLChargeRefinement._rebuild_from_payload(payload)
            reference = np.asarray(payload["reference_charges"], dtype=float)
            deviation = float(
                np.abs(model.predict(np.asarray(payload["X"], dtype=float))
                       - reference).max()
            )
            return model, payload, deviation

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            model = pickle.loads(raw)
        return model, caught, None

    @staticmethod
    def _rebuild_from_payload(payload: dict):
        """Rebuild the GPR pipeline from recorded hyperparameters and data.

        ``optimizer=None`` holds the kernel at the recorded hyperparameters, so
        the fit is the one deterministic linear-algebra step -- no restarts, no
        random state, no dependence on the optimiser a given scikit-learn ships.
        """
        _require_sklearn()
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import RBF, WhiteKernel
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler

        recorded_format = payload.get("format")
        if recorded_format != _ARTIFACT_FORMAT:
            raise ValueError(
                f"Unsupported ML charge model format {recorded_format!r}; this "
                f"biochar reads {_ARTIFACT_FORMAT!r}."
            )

        kernel = (
            RBF(length_scale=payload["kernel"]["rbf_length_scale"])
            + WhiteKernel(noise_level=payload["kernel"]["white_noise_level"])
        )
        gpr = GaussianProcessRegressor(
            kernel=kernel, optimizer=None, normalize_y=payload["normalize_y"]
        )
        pipe = Pipeline([("scaler", StandardScaler()), ("gpr", gpr)])
        pipe.fit(
            np.asarray(payload["X"], dtype=float),
            np.asarray(payload["y"], dtype=float),
        )
        return pipe

    @staticmethod
    def _write_model(pipe, X: np.ndarray, y: np.ndarray, path: Path) -> None:
        """Write *pipe* to *path* in the version-independent format.

        The reference charges are the fitting library's own predictions, and the
        rebuild is exercised here rather than left for the reader to discover:
        an artifact that does not reproduce on the machine that wrote it would
        not reproduce anywhere.
        """
        fitted = pipe.named_steps["gpr"].kernel_.get_params()
        payload = {
            "format": _ARTIFACT_FORMAT,
            "sklearn_version": sklearn.__version__,
            "normalize_y": bool(pipe.named_steps["gpr"].normalize_y),
            "kernel": {
                "rbf_length_scale": float(fitted["k1__length_scale"]),
                "white_noise_level": float(fitted["k2__noise_level"]),
            },
            "X": np.asarray(X, dtype=float).tolist(),
            "y": np.asarray(y, dtype=float).tolist(),
            "reference_charges": [float(q) for q in pipe.predict(X)],
        }

        deviation = float(
            np.abs(
                MLChargeRefinement._rebuild_from_payload(payload).predict(X)
                - np.asarray(payload["reference_charges"], dtype=float)
            ).max()
        )
        if deviation > _REPRODUCTION_TOLERANCE:
            raise RuntimeError(
                f"The model just fitted does not survive being written and "
                f"rebuilt: predictions move by up to {deviation:.2e} e against "
                f"a tolerance of {_REPRODUCTION_TOLERANCE:.0e} e. Refusing to "
                f"write {path}."
            )

        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload) + "\n")

    @staticmethod
    def _build_fallback_model():
        """
        Train and return a GPR pipeline on OPLS-AA reference charges.

        Used when the bundled artifact is absent.  Reproduces OPLS-AA charges
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
            output_path: Where to write the model.  Defaults to
                ``biochar/data/charges_gpr_cm5.json``.

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

        path = Path(output_path) if output_path is not None else _DEFAULT_MODEL_PATH
        cls._write_model(pipe, X, y, path)
        logger.info("ML charge model saved to %s", path)

        instance = cls.__new__(cls)
        instance._model = pipe
        instance.model_path = path
        instance.model_source = "bundled" if path == _DEFAULT_MODEL_PATH else "custom"
        return instance


# ---------------------------------------------------------------------------
# Training data generation
# ---------------------------------------------------------------------------

def _generate_training_data() -> tuple[np.ndarray, np.ndarray]:
    """
    Build (X, y) from representative PAH/biochar structures using OPLS charges
    as targets.  Returns arrays suitable for ``sklearn`` fit calls.
    """
    from ..pipeline.opls_typing import AtomTyper, ChargeAssigner

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
    ``biochar/data/charges_gpr_cm5.json``::

        python -c "from biochar.charges.ml_charges import build_and_save_bundled_model; build_and_save_bundled_model()"

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
