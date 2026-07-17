"""
OPLS-AA Force Field Constants

Atom types, charges, and masses for OPLS-AA force field.
Reference: Jorgensen, W. L., et al. JACS 118.45 (1996): 11225-11236.
"""


# OPLS-AA Atom Types for Biochar Systems
# Format: atom_type: (description, mass_amu, default_charge)

OPLS_ATOM_TYPES = {
    # Aromatic carbons and hydrogens
    "CA": ("Aromatic carbon", 12.01, -0.115),
    "HA": ("Aromatic hydrogen", 1.008, 0.115),

    # Aliphatic carbons
    "CT": ("Aliphatic C (sp3)", 12.01, -0.18),
    "HC": ("H on aliphatic C", 1.008, 0.06),

    # Oxygens
    "OC": ("Carbonyl oxygen", 15.999, -0.56),
    "OH": ("Hydroxyl oxygen", 15.999, -0.661),
    "OS": ("Ether oxygen", 15.999, -0.322),
    "OW": ("Water oxygen", 15.999, -0.820),

    # Hydroxyls and carboxylic acid
    "HO": ("H in hydroxyl", 1.008, 0.436),
    "C": ("Carboxylic acid carbonyl C", 12.01, 0.620),
    "O": ("Carboxylic acid carbonyl O", 15.999, -0.540),
    "OH2": ("Hydroxyl on carboxylic acid", 15.999, -0.540),
    "HO2": ("H on carboxylic hydroxyl", 1.008, 0.436),

    # Nitrogen
    "N": ("Tertiary nitrogen", 14.007, -0.70),
    "NT": ("Quaternary nitrogen", 14.007, 0.0),
    "NA": ("Aromatic amine nitrogen (aniline Ar-NH2)", 14.007, -0.60),
    "HNA": ("H on aromatic amine nitrogen", 1.008, 0.30),

    # Sulfur
    "SH_": ("Aromatic thiol sulfur (Ar-SH)", 32.06, -0.39),    # opls_734
    "HSH": ("H on thiol sulfur", 1.008, 0.21),                 # opls_204
    "SS": ("Thioether sulfur bridging two aryl C (Ar-S-Ar)", 32.06, -0.16),  # opls_222
    # Ring-substituting nitrogen (biochar N-doping)
    "NPY": ("Pyridinic N (substituted into 6-ring, no H)", 14.007, -0.36),
    "NPR": ("Pyrrolic N (substituted into 5-ring, with H)", 14.007, -0.52),
    "NGR": ("Graphitic/quaternary N (interior 6-ring, no H)", 14.007, 0.02),
    "HNPR": ("H on pyrrolic nitrogen", 1.008, 0.38),

}

# Lennard-Jones, bond, and angle parameters intentionally live in oplsaa.ff, not
# here. GROMACS resolves them from the #included forcefield by the opls_XXX name in
# GROMACS_OPLS_TYPE_MAP below, so a hand-copied table here would be dead weight that
# can only drift. A previous table did drift: its values were a mix of AMBER and OPLS
# with no single provenance. Add parameters to a local .itp instead of reviving it.

# Mapping from internal generic type names to GROMACS OPLS-AA opls_XXX names.
# Internal names (CA, HA, etc.) are used throughout the biochar generation pipeline.
# At GROMACS export time, these are translated so the topology is compatible with
# the standard oplsaa.ff forcefield shipped with GROMACS.
#
# When the exported .itp #includes a real oplsaa.ff, only the opls_XXX name reaches
# GROMACS -- mass, charge, LJ and bonded parameters are all resolved from the
# forcefield by that name. A name naming the wrong element therefore silently
# simulates the wrong chemistry. Trailing comments quote the type's description in
# oplsaa.ff/atomtypes.atp verbatim; keep them in sync when editing.
# tests/test_opls_type_map.py checks every entry against an installed oplsaa.ff.
GROMACS_OPLS_TYPE_MAP: dict[str, str] = {
    "CA":  "opls_145",   # aromatic carbon, benzene-type
    "HA":  "opls_146",   # aromatic hydrogen
    "CT":  "opls_135",   # aliphatic sp3 carbon
    "HC":  "opls_140",   # hydrogen on aliphatic carbon
    "OH":  "opls_154",   # phenol / alcohol oxygen
    "HO":  "opls_155",   # phenol / alcohol O-H hydrogen
    "OS":  "opls_467",   # aryl ether oxygen (Ar-O-Ar or Ar-O-R)
    "OC":  "opls_278",   # ketone / aldehyde C=O oxygen
    "C":   "opls_267",   # carboxylic acid carbonyl carbon
    "O":   "opls_269",   # carboxylic acid C=O oxygen
    "OH2": "opls_268",   # carboxylic acid -OH oxygen
    "HO2": "opls_270",   # carboxylic acid -OH hydrogen
    "OW":  "opls_116",   # SPC/E water oxygen
    "NA":  "opls_900",   # "N primary   amines" -- aniline Ar-NH2 nitrogen
    "HNA": "opls_909",   # "H(N)   primary   amines" -- H on that nitrogen
    "SH_": "opls_734",   # "all-atom S: thiophenol (HS is #204)" -- Ar-SH sulfur
    "HSH": "opls_204",   # "all-atom H(S): thiols" -- the HS named by opls_734
    "SS":  "opls_222",   # "S in thioanisoles" -- nearest aryl-S; see note below
    "NPY": "opls_520",   # "N   in pyridine 6-31G*" -- pyridinic ring N
    "NPR": "opls_542",   # "N   in pyrrole" -- pyrrolic ring N
    "HNPR": "opls_545",  # "H1  in pyrrole" -- the pyrrole N-H hydrogen
    "NGR": "opls_520",   # graphitic/quaternary ring N; approximated, see note below
}

# Two entries above are deliberate approximations, not exact matches:
#
#   SS  -> opls_222 is the thioanisole sulfur (Ar-S-CH3). OPLS-AA has no diaryl
#          thioether (Ar-S-Ar) type; opls_222 is the only aryl-attached sulfide S
#          and carries the matching CA-S bond (0.176 nm, "thioanisole" in
#          ffbonded.itp).
#
#          KNOWN GAP: the bond resolves but the angle does not. A thioether
#          emits a CA-S-CA angle, and ffbonded.itp has no such angletype --
#          only CA-S-CT and CA-S-CM (both 104.200 deg, 518.816 kJ/mol/rad^2).
#          grompp therefore rejects a thioether topology with "No default
#          Angle types". Supply the angle from a local .itp; the closest
#          stock value is the CA-S-CT thioanisole angle above.
#
#   NGR -> opls_520 is the pyridine N, reused for graphitic/quaternary N because
#          OPLS-AA has no substitutional 3-coordinate aromatic N type. Element and
#          ring aromaticity are right; the charge/bonded environment is only
#          approximate. NPY maps here too -- for pyridinic N it is exact.

# Functional groups definitions
# Each functional group specifies how to add atoms to the carbon skeleton
FUNCTIONAL_GROUPS = {
    "hydroxyl": {
        "description": "Hydroxyl group (-OH)",
        "atoms": [("O", "OH"), ("H", "HO")],  # (atom_type, group_code)
        "connectivity": [(0, "C", 1), (0, 1, 1)],  # (atom_idx, connect_to, bond_type)
        "composition": {"O": 1, "H": 1},
        "O_per_group": 1,
        "H_per_group": 1,
    },
    "carboxyl": {
        "description": "Carboxylic acid group (-COOH)",
        "atoms": [("C", "C"), ("O", "O"), ("O", "OH2"), ("H", "HO2")],
        "connectivity": [
            (0, "C", 2),  # C=O double bond
            (1, "C", 1),  # Single bond C-C
            (2, "C", 1),  # O-H bond
            (3, 2, 1),
        ],
        "composition": {"C": 1, "O": 2, "H": 1},
        "O_per_group": 2,
        "H_per_group": 1,
    },
    "phenolic": {
        "description": "Phenolic group (aromatic -OH)",
        "atoms": [("O", "OH"), ("H", "HO")],
        "connectivity": [(0, "CA", 1), (0, 1, 1)],
        "composition": {"O": 1, "H": 1},
        "O_per_group": 1,
        "H_per_group": 1,
    },
    "ether": {
        "description": "Ether group (-O-)",
        "atoms": [("O", "OS")],
        "connectivity": [(0, "C", 1), (0, "C", 1)],  # Two C-O bonds
        "composition": {"O": 1},
        "O_per_group": 1,
        "H_per_group": 0,
    },
    "carbonyl": {
        "description": "Carbonyl group (C=O)",
        "atoms": [("O", "OC")],
        "connectivity": [(0, "C", 2)],  # C=O double bond
        "composition": {"O": 1},
        "O_per_group": 1,
        "H_per_group": 0,
    },
    "lactone": {
        "description": "Lactone group (cyclic ester)",
        "atoms": [("C", "C"), ("O", "O"), ("O", "OS")],
        "connectivity": [
            (0, "C", 2),
            (1, "C", 1),
            (2, "C", 1),
        ],
        "composition": {"C": 1, "O": 2},
        "O_per_group": 2,
        "H_per_group": 0,
    },
    "quinone": {
        "description": "Quinone group (C=O with aromatic C)",
        "atoms": [("O", "OC"), ("O", "OC")],
        "connectivity": [(0, "CA", 2), (1, "CA", 2)],
        "composition": {"O": 2},
        "O_per_group": 2,
        "H_per_group": 0,
    },
    "amino": {
        "description": "Amino group (-NH2) on aromatic ring",
        "atoms": [("N", "NA"), ("H", "HNA"), ("H", "HNA")],
        "connectivity": [(0, "CA", 1), (1, 0, 1), (2, 0, 1)],
        "composition": {"N": 1, "H": 2},
        "O_per_group": 0,
        "H_per_group": 2,
    },
    "thiol": {
        "description": "Thiol group (-SH) on aromatic ring",
        "atoms": [("S", "SH_"), ("H", "HSH")],
        "connectivity": [(0, "CA", 1), (1, 0, 1)],
        "composition": {"S": 1, "H": 1},
        "O_per_group": 0,
        "H_per_group": 1,
    },
    "thioether": {
        "description": "Thioether group (-S-) bridging two aromatic carbons",
        "atoms": [("S", "SS")],
        "connectivity": [(0, "C", 1), (0, "C", 1)],  # Two C-S bonds
        "composition": {"S": 1},
        "O_per_group": 0,
        "H_per_group": 0,
    },
}

# Common PAH structures (SMILES notation)
# All entries validated: correct carbon count, 100% aromatic, all atoms in 6-membered rings.
# Hex-lattice entries are programmatically generated compact graphene-nanoflake topologies.
PAH_LIBRARY = {
    # --- 6 carbons ---
    "benzene": {
        "smiles": "c1ccccc1",
        "num_atoms": 6,
        "num_aromatic": 6,
        "molecular_formula": "C6H6",
        "references": "Basic aromatic ring",
    },
    # --- 10 carbons ---
    "naphthalene": {
        "smiles": "c1ccc2ccccc2c1",
        "num_atoms": 10,
        "num_aromatic": 10,
        "molecular_formula": "C10H8",
        "references": "Two fused rings",
    },
    # --- 14 carbons ---
    "anthracene": {
        "smiles": "c1ccc2cc3ccccc3cc2c1",
        "num_atoms": 14,
        "num_aromatic": 14,
        "molecular_formula": "C14H10",
        "references": "Three fused rings (linear acene)",
    },
    "phenanthrene": {
        "smiles": "c1ccc2ccc3ccccc3c2c1",
        "num_atoms": 14,
        "num_aromatic": 14,
        "molecular_formula": "C14H10",
        "references": "Three fused rings (angular)",
    },
    # --- 16 carbons ---
    "pyrene": {
        "smiles": "c1cc2ccc3cccc4ccc(c1)c2c34",
        "num_atoms": 16,
        "num_aromatic": 16,
        "molecular_formula": "C16H10",
        "references": "Four fused rings (pericondensed)",
    },
    # --- 18 carbons ---
    "chrysene": {
        "smiles": "c1ccc2c(c1)cc1ccc3ccccc3c1c2",
        "num_atoms": 18,
        "num_aromatic": 18,
        "molecular_formula": "C18H12",
        "references": "Four fused rings (chrysene topology)",
    },
    "tetracene": {
        "smiles": "c1ccc2cc3cc4ccccc4cc3cc2c1",
        "num_atoms": 18,
        "num_aromatic": 18,
        "molecular_formula": "C18H12",
        "references": "Four fused rings (linear acene, naphthacene)",
    },
    "triphenylene": {
        "smiles": "c1ccc2c(c1)c1ccccc1c1ccccc21",
        "num_atoms": 18,
        "num_aromatic": 18,
        "molecular_formula": "C18H12",
        "references": "Four fused rings (triphenylene topology)",
    },
    # --- 22 carbons ---
    "pentacene": {
        "smiles": "c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1",
        "num_atoms": 22,
        "num_aromatic": 22,
        "molecular_formula": "C22H14",
        "references": "Five fused rings (linear acene)",
    },
    "picene": {
        "smiles": "c1ccc2cc3ccc4ccc5ccccc5c4c3cc2c1",
        "num_atoms": 22,
        "num_aromatic": 22,
        "molecular_formula": "C22H14",
        "references": "Five fused rings (picene/[5]helicene topology)",
    },
    "hex_lattice_22": {
        "smiles": "c1cc2ccc3ccc4ccc5cccc6c(c1)c2c3c4c56",
        "num_atoms": 22,
        "num_aromatic": 22,
        "molecular_formula": "C22H12",
        "references": "Six-ring compact graphene nanoflake (hex-lattice seed)",
    },
    # --- 24 carbons ---
    "coronene": {
        "smiles": "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61",
        "num_atoms": 24,
        "num_aromatic": 24,
        "molecular_formula": "C24H12",
        "references": "Seven fused rings (central + 6 surrounding)",
    },
    # --- 26 carbons ---
    "hexacene": {
        "smiles": "c1ccc2cc3cc4cc5cc6ccccc6cc5cc4cc3cc2c1",
        "num_atoms": 26,
        "num_aromatic": 26,
        "molecular_formula": "C26H16",
        "references": "Six fused rings (linear acene)",
    },
    "dibenzo_bc_ef_coronene": {
        "smiles": "c1ccc2cc3cc4ccc5ccc6ccccc6c5c4cc3cc2c1",
        "num_atoms": 26,
        "num_aromatic": 26,
        "molecular_formula": "C26H14",
        "references": "Six fused rings (angular pericondensed topology)",
    },
    # --- 28 carbons ---
    "hex_lattice_28": {
        "smiles": "c1ccc2c(c1)c1ccc3ccc4ccc5ccc6ccc2c2c6c5c4c3c12",
        "num_atoms": 28,
        "num_aromatic": 28,
        "molecular_formula": "C28H14",
        "references": "Eight-ring compact graphene nanoflake (hex-lattice seed)",
    },
    # --- 30 carbons ---
    "hex_lattice_30": {
        "smiles": "c1cc2ccc3cc4ccc5ccc6ccc7ccc8c(c1)c2c3c1c4c5c6c7c81",
        "num_atoms": 30,
        "num_aromatic": 30,
        "molecular_formula": "C30H14",
        "references": "Nine-ring compact graphene nanoflake (hex-lattice seed)",
    },
    # --- 38 carbons ---
    "hex_lattice_38": {
        "smiles": "c1cc2cc3ccc4cc5cccc6c7ccc8ccc9ccc%10c(c1)c2c1c3c4c(c56)c2c7c8c9c%10c12",
        "num_atoms": 38,
        "num_aromatic": 38,
        "molecular_formula": "C38H16",
        "references": "Twelve-ring compact graphene nanoflake (hex-lattice seed)",
    },
    # --- 40 carbons ---
    "hex_lattice_40": {
        "smiles": "c1cc2cc3ccc4cc5ccc6ccc7cc8ccc9ccc%10c(c1)c2c1c3c4c2c5c6c7c3c8c9c%10c1c32",
        "num_atoms": 40,
        "num_aromatic": 40,
        "molecular_formula": "C40H16",
        "references": "Thirteen-ring compact graphene nanoflake (hex-lattice seed)",
    },
}

# Atomic masses (in amu)
ATOMIC_MASSES = {
    "H": 1.008,
    "C": 12.01,
    "N": 14.007,
    "O": 15.999,
    "S": 32.065,
    "P": 30.974,
    "Cl": 35.45,
    "Br": 79.904,
}

# Van der Waals radii (in Angstrom) for steric clash detection
VDW_RADII = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "P": 1.80,
    "Cl": 1.75,
    "Br": 1.85,
}

# Van der Waals diameter of graphitic carbon (graphite interlayer spacing).
# Used as effective sheet thickness when computing slit-pore geometry.
CARBON_VDW_DIAMETER = 3.4  # Angstroms

# ---------------------------------------------------------------------------
# Experimental-data model provenance & tunables
# ---------------------------------------------------------------------------
# Primary characterization dataset behind the temperature/feedstock composition
# model (see :mod:`biochar.temperature_model`):
UC_DAVIS_DB_URL = "https://biochar.ucdavis.edu/"   # UC Davis Biochar Database
# Methodological parent: Wood, Mašek & Erastova, Cell Reports Physical Science
# 5(7), 2024, DOI 10.1016/j.xcrp.2024.102036.
#
# Minimum aromaticity (%) the PAH skeleton builder can faithfully realise.  When
# the data-derived aromaticity (from temperature/feedstock) falls below this it
# is clamped to this floor and a warning is emitted.
MIN_BUILDABLE_AROMATICITY = 70.0

# Ring-curvature ratios of non-graphitizing carbons (Wood et al. 2024): the
# published island building blocks average roughly 10 hexagons : 2 pentagons :
# 1 heptagon.  As per-ring-addition probabilities over the 13-ring total this is
# 2/13 pentagons and 1/13 heptagons (hexagons are the ~10/13 remainder).  Pass
# these as ``defect_fraction`` / ``heptagon_fraction`` to reproduce that mix.
WOOD_PENTAGON_FRACTION = 2.0 / 13.0   # ≈ 0.154
WOOD_HEPTAGON_FRACTION = 1.0 / 13.0   # ≈ 0.077

# Covalent radii (in Angstrom) for bond length validation
COVALENT_RADII = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "S": 1.05,
    "P": 1.07,
    "Cl": 1.02,
    "Br": 1.20,
}


def get_atom_mass(atom_type: str) -> float:
    """Get atomic mass from OPLS atom type."""
    if atom_type in OPLS_ATOM_TYPES:
        return OPLS_ATOM_TYPES[atom_type][1]
    # Extract element from atom type (first character or two)
    element = atom_type[0] if atom_type[0].isupper() else atom_type[:2]
    return ATOMIC_MASSES.get(element, 12.01)


def get_vdw_radius(element: str) -> float:
    """Get Van der Waals radius for element."""
    return VDW_RADII.get(element, 1.70)


def get_covalent_radius(element: str) -> float:
    """Get covalent radius for element."""
    return COVALENT_RADII.get(element, 0.76)
