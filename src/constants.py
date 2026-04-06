"""
OPLS-AA Force Field Constants

Atom types, charges, and masses for OPLS-AA force field.
Reference: Jorgensen, W. L., et al. JACS 118.45 (1996): 11225-11236.
"""

import numpy as np

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

    # Nitrogen (for future extensions)
    "N": ("Tertiary nitrogen", 14.007, -0.70),
    "NT": ("Quaternary nitrogen", 14.007, 0.0),

    # Sulfur (for future extensions)
    "S": ("Sulfur sp3", 32.065, -0.20),
    "SH": ("H on sulfur", 1.008, 0.16),
}

# OPLS-AA Lennard-Jones Parameters
# Format: atom_type -> (sigma_nm, epsilon_kJ/mol)
# sigma: van der Waals radius (in nanometers)
# epsilon: well depth (in kJ/mol)
# Reference: GROMACS OPLS-AA forcefield

OPLS_LJ_PARAMS = {
    "CA": (0.3550, 0.2928),      # Aromatic carbon
    "HA": (0.2600, 0.0630),      # Aromatic hydrogen
    "CT": (0.3500, 0.2761),      # Aliphatic carbon (sp3)
    "HC": (0.2500, 0.0630),      # H on aliphatic C
    "OC": (0.2960, 0.5023),      # Carbonyl oxygen
    "OH": (0.3066, 0.7113),      # Hydroxyl oxygen
    "OS": (0.2960, 0.5023),      # Ether oxygen
    "OW": (0.3066, 0.6276),      # Water oxygen
    "HO": (0.0000, 0.0000),      # H in hydroxyl (no LJ)
    "C": (0.3750, 0.4392),       # Carboxylic acid carbonyl C
    "O": (0.2960, 0.5023),       # Carboxylic acid carbonyl O
    "OH2": (0.3066, 0.7113),     # Hydroxyl on carboxylic acid
    "HO2": (0.0000, 0.0000),     # H on carboxylic hydroxyl (no LJ)
    "N": (0.3250, 0.7113),       # Tertiary nitrogen
    "NT": (0.3250, 0.7113),      # Quaternary nitrogen
    "S": (0.3550, 1.0460),       # Sulfur sp3
    "SH": (0.0000, 0.0000),      # H on sulfur (no LJ)
}

# OPLS-AA Bond Parameters (k_bond, r0)
# Format: (atom_type1, atom_type2) -> (k_kcal/mol/Angstrom^2, r0_Angstrom)
# Values from GROMACS OPLS-AA forcefield

OPLS_BOND_PARAMS = {
    ("CA", "CA"): (770.0, 1.387),      # Aromatic C-C
    ("CA", "HA"): (367.0, 1.084),      # Aromatic C-H
    ("CT", "CT"): (268.0, 1.529),      # Aliphatic C-C
    ("CT", "HC"): (367.0, 1.084),      # Aliphatic C-H
    ("CA", "CT"): (268.0, 1.529),      # Aromatic-aliphatic C-C
    ("CT", "OH"): (320.0, 1.410),      # Aliphatic C-OH
    ("CA", "OH"): (450.0, 1.367),      # Aromatic C-OH (phenolic)
    ("CT", "OS"): (320.0, 1.410),      # Aliphatic C-O (ether)
    ("CA", "OS"): (450.0, 1.367),      # Aromatic C-O (ether)
    ("C", "O"): (750.0, 1.230),        # Carboxylic carbonyl
    ("C", "OH2"): (320.0, 1.364),      # C-OH in carboxylic acid
    ("CT", "OC"): (320.0, 1.410),      # Aliphatic C-carbonyl O
    ("OH", "HO"): (367.0, 0.960),      # O-H hydroxyl
    ("OH2", "HO2"): (367.0, 0.960),    # O-H carboxylic
    ("OS", "HO"): (367.0, 0.960),      # This shouldn't exist, but for safety
}

# OPLS-AA Angle Parameters (k_angle, theta0)
# Format: (atom_type1, atom_type2, atom_type3) -> (k_kcal/mol/rad^2, theta0_deg)

OPLS_ANGLE_PARAMS = {
    ("CA", "CA", "CA"): (126.0, 120.0),
    ("CA", "CA", "HA"): (70.0, 120.0),
    ("CA", "CA", "CT"): (70.0, 120.0),
    ("CA", "CA", "OH"): (70.0, 119.7),
    ("CA", "CA", "OS"): (70.0, 119.7),
    ("HA", "CA", "HA"): (35.0, 120.0),
    ("CT", "CT", "CT"): (63.0, 109.47),
    ("CT", "CT", "HC"): (48.0, 109.47),
    ("CT", "CT", "OH"): (55.0, 109.47),
    ("CT", "CT", "OS"): (55.0, 109.47),
    ("HC", "CT", "HC"): (35.0, 109.47),
    ("HC", "CT", "OH"): (35.0, 109.47),
    ("CT", "OH", "HO"): (55.0, 104.52),
    ("CA", "OH", "HO"): (70.0, 108.0),
    ("CT", "OS", "CT"): (60.0, 110.7),
    ("CA", "OS", "CA"): (70.0, 113.2),
    ("CT", "OS", "CA"): (70.0, 113.2),
}

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
