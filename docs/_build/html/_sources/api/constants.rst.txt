constants
=========

OPLS-AA force field parameters, PAH library entries, and physical constants
used throughout the package.

.. rubric:: OPLS-AA parameters

``OPLS_ATOM_TYPES``
    ``dict[str, tuple]`` — atom type → ``(description, mass_amu, default_charge)``

``OPLS_LJ_PARAMS``
    ``dict[str, tuple]`` — atom type → ``(sigma_nm, epsilon_kJ/mol)``

``OPLS_BOND_PARAMS``
    ``dict[tuple, tuple]`` — ``(type1, type2)`` → ``(k_kcal/mol/Å², r0_Å)``

``OPLS_ANGLE_PARAMS``
    ``dict[tuple, tuple]`` — ``(type1, type2, type3)`` → ``(k_kcal/mol/rad², θ₀_deg)``

.. rubric:: Functional groups

``FUNCTIONAL_GROUPS``
    ``dict[str, dict]`` — group name → definition dict with keys
    ``atoms``, ``connectivity``, ``composition``, ``O_per_group``,
    ``H_per_group``.

.. rubric:: PAH library

``PAH_LIBRARY``
    ``dict[str, dict]`` — molecule name → ``{smiles, num_atoms,
    num_aromatic, molecular_formula, references}``.
    Contains 18 validated entries from benzene (6 C) to hex_lattice_40 (40 C).

.. rubric:: Physical constants

``CARBON_VDW_DIAMETER``
    ``float`` — 3.4 Å — graphite interlayer spacing, used as effective
    sheet thickness when computing slit-pore geometry.

.. rubric:: Helper functions

.. autofunction:: biochar_simulator.constants.get_atom_mass

.. autofunction:: biochar_simulator.constants.get_vdw_radius

.. autofunction:: biochar_simulator.constants.get_covalent_radius
