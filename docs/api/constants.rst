constants
=========

OPLS-AA force field parameters, PAH library entries, and physical constants
used throughout the package.

.. rubric:: OPLS-AA parameters

``OPLS_ATOM_TYPES``
    ``dict[str, tuple]`` — atom type → ``(description, mass_amu, default_charge)``

``GROMACS_OPLS_TYPE_MAP``
    ``dict[str, str]`` — internal atom type → GROMACS ``opls_XXX`` name.

.. note::

   Lennard-Jones, bond, and angle parameters are **not** stored in this module.
   An exported topology ``#include``\ s a real ``oplsaa.ff``, and GROMACS resolves
   every one of those parameters from it by the ``opls_XXX`` name above. Earlier
   releases carried hand-copied ``OPLS_LJ_PARAMS``, ``OPLS_BOND_PARAMS`` and
   ``OPLS_ANGLE_PARAMS`` tables; they were unused by the exporter and had drifted
   from the force field, so they were removed. Supply any parameter the force
   field lacks via a local ``.itp`` rather than reinstating a table here.

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

.. autofunction:: biochar.constants.get_atom_mass

.. autofunction:: biochar.constants.get_vdw_radius

.. autofunction:: biochar.constants.get_covalent_radius
