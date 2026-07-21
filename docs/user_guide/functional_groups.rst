Functional Groups
=================

Two modes
---------

**Ratio-driven (default)**: Set ``O_C_ratio`` and leave
``functional_groups=None``.  Oxygen fills the aromatic edge sites as phenolic
(–OH) groups first; if the target still is not met — typical of low-temperature,
low-aromaticity chars, whose edges are largely consumed by aliphatic decoration
and hydrogen saturation — the remainder is placed as aliphatic hydroxyls
(–CH\ :sub:`2`\ OH) on the sp3 carbons.  This spill is on by default and can be
disabled with ``allow_aliphatic_oxygen=False`` (see
:ref:`aliphatic-oxygen`).

**Explicit dict**: Pass ``functional_groups`` as a mapping of group name to
exact placement count.  The ``O_C_ratio`` parameter is ignored.

.. code-block:: python

   from biochar import generate_biochar

   # Ratio-driven: 10 % O/C, all placed as phenolic
   mol, coords, gro, top, itp = generate_biochar(
       target_num_carbons=60,
       O_C_ratio=0.10,
   )

   # Explicit: exactly 3 phenolic + 1 carboxyl + 2 ether
   mol, coords, gro, top, itp = generate_biochar(
       target_num_carbons=60,
       functional_groups={"phenolic": 3, "carboxyl": 1, "ether": 2},
   )

Available groups
----------------

+--------------------+-------+-------------------------------------------------------+
| Name               | O     | Description                                           |
+====================+=======+=======================================================+
| phenolic           | 1     | Aromatic C–OH.  Most common; always succeeds on edge  |
|                    |       | sites.                                                |
+--------------------+-------+-------------------------------------------------------+
| hydroxyl           | 1     | Alias for phenolic on pure aromatic PAH.              |
+--------------------+-------+-------------------------------------------------------+
| aliphatic_hydroxyl | 1     | Primary alcohol –CH\ :sub:`2`\ OH on an sp3 carbon (a |
|                    |       | pendant methyl becomes hydroxymethyl).  The dominant  |
|                    |       | O group in low-temperature, cellulose-derived         |
|                    |       | biochar.  Requires aliphatic carbons in the skeleton  |
|                    |       | (see :ref:`aliphatic-oxygen`).                        |
+--------------------+-------+-------------------------------------------------------+
| carboxyl           | 2     | Ar–C(=O)(OH); adds one extra carbon and two oxygens.  |
+--------------------+-------+-------------------------------------------------------+
| ether              | 1     | Ar–O–Ar bridge linking two edge carbons (5-membered   |
|                    |       | ring).                                                |
+--------------------+-------+-------------------------------------------------------+
| carbonyl           | 1     | Falls back to phenolic (pure aromatic edge has no C=O |
|                    |       | site).                                                |
+--------------------+-------+-------------------------------------------------------+
| quinone            | 2     | Falls back to phenolic.                               |
+--------------------+-------+-------------------------------------------------------+
| lactone            | 2     | Falls back to phenolic.                               |
+--------------------+-------+-------------------------------------------------------+

Ether bridge geometry
---------------------

The ether group bridges two ring-edge carbons through one oxygen atom,
forming a cyclic C–O–C linkage.  By default the maximum allowed
carbon-skeleton distance between the two bridged carbons is 3 bonds
(``max_ether_span=3``), producing a furan-like 5-membered ring that stays
geometrically flat.

.. code-block:: python

   # Default (furan-like, always flat)
   mol, coords, gro, top, itp = generate_biochar(
       target_num_carbons=60,
       functional_groups={"ether": 2},
       max_ether_span=3,   # default
   )

   # Pyran-like 6-membered bridge (may introduce minor strain)
   mol, coords, gro, top, itp = generate_biochar(
       target_num_carbons=60,
       functional_groups={"ether": 2},
       max_ether_span=4,
   )

.. warning::

   ``max_ether_span`` values > 4 can produce long-range bridges that fold
   the aromatic sheet into a nanotube-like shape.  The default of 3 is
   recommended for flat PAH biochar.

.. _aliphatic-oxygen:

Oxygen on aliphatic carbons
---------------------------

The groups above except ``aliphatic_hydroxyl`` attach only to aromatic edge
carbons.  On low-temperature, low-aromaticity chars (roughly aromaticity below
70 %) the hydrogen-shaping stage decorates the skeleton with pendant sp3
carbons and saturates the remaining edges with hydrogen, so very few aromatic
edge sites are left for oxygen.  A high ``O_C_ratio`` then could not be reached
on aromatic edges alone — e.g. hardwood at 300 °C (target O/C ≈ 0.25) could
place only a single oxygen — and the point would be reported as failing its
composition tolerance.

To reach these targets, ratio-driven placement spills the shortfall onto the
sp3 carbons as ``aliphatic_hydroxyl`` groups: each pendant –CH\ :sub:`3`
becomes a hydroxymethyl –CH\ :sub:`2`\ OH.  This is the primary alcohol that
dominates the oxygen content of cellulose-derived biochar, and it uses only the
existing OPLS types (the carbon stays ``CT``; the O–H is the same ``OH``/``HO``
as a phenol), so no force-field change is involved.

The spill is on by default.  Set ``allow_aliphatic_oxygen=False`` to keep
oxygen on aromatic edges only, reproducing the earlier behaviour; it has no
effect on purely aromatic skeletons, which have no sp3 carbons.

.. code-block:: python

   # Reach a high O/C on a low-aromaticity char (default: spill enabled)
   mol, coords, gro, top, itp = generate_biochar(
       temperature=300, feedstock="hardwood",   # target O/C ≈ 0.25
   )

   # Keep oxygen on aromatic edges only
   mol, coords, gro, top, itp = generate_biochar(
       temperature=300, feedstock="hardwood",
       allow_aliphatic_oxygen=False,
   )

In a sweep this is a normal ``GeneratorConfig`` field, so it can go under
``fixed`` (or an axis) in the sweep config just like any other parameter.

Mixing group types
------------------

Any combination of supported groups can be used together:

.. code-block:: python

   mol, coords, gro, top, itp = generate_biochar(
       target_num_carbons=100,
       functional_groups={
           "phenolic": 4,
           "carboxyl": 2,
           "ether":    3,
       },
       seed=42,
   )
