Functional Groups
=================

Two modes
---------

**Ratio-driven (default)**: Set ``O_C_ratio`` and leave
``functional_groups=None``.  Oxygen is placed as phenolic (–OH) groups
until the ratio is met.

**Explicit dict**: Pass ``functional_groups`` as a mapping of group name to
exact placement count.  The ``O_C_ratio`` parameter is ignored.

.. code-block:: python

   from biochar_simulator.biochar_generator import generate_biochar

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

+------------+-------+-------------------------------------------------------------+
| Name       | O     | Description                                                 |
+============+=======+=============================================================+
| phenolic   | 1     | Aromatic C–OH.  Most common; always succeeds on edge sites. |
+------------+-------+-------------------------------------------------------------+
| hydroxyl   | 1     | Alias for phenolic on pure aromatic PAH.                    |
+------------+-------+-------------------------------------------------------------+
| carboxyl   | 2     | Ar–C(=O)(OH); adds one extra carbon and two oxygens.        |
+------------+-------+-------------------------------------------------------------+
| ether      | 1     | Ar–O–Ar bridge linking two edge carbons (5-membered ring).  |
+------------+-------+-------------------------------------------------------------+
| carbonyl   | 1     | Falls back to phenolic (pure aromatic edge has no C=O site).|
+------------+-------+-------------------------------------------------------------+
| quinone    | 2     | Falls back to phenolic.                                     |
+------------+-------+-------------------------------------------------------------+
| lactone    | 2     | Falls back to phenolic.                                     |
+------------+-------+-------------------------------------------------------------+

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
