# GROMACS Simulation Workflow with Biochar Generator

## Generated Files

The Biochar Generator produces three files:

1. **structure.gro** - 3D atomic coordinates in GROMACS format
2. **structure.top** - Topology file with molecule definition and forcefield references
3. **structure.itp** - Molecule definition file (included in .top)

## Atom Type Resolution

The generator **intelligently handles atom type definitions**:

- If your forcefield file (`oplsaa.ff/forcefield.itp`) already defines `[ atomtypes ]`, they are used automatically
- The generated `.top` file references the forcefield via `#include "oplsaa.ff/forcefield.itp"`
- GROMACS merges atom type definitions from the entire include chain
- **No duplicate atom type definitions** - cleaner and avoids conflicts

## Standard GROMACS Workflow

Assuming you have the OPLS-AA forcefield installed:

```bash
# 1. Energy minimization
gmx grompp -f em.mdp -c structure.gro -p structure.top -o em.tpr
gmx mdrun -deffnm em

# 2. NVT equilibration (constant volume)
gmx grompp -f nvt.mdp -c em.gro -p structure.top -o nvt.tpr
gmx mdrun -deffnm nvt

# 3. NPT equilibration (constant pressure)
gmx grompp -f npt.mdp -c nvt.gro -p structure.top -o npt.tpr
gmx mdrun -deffnm npt

# 4. Production MD
gmx grompp -f md.mdp -c npt.gro -p structure.top -o md.tpr
gmx mdrun -deffnm md
```

## Adding Other Molecules

You can add other molecules to your system using standard GROMACS tools:

### Option A: Add Solvent (Water)

```bash
gmx solvate -cp structure.gro -cs spc216.gro -p structure.top -o structure_solv.gro
```

This works because:
- Water atom types (OW, HW) are standard OPLS-AA types
- They're defined in your forcefield's `[ atomtypes ]` section
- Solvent `.top` is automatically merged with your biochar `.top`

### Option B: Add Other Molecules

If adding molecules with different atom types:

1. **Ensure atom types are defined** - They must be in the include chain:
   - Your `oplsaa.ff/forcefield.itp` file, OR
   - The other molecule's forcefield file

2. **Check atom type compatibility** - Use `gmx pdb2gmx` or manually edit if needed:
   ```bash
   gmx pdb2gmx -f other_molecule.pdb -p other.top -ff oplsaa
   ```

3. **Merge topologies** - Combine multiple `.top` files:
   ```
   [ system ]
   Biochar + Other Molecules

   [ molecules ]
   BC 1
   OTHER 2
   SOL 1000
   ```

## Troubleshooting

### "Atom type 'XX' not found"
- Check that the atom type is defined in the forcefield `[ atomtypes ]` section
- Verify the forcefield file path is correct
- Run: `grep "\[ atomtypes \]" oplsaa.ff/forcefield.itp`

### "Unresolved 'atom_type' in atom definition"
- The atom type is referenced but not defined anywhere
- Manually add it to your `.top` file's `[ atomtypes ]` section

## Forcefield Include Chain

The include chain is processed in this order:

```
your_structure.top
  └─ #include "oplsaa.ff/forcefield.itp"
      └─ #include "oplsaa.ff/atomtypes.atp"
      └─ #include "oplsaa.ff/ffbonded.itp"
      └─ #include "oplsaa.ff/ffnonbonded.itp"
```

All `[ atomtypes ]` sections in this chain are merged by GROMACS.

## Expected Atom Types for Biochar

Your generated molecules use these OPLS-AA atom types:

| Type | Description | Element |
|------|-------------|---------|
| CA | Aromatic carbon | C |
| HA | Aromatic hydrogen | H |
| CT | Aliphatic carbon (sp3) | C |
| HC | Hydrogen on aliphatic C | H |
| OH | Hydroxyl oxygen | O |
| HO | Hydrogen in hydroxyl | H |
| OS | Ether oxygen | O |
| OC | Carbonyl oxygen | O |
| C | Carboxylic acid carbonyl | C |
| O | Carboxylic acid oxygen | O |

All these must be available in `oplsaa.ff/atomtypes.atp` or defined in your `.top` file.
