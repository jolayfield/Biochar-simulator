import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path

ANA = Path('/Users/layf0001/Claude Cowork/Biochar-simulator/sim_200/analysis')
FIG = ANA / 'figures'
FIG.mkdir(exist_ok=True)

BLUE   = '#2166AC'
RED    = '#D6604D'
GREEN  = '#4DAC26'
ORANGE = '#F4A582'
GRAY   = '#636363'

def read_xvg(path):
    t, y = [], []
    with open(path) as f:
        for line in f:
            if line.startswith(('#', '@')): continue
            vals = line.split()
            if len(vals) >= 2:
                t.append(float(vals[0]))
                y.append(float(vals[1]))
    return np.array(t), np.array(y)

# ── 1. Energy panel (2×2) ────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle('Wet NPT Equilibration — Thermodynamic Properties\n(200-molecule biochar bilayer in SPC/E water)', fontsize=14, fontweight='bold')

panels = [
    ('energy_temperature.xvg', 'Temperature', 'Temperature (K)',     BLUE,   axes[0,0]),
    ('energy_pressure.xvg',    'Pressure',    'Pressure (bar)',       RED,    axes[0,1]),
    ('energy_density.xvg',     'Density',     'Density (kg m⁻³)',     GREEN,  axes[1,0]),
    ('energy_potential.xvg',   'Potential\nEnergy', 'Potential Energy (kJ mol⁻¹)', ORANGE, axes[1,1]),
]

for fname, title, ylabel, color, ax in panels:
    t, y = read_xvg(ANA / fname)
    ax.plot(t, y, color=color, lw=1.2, alpha=0.85, label='Simulation')
    mean = y.mean()
    ax.axhline(mean, color='k', lw=1.0, ls='--', label=f'Mean: {mean:.2f}')
    ax.set_xlabel('Time (ps)', fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=9)

plt.tight_layout()
plt.savefig(FIG / 'energy_panels.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved energy_panels.png')

# ── 2. SASA over time ────────────────────────────────────────────────────────
t_sasa, sasa = read_xvg(ANA / 'sasa.xvg')

fig, ax = plt.subplots(figsize=(9, 4.5))
ax.plot(t_sasa, sasa, color=BLUE, lw=1.5, label='SASA')
ax.axhline(sasa.mean(), color=RED, lw=1.2, ls='--',
           label=f'Mean: {sasa.mean():.1f} nm²')
ax.fill_between(t_sasa, sasa.mean()-sasa.std(), sasa.mean()+sasa.std(),
                color=BLUE, alpha=0.15, label=f'±1σ ({sasa.std():.1f} nm²)')
ax.set_xlabel('Time (ps)', fontsize=12)
ax.set_ylabel('Solvent-Accessible Surface Area (nm²)', fontsize=12)
ax.set_title('Biochar Solvent-Accessible Surface Area\n(200 molecules, probe radius 0.14 nm)', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
per_mol = sasa.mean() / 200
ax.text(0.02, 0.05, f'Per molecule: {per_mol:.1f} nm²', transform=ax.transAxes,
        fontsize=10, color=GRAY, bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
plt.tight_layout()
plt.savefig(FIG / 'sasa.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved sasa.png')

# ── 3. RDF ───────────────────────────────────────────────────────────────────
r, g = read_xvg(ANA / 'rdf_C_OW.xvg')

fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(r, g, color=BLUE, lw=2.0, label='g(r) — Biochar C → Water O$_W$')
ax.axhline(1.0, color=GRAY, lw=1.0, ls='--', alpha=0.7, label='Bulk (g = 1)')
ax.fill_between(r, g, 1.0, where=(g < 1.0), alpha=0.12, color=RED,
                label='Depletion zone (g < 1)')
ax.set_xlabel('Distance r (nm)', fontsize=12)
ax.set_ylabel('Radial Distribution Function g(r)', fontsize=12)
ax.set_title('Biochar Carbon → Water Oxygen RDF\n(hydrophobic character of aromatic surface)', fontsize=13, fontweight='bold')
ax.set_xlim(0, 2.0)
ax.set_ylim(0, 1.2)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.text(0.55, 0.15,
        'g(r) < 1 throughout indicates\nwater depletion near hydrophobic\naromatic biochar surface',
        transform=ax.transAxes, fontsize=9, color=GRAY,
        bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow', alpha=0.85))
plt.tight_layout()
plt.savefig(FIG / 'rdf.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved rdf.png')

# ── 4. Density profiles — water and biochar separately, then combined ────────
z_wat, rho_wat = read_xvg(ANA / 'density_water_z.xvg')
z_bc,  rho_bc  = read_xvg(ANA / 'density_biochar_z.xvg')

# Normalise biochar to kg/m³ (gmx density outputs kg/m³ already)
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.suptitle('Mass Density Profiles along Z-axis\n(200-molecule biochar bilayer in SPC/E water)',
             fontsize=13, fontweight='bold')

# Panel A: water only
ax = axes[0]
ax.plot(z_wat, rho_wat, color=BLUE, lw=2.0)
ax.axhline(rho_wat[rho_wat > 50].mean(), color=BLUE, lw=1.0, ls='--', alpha=0.6,
           label=f'Bulk mean: {rho_wat[rho_wat>50].mean():.0f} kg m⁻³')
ax.set_xlabel('Z position (nm)', fontsize=12)
ax.set_ylabel('Mass Density (kg m⁻³)', fontsize=12)
ax.set_title('Water (SPC/E)', fontsize=12, fontweight='bold', color=BLUE)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(z_wat.min(), z_wat.max())

# Panel B: biochar only
ax = axes[1]
ax.plot(z_bc, rho_bc, color=RED, lw=2.0)
ax.set_xlabel('Z position (nm)', fontsize=12)
ax.set_ylabel('Mass Density (kg m⁻³)', fontsize=12)
ax.set_title('Biochar (200 molecules)', fontsize=12, fontweight='bold', color=RED)
ax.grid(True, alpha=0.3)
ax.set_xlim(z_bc.min(), z_bc.max())
# Annotate slab extent
bc_nonzero = z_bc[rho_bc > rho_bc.max()*0.05]
if len(bc_nonzero) > 0:
    ax.axvspan(bc_nonzero.min(), bc_nonzero.max(), alpha=0.08, color=RED,
               label=f'Slab: {bc_nonzero.min():.1f}–{bc_nonzero.max():.1f} nm')
    ax.legend(fontsize=9)

# Panel C: overlay (dual y-axis)
ax = axes[2]
ax2 = ax.twinx()
l1, = ax.plot(z_wat, rho_wat, color=BLUE, lw=2.0, label='Water')
l2, = ax2.plot(z_bc,  rho_bc,  color=RED,  lw=2.0, label='Biochar')
ax.set_xlabel('Z position (nm)', fontsize=12)
ax.set_ylabel('Water Density (kg m⁻³)', fontsize=12, color=BLUE)
ax2.set_ylabel('Biochar Density (kg m⁻³)', fontsize=12, color=RED)
ax.tick_params(axis='y', colors=BLUE)
ax2.tick_params(axis='y', colors=RED)
ax.set_title('Water + Biochar Overlay', fontsize=12, fontweight='bold')
ax.legend(handles=[l1, l2], loc='upper left', fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(min(z_wat.min(), z_bc.min()), max(z_wat.max(), z_bc.max()))

plt.tight_layout()
plt.savefig(FIG / 'density_profiles.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved density_profiles.png')

print('\nAll figures saved to:', FIG)
