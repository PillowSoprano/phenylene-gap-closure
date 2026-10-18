"""Draw Figs. 1-3 from fig_data/*.csv only (no calculations).
Figures are produced at their final printed width so that no scaling happens in LaTeX."""
import csv, os
from pathlib import Path
os.chdir(Path(__file__).resolve().parent)
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D

FIG_W = 6.5            # inches; equals \textwidth of the manuscript. Change here for a journal template.
plt.rcParams.update({
    'text.usetex': True, 'font.family': 'serif',
    'text.latex.preamble': r'\usepackage{amsmath}',
    'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'legend.fontsize': 7, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.minor.width': 0.4, 'xtick.major.size': 3, 'ytick.major.size': 3, 'xtick.direction': 'out',
    'ytick.direction': 'out', 'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'savefig.bbox': None, 'savefig.pad_inches': 0, 'lines.linewidth': 1.2,
})
INK, MUTED, GRID = '#1f2937', '#8c8c8c', '#e6e6e6'
RED, BLUE, ORANGE = '#e34948', '#2a78d6', '#eb6834'
NEUTRAL_BOND = '#9aa3ad'
DIV = LinearSegmentedColormap.from_list('div', [BLUE, '#f0efec', RED])

def read(name):
    with open(f'fig_data/{name}') as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k, v in r.items():
            try: r[k] = float(v)
            except ValueError: pass
    return rows

def panel_letter(ax, letter, x=-0.02, y=1.0):
    ax.text(x, y, rf'\textbf{{{letter}}}', transform=ax.transAxes, fontsize=9, va='bottom', ha='left')

def draw_molecule(ax, atoms, bonds, colour_of, lw_of, dot=6):
    P = {int(a['atom']): (a['x'], a['y']) for a in atoms}
    for b in bonds:
        (x1, y1), (x2, y2) = P[int(b['atom1'])], P[int(b['atom2'])]
        ax.plot([x1, x2], [y1, y2], color=colour_of(b), lw=lw_of(b), solid_capstyle='round', zorder=2)
    xy = np.array(list(P.values()))
    ax.scatter(xy[:, 0], xy[:, 1], s=dot, color=INK, zorder=3, lw=0)
    ax.set_aspect('equal'); ax.axis('off')
    return xy

# ======================= Fig. 1 =======================
atoms, bonds = read('fig1_atoms.csv'), read('fig1_bonds.csv')
curves = read('fig1_curves.csv'); summ = {r['quantity']: r['value'] for r in read('fig1_summary.csv')}
rho = summ['rho_biphenylene']
fig = plt.figure(figsize=(FIG_W, 2.25))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.75], left=0.0, right=0.985, bottom=0.18, top=0.86, wspace=0.28)
axa, axb, axc = (fig.add_subplot(gs[0, i]) for i in range(3))
draw_molecule(axa, atoms, bonds, lambda b: RED if b['single_bond_route'] == 1 else NEUTRAL_BOND,
              lambda b: 3.4 if b['single_bond_route'] == 1 else 1.8)
draw_molecule(axb, atoms, bonds, lambda b: RED if b['cooperative_sign'] > 0 else BLUE, lambda b: 2.6)
for ax, letter, title in ((axa, 'a', r'one bond: $\times 4$'),
                          (axb, 'b', rf'all 14 bonds: $\pm{100*rho:.1f}\%$')):
    panel_letter(ax, letter, x=0.0, y=1.02); ax.text(0.12, 1.02, title, transform=ax.transAxes, va='bottom')
for ax in (axa, axb):
    ax.set_xlim(-1.35, 2.75); ax.set_ylim(-3.75, 1.35)
leg = [Line2D([], [], color=RED, lw=2.4, label='stronger'), Line2D([], [], color=BLUE, lw=2.4, label='weaker')]
axb.legend(handles=leg, frameon=False, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=2, handlelength=1.4, columnspacing=1.2, borderaxespad=0)
single = [(r['largest_relative_change_percent'], r['huckel_gap_abs_beta']) for r in curves if r['route'] == 'single']
coop = [(r['largest_relative_change_percent'], r['huckel_gap_abs_beta']) for r in curves if r['route'] == 'cooperative']
single, coop = np.array(single), np.array(coop)
axc.plot(*single.T, color=MUTED, lw=1.4, ls=(0, (4, 2)))
axc.plot(*coop.T, color=INK, lw=1.4)
axc.plot([single[-1, 0]], [single[-1, 1]], 'o', color=MUTED, ms=3.5)
axc.plot([coop[-1, 0]], [coop[-1, 1]], 'o', color=INK, ms=3.5)
axc.text(60, 0.47, r'one bond (a)', color=MUTED, ha='right')
axc.text(9.5, 0.12, r'all bonds (b)', color=INK, ha='right')
axc.set_xscale('log'); axc.set_xlim(0.5, 420); axc.set_ylim(-0.02, 0.95)
axc.set_xticks([1, 3, 10, 30, 100, 300]); axc.set_xticklabels(['1', '3', '10', '30', '100', '300'])
axc.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
axc.set_xlabel(r'largest relative change of a coupling (\%)'); axc.set_ylabel(r'H\"uckel gap ($|\beta|$)')
axc.yaxis.grid(True, color=GRID, lw=0.5); axc.set_axisbelow(True)
panel_letter(axc, 'c', x=-0.16, y=1.02)
fig.savefig('fig1_hook.pdf'); fig.savefig('fig1_hook.png', dpi=200)
plt.close(fig)

# ======================= Fig. 2 =======================
keys = [('linear3', 'linear [3]'), ('angular3', 'angular [3]'), ('linear4', 'linear [4]'),
        ('bent4', 'bent [4]'), ('zigzag4', 'zigzag [4]'), ('triangular4', 'triangular [4]')]
rhos = {r['skeleton']: r['rho'] for r in read('fig2_summary.csv')}
data = {k: (read(f'fig2_{k}_atoms.csv'), read(f'fig2_{k}_bonds.csv')) for k, _ in keys}
spans = []
for k, _ in keys:
    xy = np.array([(a['x'], a['y']) for a in data[k][0]]); spans.append(xy.max(0) - xy.min(0))
half = np.max(spans, axis=0)/2 + 0.45                     # one common scale for all panels
norm = TwoSlopeNorm(vcenter=0, vmin=-0.21, vmax=0.21)
WSCALE = 110                                              # marker area (pt^2) per unit closing-orbital weight
fig = plt.figure(figsize=(FIG_W, 5.0))
gs = fig.add_gridspec(2, 3, left=0.0, right=1.0, bottom=0.17, top=0.95, wspace=0.0, hspace=0.10)
for n, (k, label) in enumerate(keys):
    ax = fig.add_subplot(gs[n // 3, n % 3])
    atoms, bonds = data[k]
    xy = draw_molecule(ax, atoms, bonds, lambda b: DIV(norm(b['coupling_change'])), lambda b: 2.2, dot=4)
    w = np.array([a['closing_orbital_weight'] for a in atoms])
    ax.scatter(xy[:, 0], xy[:, 1], s=WSCALE*w*10, facecolor='none', edgecolor=INK, lw=0.6, zorder=4)
    c = (xy.max(0) + xy.min(0))/2
    ax.set_xlim(c[0]-half[0], c[0]+half[0]); ax.set_ylim(c[1]-half[1], c[1]+half[1])
    panel_letter(ax, 'abcdef'[n], x=0.04, y=1.0)
    ax.text(0.11, 1.0, rf'{label}, $\rho={rhos[k]:.3f}$', transform=ax.transAxes, va='bottom')
cax = fig.add_axes([0.06, 0.085, 0.42, 0.022])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=DIV), cax=cax, orientation='horizontal')
cb.set_label('relative change of the coupling at the closure threshold'); cb.outline.set_linewidth(0.5)
cb.set_ticks([-0.2, -0.1, 0, 0.1, 0.2])
lax = fig.add_axes([0.58, 0.05, 0.40, 0.06]); lax.axis('off')
for i, wv in enumerate((0.05, 0.10, 0.15)):
    lax.scatter([0.08 + 0.3*i], [0.5], s=WSCALE*wv*10, facecolor='none', edgecolor=INK, lw=0.6, clip_on=False)
    lax.text(0.08 + 0.3*i + 0.06, 0.5, rf'${wv:.2f}$', va='center')
lax.set_xlim(0, 1); lax.set_ylim(0, 1)
lax.text(0.0, 1.25, 'closing-orbital weight on an atom', va='bottom')
fig.savefig('fig2_closing.pdf'); fig.savefig('fig2_closing.png', dpi=200)
plt.close(fig)

# ======================= Fig. 3 =======================
lad, rnd = read('fig3_ladders.csv'), read('fig3_random.csv')
def ladder(mol, kind):
    pts = sorted([r for r in lad if r['molecule'] == mol and r['pattern'] == kind], key=lambda r: r['model_energy_eV'])
    return pts
fig = plt.figure(figsize=(FIG_W, 2.35))
gs = fig.add_gridspec(1, 3, left=0.075, right=0.985, bottom=0.30, top=0.90, wspace=0.38)
axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
names = {'biphenylene': 'biphenylene', 'angular3': r'angular [3]phenylene'}
for ax, mol, letter in ((axes[0], 'biphenylene', 'a'), (axes[1], 'angular3', 'b')):
    rr = [r for r in rnd if r['molecule'] == mol]
    ax.scatter([r['b3lyp_deformation_energy_eV'] for r in rr], [r['gap_decrease_eV'] for r in rr],
               s=9, color=MUTED, lw=0, zorder=2)
    for kind, colour, marker in (('P', BLUE, 'o'), ('BLA', ORANGE, 's')):
        L = ladder(mol, kind)
        ax.plot([r['b3lyp_deformation_energy_eV'] for r in L], [r['gap_decrease_eV'] for r in L],
                marker=marker, ms=3, color=colour, lw=1.3, zorder=3)
    ax.axhline(0, color='#bdbdbd', lw=0.5, zorder=1)
    ax.set_xlim(-0.1, 4.25); ax.set_xlabel('B3LYP deformation energy (eV)')
    ax.set_title(names[mol], pad=4); panel_letter(ax, letter, x=-0.2, y=1.03)
    ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
axes[0].set_ylabel('gap decrease (eV)')
ax = axes[2]
for mol, ls, marker in (('biphenylene', '-', 'o'), ('angular3', (0, (4, 2)), 's')):
    L = ladder(mol, 'P')
    ax.plot([r['b3lyp_deformation_energy_eV'] for r in L], [r['y0'] for r in L], marker=marker, ms=3,
            color=BLUE, lw=1.3, ls=ls)
ax.text(2.45, 0.80, 'biphenylene', ha='right', va='center')
ax.text(1.95, 0.47, r'angular [3]', ha='left', va='top')
ax.set_xlim(-0.1, 3.0); ax.set_xlabel('B3LYP deformation energy (eV)'); ax.set_ylabel(r'diradical character $y_0$')
ax.set_title('along the predicted pattern', pad=4); panel_letter(ax, 'c', x=-0.2, y=1.03)
ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
handles = [Line2D([], [], color=BLUE, marker='o', ms=3, lw=1.3, label='predicted pattern'),
           Line2D([], [], color=ORANGE, marker='s', ms=3, lw=1.3, label='bond-equalizing control'),
           Line2D([], [], color=MUTED, marker='o', ms=3, lw=0, label=r'random patterns (model energy 1.5 eV)')]
fig.legend(handles=handles, frameon=False, loc='lower left', bbox_to_anchor=(0.06, 0.0), ncol=3, handlelength=1.8)
fig.savefig('fig3_dft.pdf'); fig.savefig('fig3_dft.png', dpi=200)
plt.close(fig)
print('figures written')
