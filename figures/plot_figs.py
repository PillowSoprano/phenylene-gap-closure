"""Draw all figures of the paper (Figs. 1-7) at the A5 text width of the MATCH template (321.5 pt = 4.45 in).
Figs. 2, 5 and 7 are drawn from fig_data/*.csv (written by export_figdata.py); Fig. 6 reads
results/zigzag_helical_correspondence.json; Figs. 1, 3 and 4 are drawn from the molecular graphs in src/gapclosure.py.
Text is rendered with LaTeX (text.usetex)."""
import sys, os, csv, json
from pathlib import Path
HERE = Path(__file__).resolve().parent; os.chdir(HERE)
REPO = HERE.parent                                   # repository root (src/, results/)
sys.path.insert(0, str(REPO / 'src'))
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from gapclosure import phenylene, bipartite, kekule_terms
FIG_W = 321.5/72.27
plt.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'text.latex.preamble': r'\usepackage{amsmath}',
    'font.size': 7, 'axes.titlesize': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
    'legend.fontsize': 6.5, 'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'xtick.minor.width': 0.4, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'axes.spines.top': False,
    'axes.spines.right': False, 'pdf.fonttype': 42, 'lines.linewidth': 1.1})
INK, MUTED, GRID = '#1f2937', '#8c8c8c', '#e6e6e6'
RED, BLUE, ORANGE, GREEN = '#e34948', '#2a78d6', '#eb6834', '#16a34a'
NEUTRAL_BOND = '#9aa3ad'
DIV = LinearSegmentedColormap.from_list('div', [BLUE, '#f0efec', RED])

def read(name):
    rows = list(csv.DictReader(open(f'fig_data/{name}')))
    for r in rows:
        for k, v in r.items():
            try: r[k] = float(v)
            except ValueError: pass
    return rows
def letter(ax, l, x, y): ax.text(x, y, rf'\textbf{{{l}}}', transform=ax.transAxes, fontsize=8, va='bottom', ha='left')
def draw_molecule(ax, atoms, bonds, colour_of, lw_of, dot=4):
    P = {int(a['atom']): (a['x'], a['y']) for a in atoms}
    for b in bonds:
        (x1, y1), (x2, y2) = P[int(b['atom1'])], P[int(b['atom2'])]
        ax.plot([x1, x2], [y1, y2], color=colour_of(b), lw=lw_of(b), solid_capstyle='round', zorder=2)
    xy = np.array(list(P.values())); ax.scatter(xy[:, 0], xy[:, 1], s=dot, color=INK, zorder=3, lw=0)
    ax.set_aspect('equal'); ax.axis('off'); return xy
def rotate_horizontal(P):
    P = P - P.mean(0); w, V = np.linalg.eigh(P.T @ P); a = V[:, -1]
    return P @ np.array([[a[0], a[1]], [-a[1], a[0]]]).T

# ---------------- Fig. 2: biphenylene, two routes (fig1_hook.pdf) ----------------
atoms, bonds = read('hook_atoms.csv'), read('hook_bonds.csv')
curves = read('hook_curves.csv'); rho = {r['quantity']: r['value'] for r in read('hook_summary.csv')}['rho_biphenylene']
fig = plt.figure(figsize=(FIG_W, 2.0))
gs = fig.add_gridspec(1, 3, width_ratios=[0.72, 0.72, 1.6], left=0.0, right=0.985, bottom=0.2, top=0.86, wspace=0.42)
axa, axb, axc = (fig.add_subplot(gs[0, i]) for i in range(3))
draw_molecule(axa, atoms, bonds, lambda b: RED if b['single_bond_route'] == 1 else NEUTRAL_BOND, lambda b: 2.6 if b['single_bond_route'] == 1 else 1.3)
draw_molecule(axb, atoms, bonds, lambda b: RED if b['cooperative_sign'] > 0 else BLUE, lambda b: 2.0)
for ax, l, t in ((axa, 'a', r'one bond: $\times4$'), (axb, 'b', rf'all bonds: $\pm{100*rho:.1f}\%$')):
    letter(ax, l, 0.0, 1.02); ax.text(0.17, 1.02, t, transform=ax.transAxes, va='bottom')
    ax.set_xlim(-1.35, 2.75); ax.set_ylim(-3.75, 1.35)
axb.legend(handles=[Line2D([], [], color=RED, lw=2, label='stronger'), Line2D([], [], color=BLUE, lw=2, label='weaker')],
           frameon=False, loc='upper center', bbox_to_anchor=(0.0, 0.02), ncol=2, handlelength=1.2, columnspacing=1.0, borderaxespad=0)
single = np.array([(r['largest_relative_change_percent'], r['huckel_gap_abs_beta']) for r in curves if r['route'] == 'single'])
coop = np.array([(r['largest_relative_change_percent'], r['huckel_gap_abs_beta']) for r in curves if r['route'] == 'cooperative'])
axc.plot(*single.T, color=MUTED, lw=1.2, ls=(0, (4, 2))); axc.plot(*coop.T, color=INK, lw=1.2)
axc.plot([single[-1, 0]], [single[-1, 1]], 'o', color=MUTED, ms=3); axc.plot([coop[-1, 0]], [coop[-1, 1]], 'o', color=INK, ms=3)
axc.text(60, 0.47, r'one bond (a)', color=MUTED, ha='right'); axc.text(9.5, 0.12, r'all bonds (b)', color=INK, ha='right')
axc.set_xscale('log'); axc.set_xlim(0.5, 420); axc.set_ylim(-0.02, 0.95)
axc.set_xticks([1, 3, 10, 30, 100, 300]); axc.set_xticklabels(['1', '3', '10', '30', '100', '300'])
axc.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
axc.set_xlabel(r'largest relative change of a coupling (\%)'); axc.set_ylabel(r'H\"uckel gap ($|\beta|$)', labelpad=2)
axc.yaxis.grid(True, color=GRID, lw=0.5); axc.set_axisbelow(True); letter(axc, 'c', -0.2, 1.02)
fig.savefig('fig1_hook.pdf'); plt.close(fig)

# ---------------- Fig. 3: Kekule structures of biphenylene (fig_kekule.pdf) ----------------
g = phenylene([]); P = rotate_horizontal(np.array(g['pos']))
T, L, R, bix = bipartite(g); terms = kekule_terms(T); eix = {e: k for k, e in enumerate(bix)}
TE = [(s, frozenset(eix[(i, p[i])] for i in range(len(T)))) for s, p in terms]
Mm = [M for s, M in TE if s < 0][0]
others = sorted([M for s, M in TE if s > 0], key=lambda M: (len(M - Mm), sorted(M)))
ring_of = {}
for k, (u, v) in enumerate(g['edges']):
    hx = [r for r in g['rings'] if len(r) == 6 and u in r and v in r]
    ring_of[k] = hx[0] if hx else [r for r in g['rings'] if u in r and v in r][0]
def draw_kekule(ax, M, colour):
    for k, (u, v) in enumerate(g['edges']):
        ax.plot(*P[[u, v]].T, color=INK, lw=0.7, solid_capstyle='round', zorder=2)
        if k in M:
            c = P[ring_of[k]].mean(0); a, b = P[u], P[v]; a2, b2 = a + 0.22*(c - a), b + 0.22*(c - b)
            ax.plot([a2[0], b2[0]], [a2[1], b2[1]], color=colour(k), lw=1.2, solid_capstyle='round', zorder=3)
            ax.plot(*P[[u, v]].T, color=colour(k), lw=1.2, solid_capstyle='round', zorder=3)
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_xlim(P[:, 0].min() - 0.2, P[:, 0].max() + 0.2); ax.set_ylim(P[:, 1].min() - 0.2, P[:, 1].max() + 0.2)
fig, axs = plt.subplots(1, 5, figsize=(FIG_W, 0.9), gridspec_kw=dict(left=0.005, right=0.995, bottom=0.02, top=0.8, wspace=0.1))
draw_kekule(axs[0], Mm, lambda k: BLUE); axs[0].set_title(r'$M_-$ ($-$)', pad=2)
for a, M in enumerate(others, 1):
    draw_kekule(axs[a], M, lambda k, M=M: RED if k not in Mm else INK)
    axs[a].set_title(rf'$M_{a}$ ($+$), $k_{a}={len(M - Mm)}$', pad=2)
fig.savefig('fig_kekule.pdf'); plt.close(fig)

# ---------------- Fig. 5: closing patterns, 3 x 2, molecules turned horizontal (fig2_closing.pdf) ----------------
keys = [('linear3', 'linear [3]'), ('angular3', 'angular [3]'), ('linear4', 'linear [4]'),
        ('bent4', 'bent [4]'), ('zigzag4', 'zigzag [4]'), ('triangular4', 'triangular [4]')]
rhos = {r['skeleton']: r['rho'] for r in read('closing_summary.csv')}
data = {k: (read(f'closing_{k}_atoms.csv'), read(f'closing_{k}_bonds.csv')) for k, _ in keys}
rot = {}
for k, _ in keys:
    xy = np.array([(a['x'], a['y']) for a in data[k][0]]); rot[k] = rotate_horizontal(xy)
pad = 0.55
spanx = max(np.ptp(rot[k][:, 0]) for k, _ in keys) + 2*pad          # common scale set by the widest molecule
colw = FIG_W/2; scale = colw/spanx                                   # inches per bond length
TOP = 0.16                                                           # inches for the panel title
rowh = [max(np.ptp(rot[k][:, 1]) for k, _ in keys[2*r:2*r+2])*scale + 2*pad*scale + TOP for r in range(3)]
leg_h = 0.62; H = sum(rowh) + leg_h
norm = TwoSlopeNorm(vcenter=0, vmin=-0.21, vmax=0.21)
WSCALE = 1100*(scale/0.27)**2                                       # marker area (pt^2) per unit weight; same size relative to bonds as before
fig = plt.figure(figsize=(FIG_W, H))
y0 = H
for r in range(3):
    y0 -= rowh[r]
    for cidx in range(2):
        k, label = keys[2*r + cidx]; atoms, bonds = data[k]
        ax = fig.add_axes([cidx*0.5, y0/H, 0.5, (rowh[r] - TOP)/H])
        P = rot[k]; idx = {int(a['atom']): n for n, a in enumerate(atoms)}
        for b in bonds:
            u, v = idx[int(b['atom1'])], idx[int(b['atom2'])]
            ax.plot(*P[[u, v]].T, color=DIV(norm(b['coupling_change'])), lw=1.6, solid_capstyle='round', zorder=2)
        ax.scatter(P[:, 0], P[:, 1], s=2.5, color=INK, zorder=3, lw=0)
        w = np.array([a['closing_orbital_weight'] for a in atoms])
        ax.scatter(P[:, 0], P[:, 1], s=WSCALE*w, facecolor='none', edgecolor=INK, lw=0.45, zorder=4)
        c = (P.max(0) + P.min(0))/2; hw = colw/scale/2; hh = (rowh[r] - TOP)/scale/2
        ax.set_xlim(c[0]-hw, c[0]+hw); ax.set_ylim(c[1]-hh, c[1]+hh); ax.set_aspect('equal'); ax.axis('off')
        fig.text(cidx*0.5 + 0.01, (y0 + rowh[r] - TOP + 0.02)/H, rf'\textbf{{{"abcdef"[2*r+cidx]}}}\ \ {label}, $\rho={rhos[k]:.3f}$', va='bottom')
cax = fig.add_axes([0.05, 0.36/H, 0.42, 0.08/H])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=DIV), cax=cax, orientation='horizontal')
cb.set_label('relative change of the coupling', labelpad=2); cb.outline.set_linewidth(0.4); cb.set_ticks([-0.2, -0.1, 0, 0.1, 0.2])
lax = fig.add_axes([0.57, 0.10/H, 0.42, 0.28/H]); lax.axis('off')
for n_, wv in enumerate((0.05, 0.10, 0.15)):
    lax.scatter([0.08 + 0.3*n_], [0.5], s=WSCALE*wv, facecolor='none', edgecolor=INK, lw=0.45, clip_on=False)
    lax.text(0.08 + 0.3*n_ + 0.07, 0.5, rf'${wv:.2f}$', va='center')
lax.set_xlim(0, 1); lax.set_ylim(0, 1); lax.text(0.0, 1.3, 'closing-orbital weight on an atom', va='bottom')
fig.savefig('fig2_closing.pdf'); plt.close(fig)
print('closing figure height', round(H, 2), 'in; scale', round(scale, 3), 'in per bond')

# ---------------- Fig. 6 (Appendix B): bond correspondence, stacked (fig_corresp.pdf) ----------------
C = json.load(open(REPO / 'results' / 'zigzag_helical_correspondence.json'))
changed = sorted(((int(k), v) for k, v in C['map'].items() if set(C['zigzag_edges'][int(k)]) != set(C['helical_edges'][v])))
P0 = np.array(C['zigzag_pos']); rings = C['zig_rings']
c1 = P0[rings[0]].mean(0); c3 = P0[[r for r in rings if 12 in r and 13 in r and len(r) == 6][0]].mean(0)
ang = np.arctan2(*(c3 - c1)[::-1]); Rm = np.array([[np.cos(-ang), -np.sin(-ang)], [np.sin(-ang), np.cos(-ang)]])
Pz, Ph = (np.array(C[f'{k}_pos']) @ Rm.T for k in ('zigzag', 'helical'))
cols = [RED, BLUE, ORANGE, GREEN]
spx = max(np.ptp(Pz[:, 0]), np.ptp(Ph[:, 0])) + 0.8; scl = FIG_W/spx
hz, hh_ = (np.ptp(Pz[:, 1]) + 0.8)*scl + 0.15, (np.ptp(Ph[:, 1]) + 0.8)*scl + 0.15
fig, axs = plt.subplots(2, 1, figsize=(FIG_W, hz + hh_), gridspec_kw=dict(left=0.0, right=1.0, bottom=0.0, top=1.0, hspace=0.0, height_ratios=[hz, hh_]))
for ax, key, Pp, title in ((axs[0], 'zigzag', Pz, 'zigzag [4]phenylene'), (axs[1], 'helical', Ph, 'helical [4]phenylene')):
    for k, (u, v) in enumerate(C[f'{key}_edges']):
        col, lw = NEUTRAL_BOND, 1.0
        for ci, (kz, kh) in enumerate(changed):
            if (key == 'zigzag' and k == kz) or (key == 'helical' and k == kh): col, lw = cols[ci], 2.2
        ax.plot(*Pp[[u, v]].T, color=col, lw=lw, solid_capstyle='round', zorder=2)
    for i, (x, y) in enumerate(Pp):
        ax.scatter([x], [y], s=70, color='white', edgecolor=INK, lw=0.45, zorder=3)
        ax.text(x, y, str(i), fontsize=5.5, ha='center', va='center', zorder=4)
    ax.set_aspect('equal'); ax.axis('off'); ax.text(0.0, 0.97, title, transform=ax.transAxes, va='top')
    xw = spx/2; xc = (Pz[:, 0].min() + Pz[:, 0].max())/2
    ax.set_xlim(xc - xw, xc + xw); ax.set_ylim(Pp[:, 1].min() - 0.4, Pp[:, 1].max() + 0.4 + 0.15/scl)
fig.savefig('fig_corresp.pdf'); plt.close(fig)

# ---------------- Fig. 7 (Appendix C): B3LYP efficiencies (fig_dft.pdf) ----------------
lad, rnd = read('dft_ladders.csv'), read('dft_random.csv')
fig = plt.figure(figsize=(FIG_W, 2.2))
gs = fig.add_gridspec(1, 2, left=0.1, right=0.99, bottom=0.36, top=0.9, wspace=0.3)
axes = [fig.add_subplot(gs[0, i]) for i in range(2)]
for ax, mol, l, nm in ((axes[0], 'biphenylene', 'a', 'biphenylene'), (axes[1], 'angular3', 'b', r'angular [3]phenylene')):
    rr = [r for r in rnd if r['molecule'] == mol]
    ax.scatter([r['b3lyp_deformation_energy_eV'] for r in rr], [r['gap_decrease_eV'] for r in rr], s=6, color=MUTED, lw=0, zorder=2)
    for kind, colour, marker in (('P', BLUE, 'o'), ('BLA', ORANGE, 's')):
        Lr = sorted([r for r in lad if r['molecule'] == mol and r['pattern'] == kind], key=lambda r: r['model_energy_eV'])
        ax.plot([r['b3lyp_deformation_energy_eV'] for r in Lr], [r['gap_decrease_eV'] for r in Lr], marker=marker, ms=2.5, color=colour, lw=1.1, zorder=3)
    ax.axhline(0, color='#bdbdbd', lw=0.5, zorder=1); ax.set_xlim(-0.1, 4.25); ax.set_xlabel('B3LYP deformation energy (eV)')
    ax.set_title(nm, pad=3); letter(ax, l, -0.2, 1.03); ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
axes[0].set_ylabel('gap decrease (eV)')
fig.legend(handles=[Line2D([], [], color=BLUE, marker='o', ms=2.5, lw=1.1, label='predicted pattern $x^*$'),
                    Line2D([], [], color=ORANGE, marker='s', ms=2.5, lw=1.1, label='bond-equalizing control'),
                    Line2D([], [], color=MUTED, marker='o', ms=2.5, lw=0, label='random patterns (model energy 1.5 eV)')],
           frameon=False, loc='lower left', bbox_to_anchor=(0.06, 0.0), ncol=2, handlelength=1.6, columnspacing=1.2)
fig.savefig('fig_dft.pdf'); plt.close(fig)
print('MATCH figures written')

# ---------------- Fig. 1: benzocyclobutadiene, graph and Kekule structures (fig_bcb_graph.pdf, fig_bcb_kekule.pdf) ----------------
from gapclosure import benzocyclobutadiene
gb = benzocyclobutadiene(); Pb = np.array(gb['pos']); Eb = [tuple(e) for e in gb['edges']]
Tb, Lb, Rb, bixb = bipartite(gb)
starred = set(Lb)
fig, ax = plt.subplots(figsize=(1.45, 1.45)); fig.subplots_adjust(0, 0, 1, 1)
for u, v in Eb: ax.plot(*Pb[[u, v]].T, color=INK, lw=1.0, zorder=1)
for i, (x, y) in enumerate(Pb):
    st = i in starred
    ax.scatter([x], [y], s=95, color=INK if st else 'white', edgecolor=INK, lw=0.7, zorder=2)
    ax.text(x, y, str(i + 1), fontsize=6.5, ha='center', va='center', color='white' if st else INK, zorder=3)
ax.set_aspect('equal'); ax.axis('off')
ax.set_xlim(Pb[:, 0].min() - 0.3, Pb[:, 0].max() + 0.3); ax.set_ylim(Pb[:, 1].min() - 0.3, Pb[:, 1].max() + 0.3)
fig.savefig('fig_bcb_graph.pdf'); plt.close(fig)
termsb = kekule_terms(Tb)
order = sorted(range(len(termsb)), key=lambda k: -termsb[k][0])
fig, axs = plt.subplots(1, 3, figsize=(FIG_W*0.92, 0.95), gridspec_kw=dict(left=0.0, right=1.0, bottom=0.0, top=0.82, wspace=0.25))
ring6 = [r for r in gb['rings'] if len(r) == 6][0]; ring4 = [r for r in gb['rings'] if len(r) == 4][0]
for ax, k in zip(axs, order):
    s, p = termsb[k]; M = {frozenset((Lb[i], Rb[p[i]])) for i in range(len(Tb))}
    for u, v in Eb:
        ax.plot(*Pb[[u, v]].T, color=INK, lw=0.7, zorder=1)
        if frozenset((u, v)) in M:
            ring = ring6 if (u in ring6 and v in ring6) else ring4
            c = Pb[ring].mean(0); a, b = Pb[u] + 0.22*(c - Pb[u]), Pb[v] + 0.22*(c - Pb[v])
            colr = BLUE if s < 0 else INK
            ax.plot([a[0], b[0]], [a[1], b[1]], color=colr, lw=1.2, zorder=2); ax.plot(*Pb[[u, v]].T, color=colr, lw=1.2, zorder=2)
    ax.set_title('sign $+$' if s > 0 else 'sign $-$', pad=2)
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_xlim(Pb[:, 0].min() - 0.15, Pb[:, 0].max() + 0.15); ax.set_ylim(Pb[:, 1].min() - 0.15, Pb[:, 1].max() + 0.15)
fig.savefig('fig_bcb_kekule.pdf'); plt.close(fig)
print('BCB terms (atoms 1-based):', [(s, sorted(tuple(sorted((Lb[i]+1, Rb[p[i]]+1))) for i in range(len(Tb)))) for s, p in [termsb[k] for k in order]])

# ---------------- Fig. 4: anatomy of a phenylene, L/A modes, gluing and regluing (fig_la.pdf) ----------------
from gapclosure import phenylene as _ph
gz, gh = _ph(['L', 'L']), _ph(['L', 'R'])
P0 = np.array(gz['pos']); c1 = P0[gz['rings'][0]].mean(0); c3 = P0[gz['rings'][4]].mean(0)
ang = np.arctan2(*(c3 - c1)[::-1]); Rm = np.array([[np.cos(-ang), -np.sin(-ang)], [np.sin(-ang), np.cos(-ang)]])
def pos(g): return np.array(g['pos']) @ Rm.T
def conn(g, k):                       # square S_k (k = 1..N-1): path [a, c, d, b]
    a_, c_, d_, b_ = g['rings'][2*k - 1]; return a_, b_, c_, d_
CONN = ORANGE; TAIL = BLUE; HEAD = '#6b7280'
def draw(ax, g, P, colour_of=lambda e: HEAD, lw_of=lambda e: 1.0):
    for e in g['edges']:
        ax.plot(*P[list(e)].T, color=colour_of(tuple(e)), lw=lw_of(tuple(e)), solid_capstyle='round', zorder=2)
    ax.set_aspect('equal'); ax.axis('off')
fig = plt.figure(figsize=(FIG_W, 3.35))
# (a) zigzag with ring labels and connecting bonds
Pz = pos(gz); conn_edges = set()
for k in (1, 2, 3):
    a_, b_, c_, d_ = conn(gz, k); conn_edges |= {frozenset((a_, c_)), frozenset((b_, d_))}
axa = fig.add_axes([0.0, 0.53, 0.62, 0.42])
draw(axa, gz, Pz, lambda e: CONN if frozenset(e) in conn_edges else HEAD, lambda e: 2.0 if frozenset(e) in conn_edges else 1.0)
hexes = [r for r in gz['rings'] if len(r) == 6]; sqs = [r for r in gz['rings'] if len(r) == 4]
for n, r in enumerate(hexes, 1): axa.text(*Pz[r].mean(0), rf'$H_{n}$', ha='center', va='center', fontsize=7)
for n, r in enumerate(sqs, 1): axa.text(*Pz[r].mean(0), rf'$S_{n}$', ha='center', va='center', fontsize=5.5)
axa.set_xlim(Pz[:, 0].min() - 0.3, Pz[:, 0].max() + 0.3); axa.set_ylim(Pz[:, 1].min() - 0.3, Pz[:, 1].max() + 0.3)
fig.text(0.0, 0.96, r'\textbf{a}\ \ zigzag [4]phenylene (L/A sequence AA)', va='bottom')
# (b) one internal hexagon: incoming bond h0h1 and the three allowed outgoing bonds
axb = fig.add_axes([0.66, 0.53, 0.34, 0.42])
angs = np.deg2rad([240, 300, 0, 60, 120, 180]); Hx = np.c_[np.cos(angs), np.sin(angs)]
sq = np.array([Hx[0], Hx[1], Hx[1] + [0, -1], Hx[0] + [0, -1]])
axb.plot(*np.r_[sq[1:], sq[:1]][[0, 1, 2, 3]].T, color=HEAD, lw=0.9, ls=(0, (2, 1.5)))
for i in range(6):
    j = (i + 1) % 6; col, lw = HEAD, 1.0
    if (i, j) == (3, 4): col, lw = GREEN, 2.2
    if (i, j) in ((2, 3), (4, 5)): col, lw = ORANGE, 2.2
    if (i, j) == (0, 1): col, lw = INK, 2.2
    axb.plot(*Hx[[i, j]].T, color=col, lw=lw, solid_capstyle='round')
for i, (x, y) in enumerate(Hx): axb.text(1.28*x, 1.22*y, rf'$h_{i}$', ha='center', va='center', fontsize=6.5)
axb.text(0, 1.18, 'L', color=GREEN, ha='center', va='bottom', fontsize=7)
axb.text(0.98, 0.62, 'A', color=ORANGE, ha='left', va='center', fontsize=7); axb.text(-0.98, 0.62, 'A', color=ORANGE, ha='right', va='center', fontsize=7)
axb.text(0, -1.5, r'$S_{k-1}$', ha='center', va='center', fontsize=5.5)
axb.set_aspect('equal'); axb.axis('off'); axb.set_xlim(-1.55, 1.55); axb.set_ylim(-2.1, 1.5)
fig.text(0.66, 0.96, r'\textbf{b}\ \ hexagon $H_k$', va='bottom')
# (c), (d) gluing: head A grey, tail B blue, connecting bonds orange, labels a b c d
def panel(ax, g, title_x, letter_, title):
    P = pos(g); a_, b_, c_, d_ = conn(g, 3)
    tail = set(gz['rings'][6]) if g is gz else set(gh['rings'][6])
    ce = {frozenset((a_, c_)), frozenset((b_, d_))}
    def colour(e):
        if frozenset(e) in ce: return CONN
        return TAIL if (e[0] in tail and e[1] in tail) else HEAD
    draw(ax, g, P, colour, lambda e: 2.0 if (frozenset(e) in ce or (e[0] in tail and e[1] in tail)) else 1.0)
    H3 = g['rings'][4]; alt = [(H3[1], H3[2]), (H3[3], H3[4])]
    for u, v in alt:
        if {u, v} != {a_, b_}: ax.plot(*P[[u, v]].T, color=INK, lw=1.0, ls=(0, (1.2, 1.2)), zorder=3)
    cen = P[H3].mean(0)
    for lab, at in (('a', a_), ('b', b_)):
        d = P[at] - cen; ax.text(*(P[at] - 0.42*d/np.linalg.norm(d)), rf'${lab}$', ha='center', va='center', fontsize=6.5)
    cB = P[list(tail)].mean(0)
    for lab, at in (('c', c_), ('d', d_)):
        d = P[at] - cB; ax.text(*(P[at] - 0.42*d/np.linalg.norm(d)), rf'${lab}$', ha='center', va='center', fontsize=6.5)
    ax.text(*P[gz['rings'][0]].mean(0) if False else P[g['rings'][2]].mean(0), r'$A$', ha='center', va='center', fontsize=7, color=HEAD)
    ax.text(*cB, r'$B$', ha='center', va='center', fontsize=7, color=TAIL)
    return P
axc = fig.add_axes([0.0, 0.0, 0.5, 0.48]); axd = fig.add_axes([0.5, 0.0, 0.5, 0.48])
Pc = panel(axc, gz, 0.0, 'c', ''); Pd = panel(axd, gh, 0.5, 'd', '')
allx = np.r_[Pc[:, 0], Pd[:, 0]]; ally = np.r_[Pc[:, 1], Pd[:, 1]]
wspan = np.ptp(Pc[:, 0]) + 0.8
for ax, P in ((axc, Pc), (axd, Pd)):
    xc = (P[:, 0].min() + P[:, 0].max())/2
    ax.set_xlim(xc - wspan/2, xc + wspan/2); ax.set_ylim(ally.min() - 0.3, ally.max() + 0.3)
fig.text(0.0, 0.485, r'\textbf{c}\ \ tail $B$ glued at $ab$ (zigzag)', va='bottom')
fig.text(0.5, 0.485, r'\textbf{d}\ \ reglued at the other A bond (helical)', va='bottom')
fig.savefig('fig_la.pdf'); plt.close(fig)
