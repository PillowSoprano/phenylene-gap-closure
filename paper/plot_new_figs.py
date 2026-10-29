"""New figures for the revision: Kekulé structures of biphenylene (fig_kekule), bond correspondence of
zigzag and helical [4]phenylene (fig_corresp), and Fig. 3 without the diradical panel (fig3_dft)."""
import sys, os, json, csv
from pathlib import Path
HERE = Path(__file__).resolve().parent; os.chdir(HERE)
REPO = HERE.parent
sys.path.insert(0, str(REPO / 'src'))
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from gapclosure import phenylene, bipartite, kekule_terms
FIG_W = 6.5
plt.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'text.latex.preamble': r'\usepackage{amsmath}',
    'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'legend.fontsize': 7, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.major.size': 3, 'ytick.major.size': 3, 'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'lines.linewidth': 1.2})
INK, MUTED, GRID = '#1f2937', '#8c8c8c', '#e6e6e6'
RED, BLUE, ORANGE = '#e34948', '#2a78d6', '#eb6834'

def rotate_horizontal(P):
    P = P - P.mean(0); w, V = np.linalg.eigh(P.T @ P); ax = V[:, -1]
    R = np.array([[ax[0], ax[1]], [-ax[1], ax[0]]]); return P @ R.T

# ---------------- Kekulé structures of biphenylene ----------------
g = phenylene([]); P = rotate_horizontal(np.array(g['pos']))
T, L, R, bix = bipartite(g); terms = kekule_terms(T)
eix = {e: k for k, e in enumerate(bix)}
TE = [(s, frozenset(eix[(i, p[i])] for i in range(len(T)))) for s, p in terms]
Mm = [M for s, M in TE if s < 0][0]
others = sorted([M for s, M in TE if s > 0], key=lambda M: (len(M - Mm), sorted(M)))
ring_of = {}
for k, (u, v) in enumerate(g['edges']):
    hx = [r for r in g['rings'] if len(r) == 6 and u in r and v in r]
    ring_of[k] = hx[0] if hx else [r for r in g['rings'] if u in r and v in r][0]
def draw_kekule(ax, M, colour):
    for k, (u, v) in enumerate(g['edges']):
        ax.plot(*P[[u, v]].T, color=INK, lw=0.9, solid_capstyle='round', zorder=2)
        if k in M:
            c = P[ring_of[k]].mean(0); a, b = P[u], P[v]
            a2, b2 = a + 0.2*(c - a), b + 0.2*(c - b)        # inner line, shortened towards the ring centre
            ax.plot([a2[0], b2[0]], [a2[1], b2[1]], color=colour(k), lw=1.6, solid_capstyle='round', zorder=3)
            ax.plot(*P[[u, v]].T, color=colour(k), lw=1.6, solid_capstyle='round', zorder=3)
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_xlim(P[:, 0].min() - 0.25, P[:, 0].max() + 0.25); ax.set_ylim(P[:, 1].min() - 0.25, P[:, 1].max() + 0.25)
fig, axs = plt.subplots(1, 5, figsize=(FIG_W, 1.35), gridspec_kw=dict(left=0.005, right=0.995, bottom=0.03, top=0.80, wspace=0.12))
draw_kekule(axs[0], Mm, lambda k: BLUE)
axs[0].set_title(r'$M_-$\quad sign $-$', pad=3)
for a, M in enumerate(others, 1):
    draw_kekule(axs[a], M, lambda k, M=M: RED if k not in Mm else INK)
    axs[a].set_title(rf'$M_{a}$\quad sign $+$,\ $k_{a}={len(M - Mm)}$', pad=3)
fig.savefig('fig_kekule.pdf'); plt.close(fig)
print('k_a =', [len(M - Mm) for M in others])

# ---------------- zigzag / helical bond correspondence ----------------
C = json.load(open(REPO / 'results' / 'zigzag_helical_correspondence.json'))
changed = {k: v for k, v in C['map'].items() if set(C['zigzag_edges'][int(k)]) != set(C['helical_edges'][v])}
print('bonds whose atom labels change:', {tuple(C['zigzag_edges'][int(k)]): tuple(C['helical_edges'][v]) for k, v in changed.items()})
P0 = np.array(C['zigzag_pos']); rings = C['zig_rings']
c1, c3 = P0[rings[0]].mean(0), P0[[r for r in rings if 12 in r and 13 in r and len(r) == 6][0]].mean(0)
ang = np.arctan2(*(c3 - c1)[::-1]); Rm = np.array([[np.cos(-ang), -np.sin(-ang)], [np.sin(-ang), np.cos(-ang)]])
Pz, Ph = (np.array(C[f'{k}_pos']) @ Rm.T for k in ('zigzag', 'helical'))
allp = np.vstack([Pz, Ph]); lo, hi = allp.min(0) - 0.45, allp.max(0) + 0.45
fig, axs = plt.subplots(1, 2, figsize=(FIG_W, 2.7), gridspec_kw=dict(left=0.01, right=0.99, bottom=0.01, top=0.93, wspace=0.04))
cols = [RED, BLUE, ORANGE, '#16a34a']
pairs = sorted(changed.items(), key=lambda t: int(t[0]))
for ax, key, Pp, title in ((axs[0], 'zigzag', Pz, 'zigzag [4]phenylene'), (axs[1], 'helical', Ph, 'helical [4]phenylene')):
    E = C[f'{key}_edges']
    for k, (u, v) in enumerate(E):
        col, lw = '#9aa3ad', 1.4
        for c_i, (kz, kh) in enumerate(pairs):
            if (key == 'zigzag' and k == int(kz)) or (key == 'helical' and k == kh): col, lw = cols[c_i], 2.8
        ax.plot(*Pp[[u, v]].T, color=col, lw=lw, solid_capstyle='round', zorder=2)
    for i, (x, y) in enumerate(Pp):
        ax.scatter([x], [y], s=62, color='white', edgecolor=INK, lw=0.5, zorder=3)
        ax.text(x, y, str(i), fontsize=5.5, ha='center', va='center', zorder=4)
    ax.set_aspect('equal'); ax.axis('off'); ax.set_title(title, pad=2)
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1])
fig.savefig('fig_corresp.pdf'); plt.close(fig)

# ---------------- Fig. 3, two panels ----------------
def read(name):
    rows = list(csv.DictReader(open(f'fig_data/{name}')))
    for r in rows:
        for k, v in r.items():
            try: r[k] = float(v)
            except ValueError: pass
    return rows
lad, rnd = read('fig3_ladders.csv'), read('fig3_random.csv')
fig = plt.figure(figsize=(FIG_W, 2.35))
gs = fig.add_gridspec(1, 2, left=0.085, right=0.985, bottom=0.30, top=0.90, wspace=0.25)
axes = [fig.add_subplot(gs[0, i]) for i in range(2)]
for ax, mol, letter, nm in ((axes[0], 'biphenylene', 'a', 'biphenylene'), (axes[1], 'angular3', 'b', r'angular [3]phenylene')):
    rr = [r for r in rnd if r['molecule'] == mol]
    ax.scatter([r['b3lyp_deformation_energy_eV'] for r in rr], [r['gap_decrease_eV'] for r in rr], s=9, color=MUTED, lw=0, zorder=2)
    for kind, colour, marker in (('P', BLUE, 'o'), ('BLA', ORANGE, 's')):
        Lr = sorted([r for r in lad if r['molecule'] == mol and r['pattern'] == kind], key=lambda r: r['model_energy_eV'])
        ax.plot([r['b3lyp_deformation_energy_eV'] for r in Lr], [r['gap_decrease_eV'] for r in Lr], marker=marker, ms=3, color=colour, lw=1.3, zorder=3)
    ax.axhline(0, color='#bdbdbd', lw=0.5, zorder=1)
    ax.set_xlim(-0.1, 4.25); ax.set_xlabel('B3LYP deformation energy (eV)')
    ax.set_title(nm, pad=4); ax.text(-0.13, 1.03, rf'\textbf{{{letter}}}', transform=ax.transAxes, fontsize=9, va='bottom')
    ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
axes[0].set_ylabel('gap decrease (eV)')
handles = [Line2D([], [], color=BLUE, marker='o', ms=3, lw=1.3, label='predicted pattern $x^*$'),
           Line2D([], [], color=ORANGE, marker='s', ms=3, lw=1.3, label='bond-equalizing control'),
           Line2D([], [], color=MUTED, marker='o', ms=3, lw=0, label=r'random patterns (model energy 1.5 eV)')]
fig.legend(handles=handles, frameon=False, loc='lower left', bbox_to_anchor=(0.08, 0.0), ncol=3, handlelength=1.8)
fig.savefig('fig3_dft_2panel.pdf'); plt.close(fig)
print('done')
