"""Thresholds against independently computed values: the [N]phenylenes by Rohn's vertex criterion
(kekule_thresholds/results/thresholds.json) and the acene-cyclobutadiene-acene table."""
import pytest

import gapclosure as gc

ROHN = {  # exact rho by Rohn's vertex criterion, kekule_thresholds/results/thresholds.json
    "benzocyclobutadiene": 0.11972592295680506,
    "biphenylene": 0.1888289520324992,
    "linear[3]phenylene": 0.11972592295680512,
    "angular[3]phenylene": 0.1961605581327162,
    "linear[4]phenylene": 0.07558083370097893,
    "zigzag[4]phenylene": 0.16139548632918518,
    "helical[4]phenylene": 0.16139548632918505,
    "bent[4]phenylene": 0.12647161993651457,
    "triangular[4]phenylene": 0.20314202105540174,
}


@pytest.mark.parametrize("name,rho", ROHN.items())
def test_published(name, rho):
    assert gc.closure_threshold(name).rho == pytest.approx(rho, abs=1e-10)


@pytest.mark.parametrize("word", ["", "L", "A"])
def test_tree_equals_rohn(word):
    sk = gc.phenylene(word)
    a = gc.closure_threshold(sk, method="tree")
    b = gc.closure_threshold(sk, method="rohn")
    assert a.lower < b.rho <= a.upper
    assert a.rho == pytest.approx(b.rho, abs=1e-12)


ACENES = {(0, 1): 0.1197, (1, 1): 0.1888, (1, 2): 0.2070, (2, 2): 0.2218,
          (1, 4): 0.2153, (3, 3): 0.2312, (4, 4): 0.2343, (0, 4): 0.1674}


@pytest.mark.parametrize("units,rho", ACENES.items())
def test_acene_pairs(units, rho):
    r = gc.closure_threshold(gc.acene_chain(list(units)))
    assert r.minority == 1 and r.method.startswith("single")
    assert round(r.rho, 4) == rho


def test_biphenylene_equation_and_pattern():
    r = gc.closure_threshold("biphenylene")
    assert r.equation == "q^2 + 2q^4 + q^6 = 1"
    assert (r.K, r.ASC, r.minority) == (5, 3, 1)
    raised = [b for b, s in r.pattern.items() if s > 0]
    assert sorted(raised) == r.minority_structure and len(raised) == 6
    assert abs(sum(r.orbital) - 1) < 1e-9


def test_pattern_is_singular():
    import numpy as np
    for name in ["angular[3]phenylene", "triangular[4]phenylene"]:
        r = gc.closure_threshold(name)
        T, rows, cols, entry = r.skeleton.biadjacency()
        W = T.copy()
        for (i, j), b in zip(entry, r.skeleton.bonds):
            W[i, j] = 1 + r.pattern[b] * r.rho
        assert np.linalg.svd(W, compute_uv=False)[-1] < 1e-7


def test_no_minority():
    for sk in (gc.acene(1), gc.acene(3), gc.benzenoid([(0, 0), (1, 0), (1, 1)])):
        r = gc.closure_threshold(sk)
        assert r.rho == 1.0 and r.minority == 0


def test_closed_gap_and_errors():
    assert gc.closure_threshold(gc.annulene(8)).rho == 0.0      # cyclooctatetraene: det T = 0
    with pytest.raises(ValueError):
        gc.closure_threshold(gc.annulene(5))                     # odd ring


def test_disconnected():
    sk = gc.from_edges([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)][:0] + [(0, 1)] +
                       [(2 + u, 2 + v) for u, v in gc.biphenylene().bonds])
    r = gc.closure_threshold(sk)
    assert r.rho == pytest.approx(ROHN["biphenylene"], abs=1e-12)


def test_bounds_bracket_rohn():
    sk = gc.phenylene("A")
    exact = gc.closure_threshold(sk, method="rohn").rho
    r = gc.closure_threshold(sk, method="bounds")
    assert r.lower <= exact + 1e-12 and exact <= r.upper + 1e-12


def test_long_chain_is_fast():
    r = gc.closure_threshold(gc.phenylene("L" * 18), kekule_limit=1000)
    assert r.method.startswith("phenylene tree") and r.upper - r.lower < 1e-11


def test_recognize_relabelled():
    import random
    sk = gc.phenylene("LA")
    perm = list(range(sk.n)); random.Random(1).shuffle(perm)
    sk2 = gc.from_edges([(perm[u], perm[v]) for u, v in sk.bonds], n=sk.n)
    assert gc.closure_threshold(sk2).rho == pytest.approx(ROHN["bent[4]phenylene"], abs=1e-10)


def test_xyz_reader():
    sk = gc.biphenylene()
    xyz = f"{sk.n}\nbiphenylene\n" + "\n".join(f"C {1.4 * x:.5f} {1.4 * y:.5f} 0.0" for x, y in sk.pos)
    r = gc.closure_threshold(gc.from_xyz(xyz))
    assert r.rho == pytest.approx(ROHN["biphenylene"], abs=1e-12)


def test_smiles():
    pytest.importorskip("rdkit")
    r = gc.closure_threshold("c1ccc2c(c1)-c1ccccc1-2")
    assert r.rho == pytest.approx(ROHN["biphenylene"], abs=1e-12)
    r = gc.closure_threshold("C1=CC=C2C(=C1)C3=CC=C4C(=C23)C5=CC=CC=C54")
    assert 0 < r.rho < 1


def test_plot(tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    ax = gc.closure_threshold("angular[3]phenylene").plot(labels=True)
    ax.figure.savefig(tmp_path / "x.png")


def test_cli(capsys, tmp_path):
    from gapclosure.cli import main
    main(["biphenylene", "--plot", str(tmp_path / "b.png")])
    out = capsys.readouterr().out
    assert "0.18882" in out and (tmp_path / "b.png").exists()
    main(["--phenylene", "LA", "--json"])
    assert '"method": "phenylene tree (exact)"' in capsys.readouterr().out
