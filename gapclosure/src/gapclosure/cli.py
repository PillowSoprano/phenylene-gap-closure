"""Command line: gapclosure MOLECULE [--plot FILE] [--json]"""
from __future__ import annotations

import argparse
import sys


def _load(args):
    from . import skeleton as S
    if args.xyz:
        return S.from_xyz(args.xyz)
    if args.molfile:
        return S.from_molfile(args.molfile)
    if args.phenylene is not None:
        return S.phenylene(args.phenylene, helical=args.helical)
    if args.acenes:
        return S.acene_chain([int(x) for x in args.acenes.split(",")])
    if args.molecule is None:
        raise SystemExit("give a molecule: a name, a SMILES string, --phenylene WORD, "
                         "--acenes h1,h2,..., --xyz FILE or --molfile FILE")
    from .threshold import resolve
    return resolve(args.molecule)


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="gapclosure",
        description="Hueckel gap-closure threshold rho, closing pattern and closing orbital "
                    "of a conjugated hydrocarbon.")
    p.add_argument("molecule", nargs="?", help="name (e.g. biphenylene) or SMILES (needs RDKit)")
    p.add_argument("--phenylene", metavar="WORD", help="[N]phenylene chain from its L/A word ('' = biphenylene)")
    p.add_argument("--helical", action="store_true", help="draw angular turns on one side")
    p.add_argument("--acenes", metavar="H1,H2,...", help="acenes joined by four-membered rings (0 = single bond)")
    p.add_argument("--xyz", metavar="FILE", help="XYZ geometry")
    p.add_argument("--molfile", metavar="FILE", help="molfile / SDF (needs RDKit)")
    p.add_argument("--plot", metavar="FILE", help="save the closing-pattern figure (png, pdf, svg)")
    p.add_argument("--labels", action="store_true", help="number the atoms in the figure")
    p.add_argument("--json", action="store_true", help="print the full result as JSON")
    p.add_argument("--bits", type=int, default=40, help="bisection steps for phenylene trees")
    p.add_argument("--exact-limit", type=int, default=10, help="largest atoms/2 for Rohn's method")
    p.add_argument("--method", default="auto", choices=["auto", "theorem-a", "tree", "rohn", "bounds"])
    p.add_argument("--list", action="store_true", help="list the built-in names")
    args = p.parse_args(argv)

    if args.list:
        from .skeleton import NAMED
        print("\n".join(NAMED))
        return 0
    from .threshold import closure_threshold
    try:
        sk = _load(args)
        r = closure_threshold(sk, bits=args.bits, exact_limit=args.exact_limit, method=args.method)
    except (ValueError, ImportError, KeyError, OSError) as exc:
        print(f"gapclosure: {exc}", file=sys.stderr)
        return 1
    print(r.to_json(indent=1) if args.json else r.summary())
    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        ax = r.plot(labels=args.labels)
        ax.figure.savefig(args.plot, dpi=200, bbox_inches="tight")
        plt.close(ax.figure)
        if not args.json:
            print(f"  figure: {args.plot}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
