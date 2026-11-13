"""Quick look at the table: extremes per N, and checks of the files."""
import csv
import json
import sys
from collections import defaultdict
from fractions import Fraction

import numpy as np

tag = sys.argv[1] if len(sys.argv) > 1 else "N1-12"
rows = list(csv.DictReader(open(f"tree_phenylenes_{tag}.csv")))
mols = list(csv.DictReader(open(f"tree_phenylene_molecules_{tag}.csv")))
pats = [json.loads(l) for l in open(f"tree_phenylene_patterns_{tag}.jsonl")]

byN = defaultdict(list)
for r in rows:
    byN[int(r["N"])].append(r)
print("N  classes molecules  rho_min (class)              rho_max (class)              single-minority")
for N in sorted(byN):
    R = byN[N]
    lo = min(R, key=lambda r: float(r["rho"]))
    hi = max(R, key=lambda r: float(r["rho"]))
    nm = sum(int(r["n_molecules"]) for r in R)
    single = [r["la_tree"] for r in R if r["minority"] == "1"]
    print(f"{N:2d} {len(R):6d} {nm:9d}  {float(lo['rho']):.6f} {lo['word'] or lo['la_tree'][:18]:20s} "
          f"{float(hi['rho']):.6f} {hi['word'] or hi['la_tree'][:18]:20s} {single}")

# consistency checks
bad = 0
for r in rows:
    lo, hi, rho = Fraction(r["rho_lower"]), Fraction(r["rho_upper"]), float(r["rho"])
    if not (float(lo) - 1e-15 <= rho <= float(hi) + 1e-15):
        bad += 1
    if (int(r["K"]) - int(r["ASC"])) % 2:
        bad += 1
print("bracket / parity violations:", bad)
print("molecule rows:", len(mols), " distinct SMILES:", len({m["smiles"] for m in mols}))
print("pattern lines:", len(pats))
# chains vs branched
for N in sorted(byN):
    ch = [float(r["rho"]) for r in byN[N] if r["word"] != "" or N <= 2]
    br = [float(r["rho"]) for r in byN[N] if r["word"] == "" and N > 2]
    if br:
        print(f"N={N}: chains {min(ch):.4f}-{max(ch):.4f}   branched {min(br):.4f}-{max(br):.4f}")
