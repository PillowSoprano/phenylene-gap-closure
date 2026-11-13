"""Chunked runner for within_class.py (background jobs do not survive here). Usage: N start stop"""
import sys, csv, json
src=open('within_class.py').read()
exec(src.split("cls={r")[0].replace("NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 9\n",""))
N,a,b=map(int,sys.argv[1:4])
cls={r['id']:r for r in csv.DictReader(open('../dataset_tree_phenylenes/tree_phenylenes_N1-12.csv'))}
rows=[r for r in csv.DictReader(open('../dataset_tree_phenylenes/tree_phenylene_molecules_N1-12.csv')) if int(r['N'])==N][a:b]
with open('within_class_mols.jsonl','a') as f:
    for r in rows:
        et,g0=hssh(phenylene_tree(json.loads(r['links'])))
        f.write(json.dumps(dict(id=r['id'],N=N,smiles=r['smiles'],rho=float(cls[r['id']]['rho']),gap_eq=float(r['huckel_gap']),overlaps=r['overlaps_in_plane']=='True',ET1=et,g0=g0))+'\n')
print(N,a,b,len(rows))
