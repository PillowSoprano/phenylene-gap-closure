"""Merge dft_isomers_<name>.json files from parallel runs into dft_isomers.json."""
import glob, json, os
here = os.path.dirname(os.path.abspath(__file__))
out = os.path.join(here, 'dft_isomers.json')
res = json.load(open(out)) if os.path.exists(out) else {}
for f in sorted(glob.glob(os.path.join(here, 'dft_isomers_*.json'))):
    if f.endswith('_test.json'): continue
    for k, v in json.load(open(f)).items():
        res.setdefault(k, {}).update(v)
json.dump(res, open(out, 'w'), indent=1)
for k, v in res.items():
    print(k, {s: round(v[s].get('ET1_eV', v[s]['E_hartree']), 3) for s in v if isinstance(v[s], dict)})
