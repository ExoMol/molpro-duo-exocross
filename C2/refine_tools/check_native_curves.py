"""Cross-check analytic PEC diagnostics against native Duo's ab initio table."""
import argparse
import json
from pathlib import Path
import coupled
import curves

ap = argparse.ArgumentParser()
ap.add_argument('directory', type=Path)
a = ap.parse_args()
fields = {f.name: f for k,f in coupled.fields((a.directory/'input.inp').read_text()).items() if k[0]=='POTEN'}
current = None
errors = {name: [] for name in fields}
for line in (a.directory/'fit.pot').read_text().splitlines():
    p = line.split()
    if line.strip() in fields:
        current = fields[line.strip()]
        continue
    if p and not p[0].isdigit():
        current = None
    if current and len(p)==6 and p[0].isdigit():
        r, native = coupled.number(p[1]), coupled.number(p[3])
        errors[current.name].append(abs(curves.value(current.kind,current.values,r)-native))
summary = {name:dict(points=len(err), max_delta_cm=max(err) if err else None) for name,err in errors.items()}
print(json.dumps(summary, indent=2))
assert all(item['points'] and item['max_delta_cm'] < .006 for item in summary.values())
(a.directory/'analytic_native_curve_check.json').write_text(json.dumps(summary, indent=2)+'\n')
