"""Compare independent analytic curves with full-precision native Duo output."""
from pathlib import Path
import argparse
import sys
import json
deps=Path(__file__).resolve().parent.parent/'.pecfit_stage/plot_deps'
if deps.exists(): sys.path.insert(0,str(deps))
import numpy as np
from coupled import fields
from barrier import potential, audit

def check(directory):
    text=(directory/'input.inp').read_text(); fs=fields(text)
    table=np.loadtxt(directory/'Potential_functions.dat')
    r=np.linspace(table[0,0],table[-1,0],len(table))
    errors={str(s):float(max(abs(potential(r,fs['POTEN',s,s])-table[:,s]))) for s in (1,2)}
    report=dict(native_curve_max_difference_cm=errors,shape=audit(text))
    report['ok']=max(errors.values())<.001 and report['shape']['ok']
    (directory/'curve_check.json').write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('directory',type=Path)
    a=ap.parse_args(); result=check(a.directory)
    print(json.dumps(result,indent=2)); raise SystemExit(0 if result['ok'] else 1)
