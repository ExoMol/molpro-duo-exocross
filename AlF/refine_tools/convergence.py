"""Independent grid, contraction and radial-box checks with native Duo."""
from pathlib import Path
import argparse
import csv
import json
import re
import subprocess
import sys
import coupled

def variant(text, points, rmax, vmax):
    text=re.sub(r'(?im)^\s*npoints[^\n]*',f'  npoints {points}',text,count=1)
    text=re.sub(r'(?im)^\s*range[^\n]*',f'  range .7 {rmax}',text,count=1)
    text=re.sub(r'(?im)^\s*vmax[^\n]*',f'  vmax {vmax} {vmax}',text,count=1)
    return text

def read(directory):
    with (directory/'residuals.csv').open() as f: return list(csv.DictReader(f))

def compare(a,b):
    if len(a)!=len(b): raise ValueError('Observation coverage changed')
    diffs=[]
    for r,s in zip(a,b):
        if any(r[k]!=s[k] for k in ('J','parity','state','v','observed_cm','input_weight')):
            raise ValueError('Observation identity changed')
        if float(r['input_weight'])>0:
            diffs.append((float(r['input_weight']),float(s['calculated_cm'])-float(r['calculated_cm'])))
    return dict(max_abs_cm=max(abs(d) for w,d in diffs),weighted_rms_cm=(sum(w*d*d for w,d in diffs)/sum(w for w,d in diffs))**.5)

def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('input',type=Path)
    ap.add_argument('--output',type=Path,required=True); ap.add_argument('--duo',type=Path,required=True)
    a=ap.parse_args(); a.output.mkdir(parents=True,exist_ok=False); text=a.input.read_text()
    cases=[('grid601',601,6,60),('grid801',801,6,60),('basis90',801,6,90),('box8',1001,8,100)]
    result={}; rows={}
    for name,n,rmax,vmax in cases:
        inp=a.output/f'{name}.inp'; inp.write_text(variant(text,n,rmax,vmax))
        cmd=[sys.executable,'-B',str(Path(__file__).with_name('run_native.py')),str(inp),str(a.output/name),'--duo',str(a.duo),'--jmax','98','--precise','--timeout','600']
        proc=subprocess.run(cmd,capture_output=True,text=True)
        if proc.returncode: raise RuntimeError(proc.stdout+proc.stderr)
        metrics=coupled.save_assessment(a.output/name,98); rows[name]=read(a.output/name)
        result[name]=dict(grid_points=n,rmax_A=rmax,vmax_each_state=vmax,metrics=metrics)
        print(name,metrics['weighted_rms_cm'],flush=True)
    checks={name:compare(rows[left],rows[right]) for name,left,right in [
        ('grid','grid601','grid801'),('contraction','grid801','basis90'),('box_and_final_basis','basis90','box8')]}
    accepted=all(c['max_abs_cm']<.0001 for c in checks.values()) and all(x['metrics']['assignment_mismatches']==0 for x in result.values())
    report=dict(status='numerically_validated' if accepted else 'needs_review',tolerance_cm=.0001,comparisons=checks,runs=result,
                selected_calculation=str((a.output/'box8').resolve()))
    (a.output/'convergence.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
