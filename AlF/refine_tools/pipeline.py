"""Automatic AlF preparation, constrained refinement, validation and plots."""
from pathlib import Path
import argparse
import json
import re
import subprocess
import sys
from configure import configure

TOOLS=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input',type=Path); ap.add_argument('--duo',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--ordinary-rounds',type=int,default=3)
    ap.add_argument('--robust-rounds',type=int,default=6)
    a=ap.parse_args(); root=a.output.resolve(); root.mkdir(parents=True,exist_ok=False)
    (root/'runs').mkdir(); a.input=a.input.resolve(); a.duo=a.duo.resolve()
    def run(script,*args):
        subprocess.run([sys.executable,'-B',str(TOOLS/script),*map(str,args)],check=True)
    def evaluate(inp,out):
        run('run_native.py',inp,out,'--duo',a.duo,'--jmax',98,'--precise')
        run('coupled.py',out,'--jmax',98)
    run('audit_input.py',a.input,'--output',root)
    run('seed_model.py',root/'unique_levels.inp','--output',root)
    evaluate(root/'unique_levels.inp',root/'runs/unique_baseline')
    evaluate(root/'seed.inp',root/'runs/barrier_seed')
    initial=root/'runs/ordinary'
    run('refine.py',root/'seed.inp','--duo',a.duo,'--output',initial,'--jmax',98,
        '--rounds',a.ordinary_rounds,'--proposal-steps',3,'--scale',1,'--robust',0)
    text=configure((initial/'best.inp').read_text(),fit=['L+:2:1:4','BOBROT:1:1:5','BOBROT:1:1:6'],values=['L+:2:1:0=1.6485'])
    match=re.search(r'(?ims)^poten 1\b.*?^end',text)
    block=re.sub(r'(?im)^N[LR]\s+.*',lambda m:m[0][:2]+' 8',match[0])
    block=block[:-3]+'B7 0 fit\nB8 0 fit\nend'
    seed=root/'rotation_x8.inp'; seed.write_text(text[:match.start()]+block+text[match.end():])
    robust=root/'runs/robust'
    run('refine.py',seed,'--duo',a.duo,'--output',robust,'--jmax',98,'--rounds',a.robust_rounds,
        '--proposal-steps',3,'--scale',1,'--robust',.00001,'--criterion','cauchy')
    run('convergence.py',robust/'best.inp','--duo',a.duo,'--output',root/'convergence')
    check=json.loads((root/'convergence/convergence.json').read_text())
    selected=Path(check['selected_calculation'])
    run('check_curves.py',selected)
    ordinary=Path(json.loads((initial/'result.json').read_text())['selected_calculation'])
    run('make_report.py',selected,'--workdir',root,'--ordinary',ordinary,'--output',root/'report')
    (root/'final.inp').write_text((selected/'input.inp').read_text())
    result=dict(numerical_status=check['status'],scientific_status='needs_review',
                reason='Review flagged experimental discrepancies and theoretical long-range assumptions.',
                model=str(root/'final.inp'),selected_calculation=str(selected))
    (root/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
