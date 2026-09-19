"""Initialize the requested two-state barrier form without using CH constants."""
from pathlib import Path
import json
import sys
import os
import argparse
for key in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'): os.environ[key]='1'
ROOT=Path(__file__).resolve().parent
deps=ROOT.parent/'.pecfit_stage/plot_deps'
if deps.exists(): sys.path.insert(0,str(deps))
import numpy as np
from scipy.optimize import least_squares
import coupled
from barrier import emo, potential, audit, components

BOHR=.529177210903
HARTREE=219474.63136320
C5=8.07*HARTREE*BOHR**5
C6=75.23*HARTREE*BOHR**6
LIMIT=55564.02471  # Absolute minimum-referenced dissociation limit supplied for AlF.

def block(state,kind,labels,p,fit,extra=''):
    header=f'poten {state}\nname "'+('X1Sigma+' if state==1 else 'A1Pi')+'"\n'
    header+=f'lambda {state-1}\nmult 1\n'+('symmetry +\n' if state==1 else '')
    header+=f'type {kind}\n'+extra+'values\n'
    return header+'\n'.join(f'{k:8s} {v:.14E}'+(' fit' if i in fit else '') for i,(k,v) in enumerate(zip(labels,p)))+'\nend'

def replace_block(text,state,new):
    import re
    return re.sub(rf'(?ims)^poten {state}\b.*?^end\s*$',lambda m:new,text,count=1)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input',type=Path,nargs='?',default=ROOT/'unique_levels.inp')
    ap.add_argument('--output',type=Path,default=ROOT)
    args=ap.parse_args(); out=args.output.resolve(); out.mkdir(parents=True,exist_ok=True)
    t=args.input.read_text(); fs=coupled.fields(t)
    names=['V0','RE','AE','RREF','PL','PR','NL','NR']+[f'B{i}' for i in range(7)]
    r=np.linspace(1.20,2.35,240)
    old=fs['POTEN',1,1].values
    p=np.r_[old[:4],6,6,6,6,old[8:],0,0]
    target=emo(r,old)
    def fun_x(b):
        q=p.copy(); q[8:]=b
        return (emo(r,q)-target)/(1+target/3000)
    res=least_squares(fun_x,p[8:],max_nfev=1000)
    p[8:]=res.x
    t=replace_block(t,1,block(1,'EMO',names,p,{1,*range(8,15)}))
    rx=np.linspace(1.35,1.98,220); original=emo(rx,fs['POTEN',2,2].values)
    gamma=8.; delta=.3
    rep=np.r_[LIMIT,1.5e7,delta,gamma,0,0,0,0,-C5,-C6,0,0]
    cp=np.array([3000,2.55,0,-1,6,6,4,4,.8,0,0,0,0.])
    base=np.array([43949.2,1.6485,80000,-1,6,6,6,6,1.72,.2,.2,0,0,0,0.])
    rb=4.8*BOHR; vb=LIMIT+.763*8065.544005
    def unpack(x):
        a=base.copy(); a[0]=x[0]*1000; a[1]=x[1]; a[8:]=x[2:9]
        b=rep.copy(); b[1]=np.exp(x[9])
        return np.r_[a,b,cp,1.]
    def val(r,p):
        a,b,c=components(r,p)
        return .5*(a+b-np.hypot(a-b,2*c))
    def fun_a(x):
        q=unpack(x)
        rr=np.array([rb-.001,rb,rb+.001]); vv=val(rr,q)
        return np.r_[(val(rx,q)-original)/(1+(original-43949.2)/2000),
                     (vv[1]-vb)/2,(vv[2]-vv[0])/.004,x[2:9]*.03]
    x0=np.r_[43.9492,1.6485,base[8:],np.log(rep[1])]
    lo=np.r_[43,1.60, .3,np.full(6,-10),np.log(1e5)]
    hi=np.r_[45,1.70, 4.,np.full(6,10),np.log(1e9)]
    res_a=least_squares(fun_a,x0,bounds=(lo,hi),max_nfev=5000,x_scale='jac',ftol=1e-12,gtol=1e-10,xtol=1e-12)
    q=unpack(res_a.x)
    labels=names+['V0','A','DELTA','GAMMA']+[f'B{i}' for i in range(1,9)]+['V0','RE','AE','RREF','PL','PR','NL','NR']+[f'B{i}' for i in range(5)]+['COMPON']
    t=replace_block(t,2,block(2,'coupled-pec',labels,q,{0,1,*range(8,15)},'sub-types EMO repulsive_exp EMO\nNparameters 15 12 13\n'))
    # Discard the copied AlCl-centred rotational correction from the seed only.
    f=coupled.fields(t)['BOBROT',1,1]
    edits=[(*span, f'{(1.654368 if i==0 else 0 if i>=4 else v):.14E}') for i,(span,v) in enumerate(zip(f.spans,f.values))]
    for start,stop,v in reversed(edits): t=t[:start]+v+t[stop:]
    import re
    m=re.search(r'(?ims)^bob-rot\b.*?^end\s*$',t)
    b=re.sub(r'\bfit\b','',m[0],flags=re.I)
    t=t[:m.start()]+b+t[m.end():]
    t='PRINT_PECS_AND_COUPLINGS_TO_FILE\n'+t
    (out/'seed.inp').write_text(t)
    info=dict(C5_au=8.07,C6_au=75.23,B5_cm_A5=-C5,B6_cm_A6=-C6,
              dissociation_absolute_cm=LIMIT,barrier_target_r_A=rb,barrier_target_above_limit_cm=vb-LIMIT,
              repulsive_gamma=gamma,shape=audit(t),inner_seed_max_error_cm=float(max(abs(val(rx,q)-original))),
              seed_optimization_success=bool(res_a.success),seed_optimization_message=res_a.message)
    (out/'seed.json').write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps(info,indent=2))

if __name__=='__main__': main()
