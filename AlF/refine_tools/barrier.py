"""Vectorized Duo curves and explicit AlF shape constraints."""
import numpy as np
import re as regex
from coupled import fields

def emo(r, p):
    te, re, ae, ref, pl, pr, nl, nr = p[:8]
    ref = re if ref <= 0 else ref
    power = np.where(r<=ref, pl, pr)
    y = (r**power-ref**power)/(r**power+ref**power)
    beta = sum(p[8+k]*y**k*np.where(r<=ref, k<=nl, k<=nr) for k in range(len(p)-8))
    return te+(ae-te)*np.expm1(-beta*(r-re))**2

def repulsive(r,p):
    return p[0]+p[1]*np.exp(-p[2]/r)/r**p[3]+sum(b/r**(k+1) for k,b in enumerate(p[4:]))

def components(r,p):
    n1=int(p[7])+9
    n2=12
    a=emo(r,p[:n1]); b=repulsive(r,p[n1:n1+n2]); c=emo(r,p[n1+n2:-1])
    return a,b,c

def potential(r,f):
    if f.kind=='EMO': return emo(r,f.values)
    if f.kind=='COUPLED-PEC':
        a,b,c=components(r,f.values)
        return .5*(a+b+(-1 if int(f.values[-1])==1 else 1)*np.hypot(a-b,2*c))
    raise ValueError(f.kind)

def bob(r,p):
    re,beta,gamma,power=p[:4]
    d=r-re; z=d*np.exp(-beta*d*d-gamma*d**4); y=(r**power-re**power)/(r**power+re**power)
    return (1-y)*np.polynomial.polynomial.polyval(z,p[4:-1])+y*p[-1]

def extrema(r,v,tol=1e-7):
    d=np.diff(v)
    return [(float(r[i+1]),float(v[i+1]),'min' if d[i]<0 else 'max') for i in range(len(d)-1) if d[i]*d[i+1]<0 and max(abs(d[i]),abs(d[i+1]))>tol]

def audit(text):
    fs=fields(text); r=np.unique(np.r_[np.linspace(.7,8,5001),np.geomspace(8,1000,1001)])
    # This implementation deliberately supports the requested three sub-functions.
    # Fail closed if another layout is supplied instead of guessing the slices.
    for key,f in fs.items():
        if f.kind=='COUPLED-PEC':
            block=regex.search(rf'(?ims)^poten {key[1]}\b.*?^end\s*$',text)[0]
            types=regex.search(r'(?im)^sub-types\s+([^\n]+)',block)
            sizes=regex.search(r'(?im)^Nparameters\s+([^\n]+)',block)
            expected=[int(f.values[7])+9,12,13]
            if not types or types[1].upper().split()!=['EMO','REPULSIVE_EXP','EMO'] or not sizes or list(map(int,sizes[1].split()))!=expected or sum(expected)+1!=len(f.values):
                return dict(ok=False,error='Unsupported coupled-pec layout')
    report={}; ok=True
    for state in (1,2):
        f=fs['POTEN',state,state]
        with np.errstate(over='ignore',invalid='ignore'):
            v=potential(r,f)
        ex=extrema(r,v)
        good=bool(np.all(np.isfinite(v)) and np.max(v)<1e9)
        inner=[e for e in ex if e[2]=='min' and e[0]<2.2]
        good=good and len(inner)==1 and 1.55<inner[0][0]<1.75
        if state==1:
            good=good and len(ex)==1 and abs(min(v))<1 and f.values[4]==f.values[5] and f.values[4] in (6,8)
        elif f.kind=='COUPLED-PEC':
            peaks=[e for e in ex if e[2]=='max']
            n1=int(f.values[7])+9; limit=f.values[n1]
            good=good and len(peaks)==1 and 2.25<peaks[0][0]<2.9 and 2500<peaks[0][1]-limit<9000
            good=good and len(ex)<=3 and all(e[1]>limit-20 for e in ex if e[2]=='min' and e[0]>2.2)
            good=good and f.values[4]==f.values[5] and f.values[4] in (6,8)
            coupling=f.values[n1+12:-1]
            good=good and coupling[4]==coupling[5] and coupling[4] in (6,8)
        report[str(state)]=dict(ok=bool(good), extrema=ex, value_at_1000=float(v[-1]))
        ok=ok and good
    for key,f in fs.items():
        if key[0]=='BOBROT':
            vals=bob(r,f.values); good=bool(max(abs(vals))<.02)
            report[str(key)]=dict(ok=good,max_abs=float(max(abs(vals))))
            ok=ok and good
        if key[0]=='L+':
            re,_,power,_=f.values[:4]; y=(r**power-re**power)/(r**power+re**power)
            vals=(1-y)*np.polynomial.polynomial.polyval(y,f.values[4:-1])+y*f.values[-1]
            good=bool(max(abs(vals))<3)
            report[str(key)]=dict(ok=good,max_abs_analytic=float(max(abs(vals))))
            ok=ok and good
    return dict(ok=bool(ok),range_angstrom=[.7,1000],fields=report)
