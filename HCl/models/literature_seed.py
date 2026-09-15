"""Numerical MLJ seeds, not published MLJ coefficients.
RKR points quoted in Table I of https://arxiv.org/abs/2208.07895, attributed
there to Coxon & Hajigeorgiou (2015). Radii in bohr, energies in Hartree.
"""
import json,math
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
import numpy as np
from pecfit.forms import potential
from pecfit.model import shape_check

data="""1.717716 8.825021 -.000082
1.719071 6.585696 -.001031
1.722427 5.855864 -.003364
1.727531 5.437063 -.006866
1.734145 5.134002 -.011318
1.742171 4.890604 -.016584
1.751549 4.683618 -.022571
1.762308 4.500867 -.029214
1.774500 4.335109 -.036461
1.788227 4.181679 -.044276
1.803636 4.037339 -.052627
1.820931 3.899720 -.061495
1.840391 3.766980 -.070860
1.862397 3.637572 -.080712
1.887477 3.510072 -.091041
1.916400 3.383033 -.101842
1.950315 3.254758 -.113112
1.991106 3.122935 -.124850
2.042201 2.983789 -.137055
2.111183 2.829445 -.149731
2.224536 2.633188 -.162880"""
bohr=.529177210903; hartree=219474.6313708
rr=[]; energies=[]
for line in data.splitlines():
    a,b,e=map(float,line.split());rr.extend([a*bohr,b*bohr]);energies.extend([(e+.169641)*hartree]*2)
r=np.array(rr);ev=np.array(energies)
re=2.408544*bohr; de=.169641*hartree
z=2*(r-re)/(r+re)
sign=np.sign(r-re)
phi=(-np.log1p(-sign*np.sqrt(ev/de))+6*np.log(re/r))/z
phiinf=.5*math.log(2*de*re**6/(23.41*hartree*bohr**6))
grid=np.linspace(.5,12,4001)
mu=1.00782503223*34.968852682/(1.00782503223+34.968852682)
accepted=[]
for degree in (4,6,8,10):
    for r12 in (2.5,3.,3.5,4.):
        for delta in (2.,3.,4.):
            sw=1/(1+np.exp(delta*(r-r12)))
            design=np.polynomial.polynomial.polyvander(z,degree)*sw[:,None]
            coeff=np.linalg.lstsq(design,phi-(1-sw)*phiinf,rcond=None)[0]
            p={"TE":0.,"RE":re,"AE":de,"R12":r12,"DELTA":delta,"PHIINF":phiinf,"N":6.,**{f"PHI{i}":float(c) for i,c in enumerate(coeff)}}
            if not shape_check(p,(.5,12),max_potential=1e8)["ok"]:continue
            vs=np.array([potential(float(x),p) for x in grid])
            maxturns=0
            for j in (0,10,20,30,40,41):
                eff=vs+16.857629206/mu*j*(j+1)/grid**2
                signs=np.sign(np.diff(eff)); maxturns=max(maxturns,int(np.sum(np.diff(signs)!=0)))
            if maxturns>2:continue
            error=np.array([potential(float(x),p) for x in r])-ev
            rms=float(np.sqrt(np.mean(error*error)))
            accepted.append((rms,p))
            print(degree,r12,delta,rms,flush=True)
accepted.sort(key=lambda x:x[0])
(ROOT/"HCl/pecfit_runs").mkdir(parents=True, exist_ok=True)
(ROOT/"HCl/pecfit_runs/literature_mlj_seed_candidates.json").write_text(json.dumps([{"rms_pec_cm":score,"parameters":p} for score,p in accepted],indent=2))
print("best",accepted[:1])
