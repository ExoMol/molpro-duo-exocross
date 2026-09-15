"""Optional single-state sinc-DVR refinement with localized-state assignments.

Requires NumPy/SciPy. The primitive-grid Hamiltonian follows Duo's SINC
implementation and CODATA-2010 constants. Final energies are independently
recomputed by the supplied Duo executable. Quasibound levels are *not* certified
resonance poles: radial-box sensitivity remains an explicit validation failure.
"""
from __future__ import annotations

from dataclasses import replace
import json
import math
from pathlib import Path
import re
import time

import numpy as np
from scipy.linalg import eigh

from .duo import Runner, assess_energies, digest, parse_en, write_csv, write_json
from .forms import asymptote, bob_keys, bob_shape, centrifugal_shape, kind
from .model import FitError, render, shape_check
from .pipeline import load

# Exact constants in Duo/accuracy.f90 (CODATA 2010); native agreement is checked.
ASTON = 6.62606957e-27/(8*math.pi**2*2.99792458e10*1.660538921e-24*1e-16)


def curves(r, p):
    """Vectorized forms, with the same conventions as forms.potential/bob."""
    r = np.asarray(r, dtype=float)
    with np.errstate(over='raise', invalid='raise', divide='raise'):
        if kind(p) == 'EMO':
            ref = p['RREF'] if p['RREF'] > 0 else p['RE']
            power = np.where(r <= ref, p['PL'], p['PR'])
            y = (r**power-ref**power)/(r**power+ref**power)
            beta = sum(p[f'B{i}']*y**i*np.where(r <= ref, i <= p['NL'], i <= p['NR'])
                       for i in range(int(max(p['NL'],p['NR']))+1))
            v = p['V0']+(p['DE']-p['V0'])*np.expm1(-beta*(r-p['RE']))**2
        else:
            z = 2*(r-p['RE'])/(r+p['RE'])
            switch = 1/(1+np.exp(np.clip(p['DELTA']*(r-p['R12']),-700,700)))
            coeff = [p[f'PHI{i}'] for i in range(sum(bool(re.fullmatch(r'PHI\d+',k)) for k in p))]
            phi = switch*np.polynomial.polynomial.polyval(z,coeff)+(1-switch)*p['PHIINF']
            v = p['TE']+(p['AE']-p['TE'])*np.expm1(p['N']*np.log(p['RE']/r)-phi*z)**2
        q = np.zeros_like(r)
        if bob_keys(p):
            d = r-p['BR_RE']
            z = d*np.exp(-p['BR_BETA']*d*d-p['BR_GAMMA']*d**4)
            y = (r**p['BR_P']-p['BR_RE']**p['BR_P'])/(r**p['BR_P']+p['BR_RE']**p['BR_P'])
            coeff = [p[f'BR_B{i}'] for i in range(sum(bool(re.fullmatch(r'BR_B\d+',k)) for k in p))]
            q = (1-y)*np.polynomial.polynomial.polyval(z,coeff)+y*p['BR_BINF']
    return v,q


def select_localized(energies, vectors, effective, radii, dissociation, needed,
                     minimum_probability=0.5):
    """Assign inner-well node counts, choosing the most localized representative.

    A continuum state is never assigned solely because its energy is close to
    an observation. Above dissociation a candidate must lie below a resolved
    centrifugal barrier and have sufficient probability inside that barrier.
    """
    slopes = np.diff(effective)
    maxima = np.where((slopes[:-1] > 0)&(slopes[1:] < 0))[0]+1
    barrier = int(maxima[-1]) if len(maxima) else len(radii)-1
    probabilities = np.sum(vectors[:barrier+1]**2,axis=0)
    candidates = np.where((energies < dissociation) |
                          (bool(len(maxima))&(energies < effective[barrier])&
                           (probabilities >= minimum_probability)))[0]
    assigned = {}
    for col in candidates:
        wave = vectors[:barrier+1,col]
        # Remove only the numerical tails when counting inner-well nodes.
        wave = wave[np.abs(wave) > np.max(np.abs(wave))*1e-5]
        nodes = int(np.sum(wave[1:]*wave[:-1] < 0))
        if nodes not in assigned or probabilities[col] > probabilities[assigned[nodes]]:
            assigned[nodes] = int(col)
    missing = set(needed)-set(assigned)
    if missing:
        raise FitError(f'No unambiguous localized levels for v={sorted(missing)}')
    cols = np.array([assigned[v] for v in needed],dtype=int)
    info = [{'v':int(v),'rank':int(col+1),'inner_probability':float(probabilities[col]),
             'quasibound':bool(energies[col] >= dissociation),
             'barrier_r_angstrom':float(radii[barrier]) if len(maxima) else None,
             'barrier_absolute_cm':float(effective[barrier]) if len(maxima) else None}
            for v,col in zip(needed,cols)]
    return cols,info


class Spectrum:
    def __init__(self, model, grid, minimum_probability=0.5, extra_roots=40):
        if not model.masses:
            raise FitError('Localized refinement requires explicit isotopic masses')
        masses = [float(x) for x in model.masses.replace(',',' ').split()]
        if len(masses)!=2 or min(masses)<=0:
            raise FitError('Expected two positive isotopic masses')
        self.rot = ASTON*(masses[0]+masses[1])/(masses[0]*masses[1])
        self.grid = tuple(grid)
        lo,hi,n = grid
        if not 0<lo<hi or int(n)!=n or n<3:
            raise FitError('Invalid sinc-DVR grid')
        self.r = np.linspace(lo,hi,n)
        delta = np.arange(n)[:,None]-np.arange(n)
        self.kinetic = 2.*(-1.)**delta/np.maximum(delta*delta,1)
        np.fill_diagonal(self.kinetic,np.pi**2/3)
        self.kinetic *= self.rot/(self.r[1]-self.r[0])**2
        self.observations = sorted(model.observations,key=lambda o:o.key)
        self.byj = {j:[o for o in self.observations if o.j==j]
                    for j in sorted({0}|{o.j for o in self.observations})}
        self.nroots = min(n,max(o.v for o in self.observations)+extra_roots+1)
        self.minimum_probability = minimum_probability

    def evaluate(self,p,derivatives=None):
        v,q = curves(self.r,p)
        energies,metadata,jac = [],[],[]
        zpe,dzpe = None,None
        for j,observations in self.byj.items():
            centrifugal = self.rot*j*(j+1)/self.r**2
            effective = v+centrifugal*(1+q)
            vals,vec = eigh(self.kinetic+np.diag(effective),subset_by_index=(0,self.nroots-1),
                            check_finite=False,driver='evr')
            if j==0:
                zpe = float(vals[0])
                if derivatives is not None:dzpe = vec[:,0]**2@derivatives[0]
            if not observations:continue
            try:
                cols,info = select_localized(vals,vec,effective,self.r,asymptote(p),
                                             [o.v for o in observations],self.minimum_probability)
            except FitError as exc:
                raise FitError(f'J={j}: {exc}') from exc
            energies.extend(vals[cols]-zpe)
            metadata.extend(dict(row,J=j) for row in info)
            if derivatives is not None:
                dv,dq = derivatives
                jac.extend((vec[:,cols]**2).T@(dv+centrifugal[:,None]*dq)-dzpe)
        return np.asarray(energies),metadata,np.asarray(jac) if derivatives is not None else None,zpe


def physical(model,p,config):
    replace(model,parameters=p).validate()
    limits = config['validation']
    grid = limits['shape_range_angstrom'] or model.grid
    pec = shape_check(p,grid,limits['shape_tolerance_cm'],max_potential=limits.get('potential_ceiling_cm'))
    br = bob_shape(p,grid,max_abs=limits['bob_max_abs'],max_slope=limits['bob_max_slope'],
                   max_turns=limits['bob_max_turns'],turning_tolerance=limits.get('bob_turning_tolerance',0.0))
    centrifugal = centrifugal_shape(p,model.masses,max(o.j for o in model.observations),grid)
    if not pec['ok'] or not br['ok'] or not centrifugal['ok']:
        raise FitError('PEC, BOB or centrifugal shape guard failed')
    return {**pec,'bob':br,'centrifugal':centrifugal}


def stats(residual,observations):
    w = np.array([o.weight for o in observations]);active=w>0
    return {'count':len(w),'positive_weight_count':int(np.sum(active)),
            'rms_cm':float(np.sqrt(np.mean(residual[active]**2))),
            'weighted_rms_cm':float(np.sqrt(np.sum(w*residual**2)/sum(w))),
            'max_abs_cm':float(np.max(np.abs(residual[active])))}


def optimize(model,config,root,max_iterations=100,preserve_c6=False,robust=0.0,fit_below=None):
    """Damped Gauss-Newton steps; every accepted iterate passes physical guards."""
    fitting_model=replace(model,observations=[o for o in model.observations if fit_below is None or o.energy<=fit_below])
    if not fitting_model.observations:raise FitError('No levels below the requested calibration cutoff')
    spec = Spectrum(fitting_model,model.grid)
    p0 = dict(model.parameters);keys=list(model.fitted)
    if preserve_c6 and kind(p0)!='MLJ':raise FitError('C6 preservation requires MLJ')
    scales = np.array([1000. if k in ('DE','AE') else .01 if k.startswith('BR_') else .1 if k=='RE' else 1. for k in keys])
    x = np.array([p0[k] for k in keys])/scales
    observed = np.array([o.energy for o in spec.observations])
    w = np.sqrt(np.array([o.weight for o in spec.observations]))
    w /= np.sqrt(np.mean(w*w))
    # Duo refinement.f90: sigma=robust/sqrt(normalized input weight), alpha=.001.
    base_weights=np.array([o.weight for o in spec.observations]);base_weights/=sum(base_weights)
    sigma2=robust**2/np.maximum(base_weights,1e-300)
    def score(residual):
        if robust>0:return float(np.sum(np.log1p(.001*residual**2/sigma2)))
        return float(np.sum((residual*w)**2))
    def unpack(xx):
        p={**p0,**{k:float(v) for k,v in zip(keys,xx*scales)}}
        if preserve_c6:
            p['PHIINF']=p0['PHIINF']+.5*math.log(p['AE']*p['RE']**p['N']/(p0['AE']*p0['RE']**p0['N']))
        return p
    physical(model,p0,config)
    calculated,info,_,_ = spec.evaluate(p0)
    objective = score(calculated-observed)
    history=[];damping=1e-8;start=time.monotonic()
    for iteration in range(max_iterations):
        p=unpack(x);derivs=[[],[]]
        for k in range(len(x)):
            h=1e-5;dx=np.zeros(len(x));dx[k]=h
            plus=curves(spec.r,unpack(x+dx));minus=curves(spec.r,unpack(x-dx))
            for i in range(2):derivs[i].append((plus[i]-minus[i])/(2*h))
        calculated,info,jac,_=spec.evaluate(p,tuple(np.array(d).T for d in derivs))
        fit_w=w
        if robust>0:
            fit_w=1/np.sqrt(sigma2+.001*(calculated-observed)**2)
            fit_w/=np.sqrt(np.mean(fit_w*fit_w))
        residual=(calculated-observed)*fit_w;jac*=fit_w[:,None]
        norms=np.maximum(np.linalg.norm(jac,axis=0),1e-12)
        u,s,vt=np.linalg.svd(jac/norms,full_matrices=False)
        accepted=False;old=objective;reason='No improving physical step'
        for retry in range(7):
            step=-(vt.T@(s/(s*s+damping)*(u.T@residual)))/norms
            for backtrack in range(12):
                fraction=2.**(-backtrack);candidate=x+fraction*step
                try:
                    proposed=unpack(candidate);physical(model,proposed,config)
                    energies,new_info,_,_=spec.evaluate(proposed)
                    new_objective=score(energies-observed)
                    if new_objective < objective:
                        x=candidate;objective=new_objective;calculated=energies;info=new_info
                        accepted=True;damping=max(1e-14,damping*.3);break
                except (FitError,FloatingPointError,OverflowError,ValueError) as exc:
                    reason=str(exc)
            if accepted:break
            damping*=10
        current=stats(observed-calculated,spec.observations)
        row={'iteration':iteration+1,'accepted':accepted,'damping':damping,'objective':objective,
             'parameter_step_fraction':fraction if accepted else 0.,'last_rejection_reason':reason,**current}
        if not accepted:row['reason']=reason
        history.append(row)
        write_json(root/'progress.json',{'parameters':unpack(x),'history':history})
        print(f"Localized {iteration+1}: RMS {current['rms_cm']:.8g}; weighted {current['weighted_rms_cm']:.8g} cm-1; {time.monotonic()-start:.0f}s",flush=True)
        if not accepted or (old-objective)/max(old,1e-30)<1e-8:break
    return unpack(x),history


def precise_duo_energies(path):
    result={}
    for line in Path(path).read_text().splitlines():
        if '||' not in line:continue
        row=line.split('||')[0].split()
        if len(row)!=10:raise FitError('Malformed precise Duo energy record')
        j,rank=int(float(row[0])),int(row[1]);value=float(row[2])
        if (j,rank) in result or not math.isfinite(value):raise FitError('Duplicate/nonfinite precise Duo level')
        result[j,rank]=value
    if not result:raise FitError('No precise Duo rovibronic energies')
    return result


def optimize_penalized(model,config,root,max_iterations=150,preserve_c6=False,robust=0.,fit_below=None):
    """Trust-region least squares with shape penalties and strict final screening.

    Penalties guide trial steps away from boundaries. Only fully physical
    candidates can become the returned model; failed trial spectra are rejected.
    """
    from scipy.optimize import least_squares
    subset=replace(model,observations=[o for o in model.observations if fit_below is None or o.energy<=fit_below])
    if not subset.observations:raise FitError('No calibration levels')
    spec=Spectrum(subset,model.grid);keys=list(model.fitted);p0=dict(model.parameters)
    physical(model,p0,config)
    if preserve_c6 and kind(p0)!='MLJ':raise FitError('C6 preservation requires MLJ')
    scales=np.array([1000. if k in ('AE','DE') else .01 if k.startswith('BR_') else .1 if k=='RE' else 1. for k in keys])
    x0=np.array([p0[k] for k in keys])/scales
    observed=np.array([o.energy for o in spec.observations]);weights=np.array([o.weight for o in spec.observations])
    wn=weights/weights.sum();w=np.sqrt(weights/weights.mean())
    limits=config['validation'];rg=np.linspace(*(limits['shape_range_angstrom'] or model.grid[:2]),2001)
    vref=curves(rg,p0)[0];tailmask=rg>2.2*p0['RE'];kmask=rg[1:]>2*p0['RE']
    def unpack(x):
        p={**p0,**{k:float(v) for k,v in zip(keys,x*scales)}}
        if preserve_c6:p['PHIINF']=p0['PHIINF']+.5*math.log(p['AE']*p['RE']**p['N']/(p0['AE']*p0['RE']**p0['N']))
        return p
    def penalties(x):
        p=unpack(x);v,q=curves(rg,p);dv=np.gradient(v,rg);dq=np.gradient(q,rg)
        k=dv*rg**3/(2*(1+q)-rg*dq)
        ceiling=limits.get('potential_ceiling_cm')
        return np.r_[np.maximum(-dv*np.sign(rg-p['RE']),0)/100,
                     np.maximum(np.abs(q)-limits['bob_max_abs'],0)*1e5,
                     np.maximum(np.abs(dq)-limits['bob_max_slope'],0)*1e4,
                     np.maximum(np.diff(k)[kmask],0)/100,(v-vref)[tailmask]*1e-5,
                     np.maximum(v/ceiling-1,0)*1e4 if ceiling else np.zeros(len(rg)),
                     max(.001-sum(p[f'B{i}'] for i in range(int(p['NR'])+1)),0)*1e4 if kind(p)=='EMO' else 0.]
    count=len(observed)+len(penalties(x0));cache={};history=[];best={};start=time.monotonic()
    def compute(x):
        if 'x' in cache and np.array_equal(cache['x'],x):return cache['f'],cache['jac']
        try:
            p=unpack(x);dv=[];dq=[];dp=[];h=1e-5
            for i in range(len(x)):
                step=np.zeros(len(x));step[i]=h
                plus=curves(spec.r,unpack(x+step));minus=curves(spec.r,unpack(x-step))
                dv.append((plus[0]-minus[0])/(2*h));dq.append((plus[1]-minus[1])/(2*h))
                dp.append((penalties(x+step)-penalties(x-step))/(2*h))
            energies,_,jac,_=spec.evaluate(p,(np.array(dv).T,np.array(dq).T));residual=energies-observed
            penalty=penalties(x)
            if robust>0:
                a=.001*wn/robust**2;u=a*residual**2;logs=np.log1p(u)
                transformed=np.sign(residual)*np.sqrt(logs);factor=np.sqrt(a)
                nz=np.abs(residual)>1e-12
                factor[nz]=a[nz]*np.abs(residual[nz])/((1+u[nz])*np.sqrt(logs[nz]))
            else:transformed=residual*w;factor=w
            f=np.r_[transformed,penalty];jj=np.vstack([jac*factor[:,None],np.array(dp).T])
            if not np.all(np.isfinite(f)) or not np.all(np.isfinite(jj)):raise FitError('Nonfinite trial')
            objective=float(f@f);current=stats(residual,spec.observations)
            acceptable=False
            try:physical(model,p,config);acceptable=True
            except (FitError,OverflowError,ValueError):pass
            if acceptable and (not best or objective<best['objective']):
                best.update(parameters=p,objective=objective,**current)
                write_json(root/'best_candidate.json',best)
            history.append({'evaluation':len(history)+1,'physical':acceptable,'objective':objective,**current})
            write_json(root/'progress.json',{'parameters':best.get('parameters',p0),'history':history})
            if len(history)%5==1:print(f"Penalized {len(history)}: RMS {current['rms_cm']:.8g}; weighted {current['weighted_rms_cm']:.8g}; physical={acceptable}; {time.monotonic()-start:.0f}s",flush=True)
        except (FitError,FloatingPointError,OverflowError,ValueError):
            if not best:raise
            f=np.full(count,1e10);jj=np.zeros((count,len(x)))
        cache.update(x=x.copy(),f=f,jac=jj)
        return f,jj
    lo=[];hi=[]
    for k in keys:
        if k in ('AE','DE'):lo.append(1.);hi.append(max(1e6,2*p0[k]))
        elif k=='RE':lo.append(model.grid[0]+.01);hi.append(model.grid[1]-.01)
        else:
            bound=max(100. if k.startswith('BR_') else 10000.,abs(p0[k])*2+1)
            lo.append(-bound);hi.append(bound)
    fit=least_squares(lambda x:compute(x)[0],x0,jac=lambda x:compute(x)[1],x_scale='jac',
                      bounds=(np.array(lo)/scales,np.array(hi)/scales),max_nfev=max_iterations,
                      ftol=1e-10,xtol=1e-10,gtol=1e-9)
    compute(fit.x)
    if not best:raise FitError('No physical candidate retained')
    write_json(root/'optimizer.json',{'message':str(fit.message),'nfev':fit.nfev,'best_objective':best['objective']})
    return best['parameters'],history


def native_check(model,p,grid,config,runner,label):
    # A mixed resonance can share its inner probability across two box levels.
    # Keep this permissive diagnostic separate from numerical convergence.
    spec=Spectrum(model,grid,minimum_probability=.25)
    energies,info,_,zpe=spec.evaluate(p)
    mapped=[replace(o,v=m['rank']-1,n=m['rank']) for o,m in zip(spec.observations,info)]
    mapped_model=replace(model,observations=mapped,parameters=p,grid=tuple(grid))
    vmax=max(config['vmax'],max(m['rank'] for m in info)+10)
    for attempt in range(5):
        text='PRINT_ROVIBRONIC_ENERGIES_TO_FILE\n'+render(mapped_model,p,jmax=max(o.j for o in mapped),vmax=vmax,iterations=0)
        run=runner.run(label+(f'_v{vmax}' if attempt else ''),text)
        assess_energies(parse_en(run.directory/'fit.en'),mapped,max(o.j for o in mapped))
        native=precise_duo_energies(run.directory/'rovibronic_energies.dat')
        rows=[]
        for o,m,computed in zip(spec.observations,info,energies):
            value=native[o.j,m['rank']]
            rows.append({**m,'observed_cm':o.energy,'calculated_cm':value,'residual_cm':o.energy-value,
                         'input_weight':o.weight,'python_duo_delta_cm':value-float(computed)})
        delta=max(abs(r['python_duo_delta_cm']) for r in rows)
        write_csv(run.directory/'state_map_residuals.csv',rows)
        write_json(run.directory/'primitive_contracted_comparison.json',{'vmax':vmax,'max_delta_cm':delta})
        if delta<=1e-5:return rows,run,zpe
        next_vmax=min(grid[2],max(vmax+40,int(vmax*1.5)))
        if next_vmax==vmax:break
        vmax=next_vmax
    raise FitError(f'Independent primitive solver disagrees with converged Duo contraction by {delta:.6g} cm-1')


def execute(config_path,executable,*,output=None,max_iterations=100,evaluate_only=False,preserve_c6=False,robust=None,fit_below=None,optimizer='guarded'):
    path,config,model,sources=load(config_path)
    robust=config['experiment']['robust'] if robust is None else robust
    if not math.isfinite(robust) or robust<0:raise FitError('Robust must be finite and nonnegative')
    if fit_below is not None and (not math.isfinite(fit_below) or fit_below<=0):raise FitError('Calibration energy cutoff must be positive')
    if optimizer not in ('guarded','penalized'):raise FitError('Unknown localized optimizer')
    if not config['refine_only']:
        raise FitError('Localized refinement requires refine_only=true and a prepared potential seed')
    root=Path(output).resolve() if output else (path.parent/config['output_directory']).resolve()
    if root.exists() and any(root.iterdir()):raise FitError('Localized output directory must be empty')
    root.mkdir(parents=True,exist_ok=True)
    runner=Runner(executable,root/'attempts',threads=config['threads'],timeout=config['timeout_seconds'])
    write_json(root/'manifest.json',{'config':config,'sources':sources,'solver':'localized_sinc_dvr',
               'code':{p.name:digest(p) for p in Path(__file__).parent.glob('*.py')},
               'executable':str(runner.executable),'executable_sha256':runner.executable_hash,
               'max_iterations':max_iterations,'evaluate_only':evaluate_only,'preserve_c6':preserve_c6,'robust':robust,'fit_below_cm':fit_below,'optimizer':optimizer,
               'numpy':np.__version__,'aston':ASTON,'minimum_inner_probability_fit':.5,'minimum_inner_probability_diagnostic':.25})
    fitter=optimize_penalized if optimizer=='penalized' else optimize
    p,history=(model.parameters,[]) if evaluate_only else fitter(model,config,root,max_iterations,preserve_c6,robust,fit_below)
    p={k:float(v) for k,v in p.items()}
    shape=physical(model,p,config)
    rows,run,zpe=native_check(model,p,model.grid,config,runner,'final')
    tol=config['validation']['energy_tolerance_cm'];lo,hi,n=model.grid;v=config['validation']
    scans=[];differences=[]
    for i in (1,2):
        grid=(lo-i*v['range_extension_angstrom'][0],hi+i*v['range_extension_angstrom'][1],n+i*v['npoints_increment'])
        if grid[0]<=0:raise FitError('Expanded radial grid extends below zero')
        expanded,_,_=native_check(model,p,grid,config,runner,f'grid_check_{i}')
        delta=[b['calculated_cm']-a['calculated_cm'] for a,b in zip(rows,expanded)]
        scans.append({'grid':grid,'max_delta_cm':max(map(abs,delta)),
                      'max_bound_delta_cm':max((abs(d) for r,d in zip(rows,delta) if not r['quasibound']),default=0.)})
        differences.extend({'J':a['J'],'v':a['v'],'quasibound':a['quasibound'],'grid_index':i,
                            'original_rank':a['rank'],'expanded_rank':b['rank'],'difference_cm':d}
                           for a,b,d in zip(rows,expanded,delta))
    residual=np.array([r['residual_cm'] for r in rows]);observations=sorted(model.observations,key=lambda o:o.key)
    total=stats(residual,observations);bound=[i for i,r in enumerate(rows) if not r['quasibound']]
    bound_stats=stats(residual[bound],[observations[i] for i in bound]) if bound else None
    calibration=[i for i,o in enumerate(observations) if fit_below is None or o.energy<=fit_below]
    calibration_stats=stats(residual[calibration],[observations[i] for i in calibration])
    basis_delta=max(s['max_delta_cm'] for s in scans);basis_ok=basis_delta<=tol
    quality=total[config['experiment']['target_metric']]<=config['experiment']['target_rms_cm']
    result={'status':'success' if quality and basis_ok else 'needs_review','solver':'localized_sinc_dvr','robust':robust,
            'experiment':total,'bound_experiment':bound_stats,'quasibound_count':len(rows)-len(bound),
            'fit_below_cm':fit_below,'calibration':calibration_stats,'calibration_count':len(calibration),
            'shape':shape,'basis_converged':basis_ok,'basis_max_delta_cm':basis_delta,
            'bound_basis_converged':all(s['max_bound_delta_cm']<=tol for s in scans),'grid_checks':scans,
            'basis_tolerance_cm':tol,'quality_pass':quality,'target_metric':config['experiment']['target_metric'],
            'target_cm':config['experiment']['target_rms_cm'],'grid':model.grid,'parameters':p,
            'Jmax':max(o.j for o in observations),'zero_point_energy_cm':zpe,
            'native_max_delta_cm':max(abs(r['python_duo_delta_cm']) for r in rows),
            'minimum_inner_probability':min(r['inner_probability'] for r in rows),'fit_history':history}
    write_json(root/'result.json',result);write_json(root/'final_parameters.json',p)
    write_csv(root/'final_residuals.csv',rows);write_csv(root/'basis_convergence.csv',differences)
    (root/'final.inp').write_bytes((run.directory/'run.inp').read_bytes())
    rr=np.linspace(*(config['validation']['shape_range_angstrom'] or model.grid[:2]),2001);vv,qq=curves(rr,p)
    write_csv(root/'potential.csv',[{'r_angstrom':float(r),'potential_cm':float(value)} for r,value in zip(rr,vv)])
    if bob_keys(p):write_csv(root/'bob_rot.csv',[{'r_angstrom':float(r),'bob_dimensionless':float(value)} for r,value in zip(rr,qq)])
    (root/'report.md').write_text(f"# Localized-state fit: {result['status']}\n\n"
        f"All {len(rows)} input levels retained. RMS {total['rms_cm']:.8g} cm-1; fixed-weight RMS {total['weighted_rms_cm']:.8g} cm-1.\n\n"
        f"Calibration: {len(calibration)} levels; energy cutoff {fit_below} cm-1 (None means all levels). Calibration RMS {calibration_stats['rms_cm']:.8g} cm-1.\n\n"
        f"{result['quasibound_count']} fitted levels lie above dissociation. Maximum radial-box change: {basis_delta:.8g} cm-1. "
        f"Bound-state convergence: {result['bound_basis_converged']}. Full-data convergence: {basis_ok}.\n\n"
        "The internal v label is an inner-well node count. `final_residuals.csv` preserves the original physical v alongside the current energy rank. "
        "The Duo `final.inp` explicitly uses those ranks; it must be regenerated after changing the potential or grid. "
        "All final energies were checked against Duo's 12-decimal rovibronic output. "
        "Quasibound values are finite-box stabilization estimates, not complex resonance energies; failed box convergence is never marked successful.\n",encoding='utf-8')
    return result
