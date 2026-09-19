"""Publish fixed-data diagnostics and scientific plots for a selected native run."""
from pathlib import Path
import argparse
import csv
import json
import os
import sys
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'): os.environ[k]='1'
ROOT=Path(__file__).resolve().parent
deps=ROOT.parent/'.pecfit_stage/plot_deps'
if deps.exists(): sys.path.insert(0,str(deps))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import coupled
from barrier import potential, components, bob, audit

def key(r): return int(r['state']),int(r['v']),float(r['J']),r['parity']

def raw_rows(directory, original):
    with (directory/'residuals.csv').open() as f: unique={key(r):r for r in csv.DictReader(f)}
    result=[]
    for source_index,o in enumerate(coupled.observations(original),1):
        c=unique[key(o)]
        row=dict(source_index=source_index,**o,calculated_cm=float(c['calculated_cm']))
        row['residual_cm']=o['observed_cm']-row['calculated_cm']
        row['ef']='e' if (o['parity']=='+')==(int(o['J'])%2==0) else 'f'
        result.append(row)
    return result

def stats(rows):
    active=[r for r in rows if r['input_weight']>0]
    x=np.array([r['residual_cm'] for r in active]); w=np.array([r['input_weight'] for r in active])
    return dict(count=len(rows),positive_weight_count=len(active),
                rms_cm=float(np.sqrt(np.mean(x*x))),weighted_rms_cm=float(np.sqrt(w@(x*x)/sum(w))),
                median_abs_cm=float(np.median(abs(x))),p90_abs_cm=float(np.quantile(abs(x),.9)),
                p95_abs_cm=float(np.quantile(abs(x),.95)),max_abs_cm=float(max(abs(x))),
                above_1_cm=int(sum(abs(x)>1)),above_2_cm=int(sum(abs(x)>2)))

def describe(rows):
    return dict(all=stats(rows),by_state={str(s):stats([r for r in rows if r['state']==s]) for s in (1,2)})

def write_csv(path, rows):
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('selected',type=Path); ap.add_argument('--output',type=Path,default=ROOT/'report')
    ap.add_argument('--workdir',type=Path,default=ROOT)
    ap.add_argument('--ordinary',type=Path,default=ROOT/'runs/least_squares_rotation/round_05_fraction_1')
    a=ap.parse_args(); out=a.output; out.mkdir(parents=True,exist_ok=True)
    original=(a.workdir/'prepared.inp').read_text()
    comparisons={name:raw_rows(p,original) for name,p in dict(
        repaired_original=a.workdir/'runs/unique_baseline',barrier_seed=a.workdir/'runs/barrier_seed',
        ordinary_fit=a.ordinary,selected=a.selected).items()}
    report={name:describe(rows) for name,rows in comparisons.items()}
    report['selected_calculation']=str(a.selected.resolve())
    report['shape']=audit((a.selected/'input.inp').read_text())
    (out/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    rows=comparisons['selected']; write_csv(out/'all_original_observation_residuals.csv',rows)
    flagged=[r for r in rows if r['input_weight']>0 and abs(r['residual_cm'])>1]
    if flagged: write_csv(out/'levels_needing_review.csv',sorted(flagged,key=lambda r:-abs(r['residual_cm'])))
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    panels=[(1,'X residuals',None),(2,'A residuals: every positive-weight observation',None),
            (2,'A residuals: expanded vertical scale',(-.7,.7)),(2,'A residuals by vibrational level',(-.7,.7))]
    for ax,(s,title,ylim) in zip(axs.flat,panels):
        group=[r for r in rows if r['state']==s and r['input_weight']>0]
        for ef,marker in [('e','o'),('f','x')]:
            g=[r for r in group if r['ef']==ef]
            if not g: continue
            x=[r['v']+.1*(r['J']/100-.5) if ax==axs[1,1] else r['J'] for r in g]
            sc=ax.scatter(x,[r['residual_cm'] for r in g],c=[r['v'] for r in g],cmap='viridis',vmin=0,vmax=9 if s==1 else 6,s=12,alpha=.65,marker=marker,label=ef)
        ax.axhline(0,color='.4',lw=.8); ax.set_title(title); ax.set_ylabel('Observed − calculated / cm⁻¹'); ax.set_xlabel('v (small J offset)' if ax==axs[1,1] else 'J')
        if ylim:
            ax.set_ylim(*ylim)
            n=sum(not ylim[0]<=r['residual_cm']<=ylim[1] for r in group)
            ax.text(.02,.97,f'{n} observations outside this view',transform=ax.transAxes,va='top',fontsize=9)
        ax.legend(title='Parity',loc='lower left'); fig.colorbar(sc,ax=ax,label='v')
    fig.suptitle('AlF: retained source observations, including repeated assignments')
    fig.savefig(out/'residuals.png',dpi=180); plt.close(fig)
    fs=coupled.fields((a.selected/'input.inp').read_text()); old=coupled.fields(original)
    r=np.linspace(.9,10,4000); x=potential(r,fs['POTEN',1,1]); upper=potential(r,fs['POTEN',2,2]); c1,c2,c3=components(r,fs['POTEN',2,2].values)
    fig,axs=plt.subplots(1,3,figsize=(15,4.6),layout='constrained')
    axs[0].plot(r,x,label='X¹Σ⁺'); axs[0].plot(r,upper,label='A¹Π'); axs[0].axhline(55564.02471,ls=':',color='k',label='Dissociation limit'); axs[0].set(xlim=(1.15,8),ylim=(-1000,72000),title='Two physical PECs',ylabel='Energy above X minimum / cm⁻¹')
    axs[1].plot(r,upper,label='A: lower coupled branch',lw=2); axs[1].plot(r,c1,'--',label='Diabatic EMO'); axs[1].plot(r,c2,'--',label='repulsive_exp'); axs[1].plot(r,potential(r,old['POTEN',2,2]),':',label='Starting A EMO'); axs[1].set(xlim=(1.35,4),ylim=(43000,72000),title='Barrier construction')
    axs[2].plot(r,upper-55564.02471,label='A − dissociation limit'); axs[2].plot(r,-73495.9702828864/r**5-362561.9198027829/r**6,ls='--',label='−C₅/r⁵ − C₆/r⁶'); axs[2].axhline(0,color='.4',lw=.8); axs[2].set(xlim=(3.7,10),ylim=(-45,35),title='Long-range approach and shallow outer well')
    for ax in axs: ax.set_xlabel('r / Å'); ax.legend(fontsize=8)
    fig.savefig(out/'potential_curves.png',dpi=180); plt.close(fig)
    rr=np.linspace(.7,8,1600); f=fs['BOBROT',1,1]; vals=bob(rr,f.values)
    fig,ax=plt.subplots(figsize=(8,3.5),layout='constrained'); ax.plot(rr,vals); ax.axhline(0,color='.4',lw=.8); ax.set(xlabel='r / Å',ylabel='Dimensionless BOB-rot',title='AlF X rotational correction')
    fig.savefig(out/'bob_rot.png',dpi=180); plt.close(fig)
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
