"""Run/compare PEC models on identical observations and plot their diagnostics.

Example: python compare_pec_models.py HCl/models/*.json --output HCl/pecfit_runs/comparison --plots
Use --run --resume to execute configurations and reuse verified Duo attempts.
"""
import argparse
import csv
import glob
import json
from pathlib import Path
import shutil

from pecfit.duo import write_csv, write_json
from pecfit.forms import asymptote, bob, kind, potential, centrifugal_shape
from pecfit.model import FitError,shape_check
from pecfit.pipeline import Pipeline, load, merge_config


def compare(configs, output, *, run=False, resume=False, executable=None, plots=False,calibration_max_energy=None,potential_ceiling_cm=None):
    if not configs:raise FitError('No configurations to compare')
    output=Path(output).resolve()
    loaded=[load(path) for path in configs]
    reference=[(o.j,o.v,o.state,o.energy,o.weight) for o in loaded[0][2].observations]
    calibration_keys={(j,v) for j,v,_,energy,weight in reference if calibration_max_energy is not None and energy<=calibration_max_energy and weight>0}
    if any([(o.j,o.v,o.state,o.energy,o.weight) for o in m.observations]!=reference for _,_,m,_ in loaded[1:]):
        raise FitError("Model comparisons require identical energies, assignments and original weights")
    output.mkdir(parents=True,exist_ok=True)
    rows=[]; results={}; roots={}
    for path,c,model,sources in loaded:
        name=path.stem
        if name in roots: raise FitError("Configuration stems must be unique")
        root=(path.parent/c["output_directory"]).resolve()
        roots[name]=root
        error=None
        if run:
            try:
                Pipeline(path,executable=executable,resume=resume).execute()
            except (FitError,OSError) as exc:
                error=str(exc)
        result_path=root/"result.json"
        if not result_path.exists() and (root/"audited_result.json").exists():
            result_path=root/"audited_result.json"
        if error or not result_path.exists():
            rows.append({"model":name,"status":"failed_or_incomplete","rms_cm":None,"weighted_rms_cm":None,
                         "max_abs_cm":None,"basis_delta_cm":None,"shape_ok":False,"eligible":False,"reason":error or "No completed result"})
            continue
        result=json.loads(result_path.read_text())
        provenance=json.loads((root/"manifest.json").read_text())
        if provenance["sources"]!=sources or merge_config(provenance["config"])!=c:
            raise FitError(f"Inputs changed after the run in {root}; use a fresh run")
        results[name]=result
        shape_ok=result["shape"]["ok"] and result["shape"].get("bob",{}).get("ok",True)
        audit=shape_check(result['parameters'],c['validation']['shape_range_angstrom'] or result['grid'],
                          max_potential=potential_ceiling_cm or c['validation'].get('potential_ceiling_cm'))
        shape_ok=shape_ok and audit['ok']
        write_json(output/f'{name}_potential_shape.json',audit)
        if model.masses:
            centrifugal=centrifugal_shape(result["parameters"],model.masses,result["Jmax"],c["validation"]["shape_range_angstrom"] or result["grid"])
            shape_ok=shape_ok and centrifugal["ok"]
            write_json(output/f"{name}_centrifugal.json",centrifugal)
        if result["experiment"]["count"]!=len(reference):
            raise FitError(f"Incomplete observation coverage in {root}")
        rows.append({"model":name,"status":result["status"],"rms_cm":result["experiment"]["rms_cm"],
                     "weighted_rms_cm":result["experiment"]["weighted_rms_cm"],"max_abs_cm":result["experiment"]["max_abs_cm"],
                     "basis_delta_cm":result["basis_max_delta_cm"],"shape_ok":shape_ok,
                     "eligible":shape_ok and result["basis_converged"],"reason":"" if shape_ok else audit.get('reason','Shape check failed')})
        if calibration_keys:
            with (root/result.get('residual_file','final_residuals.csv')).open() as stream:
                residuals=[r for r in csv.DictReader(stream) if (int(r['J']),int(r['v'])) in calibration_keys]
            if {(int(r['J']),int(r['v'])) for r in residuals}!=calibration_keys:raise FitError('Incomplete common calibration coverage')
            sq=sum(float(r['residual_cm'])**2 for r in residuals)
            wsq=sum(float(r['input_weight'])*float(r['residual_cm'])**2 for r in residuals)
            basis=[]
            if (root/'basis_convergence.csv').exists():
                with (root/'basis_convergence.csv').open() as stream:
                    basis=[r for r in csv.DictReader(stream) if (int(r['J']),int(r['v'])) in calibration_keys]
            covered={(int(r['J']),int(r['v'])) for r in basis}==calibration_keys
            delta=max((abs(float(r['difference_cm'])) for r in basis),default=float('inf'))
            rows[-1].update(calibration_count=len(residuals),calibration_rms_cm=(sq/len(residuals))**.5,
                            calibration_weighted_rms_cm=(wsq/sum(float(r['input_weight']) for r in residuals))**.5,
                            calibration_basis_delta_cm=delta if covered else None,
                            calibration_eligible=shape_ok and covered and delta<=result['basis_tolerance_cm'])
    for row in rows:
        for k in ('calibration_count','calibration_rms_cm','calibration_weighted_rms_cm','calibration_basis_delta_cm'):
            row.setdefault(k,None)
        row.setdefault('calibration_eligible',False)
    rows.sort(key=lambda r:(not r["eligible"],r["rms_cm"] if r["rms_cm"] is not None else float("inf")))
    best=next((r["model"] for r in rows if r["eligible"]),None)
    best_calibration=min((r for r in rows if r['calibration_eligible']),key=lambda r:r['calibration_weighted_rms_cm'],default={}).get('model')
    write_csv(output/"comparison.csv",rows)
    write_json(output/"comparison.json",{"selection":"lowest unweighted RMS among shape/basis-validated candidates",
              "best":best,"best_calibration":best_calibration,"calibration_max_energy_cm":calibration_max_energy,
              "potential_ceiling_cm":potential_ceiling_cm,"levels":len(reference),"configs":[str(p) for p,_,_,_ in loaded],"models":rows})
    for filename in ("final.inp","final_parameters.json","final_residuals.csv","potential.csv","basis_convergence.csv","bob_rot.csv"):
            if best and (roots[best]/filename).exists():
                shutil.copy2(roots[best]/filename,output/("best_"+filename))
            else:
                (output/("best_"+filename)).unlink(missing_ok=True)
    table="\n".join(f"| {r['model']} | {r['rms_cm']} | {r['weighted_rms_cm']} | {r['basis_delta_cm']} | {r['eligible']} |" for r in rows)
    calibration_table='\n'.join(f"| {r['model']} | {r['calibration_rms_cm']} | {r['calibration_weighted_rms_cm']} | {r['calibration_basis_delta_cm']} | {r['calibration_eligible']} |" for r in rows)
    (output/"comparison.md").write_text(f"""# PEC model comparison

All models are evaluated on the same {len(reference)} term energies and original weights.
Best shape/basis-validated candidate: **{best or 'none'}**. Selection uses full-data
unweighted RMS; a `needs_review` model still misses its configured accuracy target.

| Model | RMS / cm-1 | Weighted RMS / cm-1 | Basis change / cm-1 | Shape and basis passed |
|---|---:|---:|---:|---|
{table}

Common calibration comparison: {len(calibration_keys)} levels with E <= {calibration_max_energy} cm-1.
Best validated model for this subset: **{best_calibration or 'none'}**, selected
using the original weighted RMS. This does not certify the remaining levels.
Some trials were fitted to all levels and others only to this explicit subset;
the individual manifests record that distinction.

| Model | Subset RMS / cm-1 | Subset weighted RMS / cm-1 | Subset grid change / cm-1 | Subset validated |
|---|---:|---:|---:|---|
{calibration_table}

This comparison measures agreement with the supplied energies. It does not
estimate parameter uncertainties or accuracy on independent data. The fitted
BOB function is an effective correction for one isotopologue, not a
mass-dependent BOB model transferable to other isotopologues.
""",encoding="utf-8")
    if plots and results:
        plot_models(output,results,roots)
    print(f"Best validated candidate: {best}; report: {output/'comparison.md'}")
    return {"best":best,"best_calibration":best_calibration,"rows":rows}


def plot_models(output,results,roots):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        raise FitError("Plots require matplotlib and numpy (pip install matplotlib)") from exc
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout="constrained")
    radii=np.linspace(.75,5.,1501)
    for name,result in results.items():
        p=result["parameters"]
        with (roots[name]/result.get("residual_file","final_residuals.csv")).open() as stream:
            rows=list(csv.DictReader(stream))
        js=np.array([int(r["J"]) for r in rows]); vs=np.array([int(r["v"]) for r in rows])
        residual=np.array([float(r["residual_cm"]) for r in rows])
        axes[0,0].scatter(vs,residual,s=6,alpha=.35,label=name)
        axes[0,1].scatter(js,residual,s=6,alpha=.35,label=name)
        axes[1,0].plot(radii,[potential(float(r),p) for r in radii],label=name)
        axes[1,1].plot(radii,[bob(float(r),p) for r in radii],label=name)
    axes[0,0].set(xlabel="v",ylabel="E(states) - E(Duo) / cm$^{-1}$",title="Residuals by vibration")
    axes[0,1].set(xlabel="J",ylabel="E(states) - E(Duo) / cm$^{-1}$",title="Residuals by rotation")
    axes[1,0].set(xlabel="r / Angstrom",ylabel="V / cm$^{-1}$",ylim=(0,60000),title="Potential energy curves")
    axes[1,1].set(xlabel="r / Angstrom",ylabel="Dimensionless BOB correction",title="Rotational BOB")
    for ax in axes.flat:
        ax.grid(alpha=.2); ax.legend(fontsize=7)
    fig.savefig(output/"model_comparison.png",dpi=180)
    fig.savefig(output/"model_comparison.pdf")
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("configs",nargs="+")
    parser.add_argument("--output",required=True)
    parser.add_argument("--run",action="store_true")
    parser.add_argument("--resume",action="store_true")
    parser.add_argument("--duo")
    parser.add_argument("--plots",action="store_true")
    parser.add_argument('--calibration-max-energy',type=float)
    parser.add_argument('--potential-ceiling',type=float)
    args=parser.parse_args()
    paths=sorted({str(Path(p).resolve()) for pattern in args.configs for p in glob.glob(pattern)})
    # Metadata files in model directories are not pipeline configurations.
    paths=[p for p in paths if isinstance(data:=json.loads(Path(p).read_text()),dict) and "template" in data]
    if not paths: parser.error("No model configurations found")
    try:
        compare(paths,args.output,run=args.run,resume=args.resume,executable=args.duo,plots=args.plots,
                calibration_max_energy=args.calibration_max_energy,potential_ceiling_cm=args.potential_ceiling)
    except (FitError,OSError,ValueError) as exc:
        parser.exit(1,f"compare_pec_models: {exc}\n")


if __name__=="__main__": main()
