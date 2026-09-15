"""Single-state PEC fitting; core Duo workflow uses the Python standard library."""
import argparse
import json
import os
from pathlib import Path
import sys

from .duo import parse_output, write_json
from .model import FitError, Model, observations, read_pec, render
from .pipeline import DEFAULTS, Pipeline, load
from .states import read_states


def initialize(args):
    dest = Path(args.config).resolve()
    if dest.exists():
        raise FitError(f"Config already exists: {dest}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    if args.template:
        template = Path(args.template).resolve()
        from .model import import_template
        import_template(template)
    else:
        if not args.pec or not (args.energies or args.states) or not args.atoms:
            raise FitError("Use --template, or provide --pec, --atoms and --energies/--states")
        pec = read_pec(Path(args.pec), args.units, True)
        obs = (read_states(args.states, j_column=args.j_column, v_column=args.v_column,
                           weight_scale=args.weight_scale) if args.states else
               observations(Path(args.energies).read_text().splitlines()))
        state = obs[0].state
        re = min(pec, key=lambda x: x[1])[0]
        p = {"V0": 0.0, "RE": re, "DE": args.de or max(e for r,e in pec if r >= re)*1.2,
             "RREF": -1.0, "PL": 6., "PR": 6., "NL": 4., "NR": 8.,
             **{f"B{k}": 2.0 if k == 0 else 0.0 for k in range(9)}}
        grid = tuple(args.grid) if args.grid else (pec[0][0], pec[-1][0], 301)
        if grid[2] != int(grid[2]):
            raise FitError("NPOINTS must be an integer")
        grid = (grid[0], grid[1], int(grid[2]))
        model = Model(" ".join(args.atoms), " ".join(map(str,args.masses)) if args.masses else None,
                      state, "X1Sigma+", grid, p, ["RE", "DE"]+[f"B{k}" for k in range(9)],
                      pec, "PS1997 1e-5 80000.0", obs)
        model.validate()
        template = dest.with_suffix(".template.inp")
        if template.exists():
            raise FitError(f"Template already exists: {template}")
        template.write_text(render(model, p, jmax=max(o.j for o in obs), vmax=max(80,max(o.v for o in obs)+1), iterations=0), encoding="ascii")
    config = json.loads(json.dumps(DEFAULTS))
    config["template"] = os.path.relpath(template, dest.parent).replace("\\", "/")
    from .model import import_template
    imported = import_template(template)
    config["vmax"] = max(80, max(o.v for o in imported.observations)+20)
    if args.states:
        config.update(states_file=os.path.relpath(Path(args.states).resolve(), dest.parent).replace("\\", "/"),
                      states_j_column=args.j_column, states_v_column=args.v_column,
                      states_weight_scale=args.weight_scale)
    if args.pec:
        config["pec_file"] = os.path.relpath(Path(args.pec).resolve(), dest.parent).replace("\\", "/")
        config["pec_units"] = args.units
    if args.duo:
        config["duo_executable"] = os.path.relpath(Path(args.duo).resolve(), dest.parent).replace("\\", "/")
    write_json(dest, config)
    print(f"Created {dest}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    init = commands.add_parser("init", help="Create a reusable molecule configuration")
    init.add_argument("--config", required=True)
    init.add_argument("--template")
    init.add_argument("--pec")
    init.add_argument("--units", choices=["hartree", "Eh", "cm-1"], default="hartree")
    energies = init.add_mutually_exclusive_group()
    energies.add_argument("--energies", help="Ten-column Duo experimental energy table, no header")
    energies.add_argument("--states", help="ExoMol states file for one X1Sigma+ state (also .bz2)")
    init.add_argument("--j-column", type=int, default=4, help="One-based J column in states file")
    init.add_argument("--v-column", type=int, default=5, help="One-based v column in states file")
    init.add_argument("--weight-scale", type=float, default=100., help="States weight = scale/sqrt(v+1)/sqrt(J+1)")
    init.add_argument("--atoms", nargs=2)
    init.add_argument("--masses", nargs=2, type=float)
    init.add_argument("--grid", nargs=3, type=float, metavar=("RMIN", "RMAX", "NPOINTS"))
    init.add_argument("--de", type=float, help="Initial well depth in cm-1 (not a fixed experimental value)")
    init.add_argument("--duo")
    plan = commands.add_parser("plan", help="Validate inputs and show the stage schedule without running Duo")
    plan.add_argument("config")
    run = commands.add_parser("run", help="Fit and validate automatically")
    run.add_argument("config")
    run.add_argument("--duo", help="Duo executable; alternatively set DUO_EXE")
    run.add_argument("--output", help="New output directory")
    run.add_argument("--resume", action="store_true", help="Replay decisions using verified cached attempts")
    local = commands.add_parser('refine-localized',help='Optional NumPy/SciPy refinement with explicit quasibound diagnostics')
    local.add_argument('config')
    local.add_argument('--duo',required=True)
    local.add_argument('--output')
    local.add_argument('--max-iterations',type=int,default=100)
    local.add_argument('--evaluate-only',action='store_true')
    local.add_argument('--preserve-c6',action='store_true')
    local.add_argument('--robust',type=float)
    local.add_argument('--fit-below',type=float,help='Explicit calibration cutoff; final diagnostics retain all source levels')
    local.add_argument('--optimizer',choices=('guarded','penalized'),default='guarded')
    inspect = commands.add_parser("inspect-output", help="Extract the final full-precision parameters from a Duo output")
    inspect.add_argument("output")
    convert = commands.add_parser("convert-states", help="Convert ExoMol single-state energies to a Duo table")
    convert.add_argument("states")
    convert.add_argument("output")
    convert.add_argument("--j-column", type=int, default=4)
    convert.add_argument("--v-column", type=int, default=5)
    convert.add_argument("--weight-scale", type=float, default=100.)
    try:
        if args := parser.parse_args(argv):
            if args.command == "init":
                initialize(args)
            elif args.command == "convert-states":
                obs = read_states(args.states, j_column=args.j_column, v_column=args.v_column,
                                  weight_scale=args.weight_scale)
                dest = Path(args.output)
                with dest.open("x", encoding="ascii") as stream:
                    stream.write("\n".join(o.line() for o in obs)+"\n")
                print(f"Wrote {len(obs)} levels to {dest}")
            elif args.command == "plan":
                _, c, m, _ = load(args.config)
                print(json.dumps({"atoms": m.atoms, "state": m.state, "pec_points": len(m.pec),
                                  "experimental_levels": len(m.observations),
                                  "Jmax": max(o.j for o in m.observations), "vmax_observed": max(o.v for o in m.observations),
                                  "jmax_schedule": c["jmax_schedule"], "vibrational_functions": c["vmax"],
                                  "grid": m.grid, "target_metric": c["experiment"]["target_metric"],
                                  "target_cm": c["experiment"]["target_rms_cm"],
                                  "de_fixed_to_experiment": c["experimental_de_cm"] is not None}, indent=2))
            elif args.command == 'refine-localized':
                if args.max_iterations<1:raise FitError('Iterations must be positive')
                _,cfg,_,_=load(args.config)
                for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]=str(cfg['threads'])
                try:
                    from .localized import execute
                except ImportError as exc:
                    raise FitError('Install numerical extras with pip install ".[localized]"') from exc
                result=execute(args.config,args.duo,output=args.output,max_iterations=args.max_iterations,
                               evaluate_only=args.evaluate_only,preserve_c6=args.preserve_c6,
                               robust=args.robust,fit_below=args.fit_below,optimizer=args.optimizer)
                print(json.dumps({k:result[k] for k in ('status','experiment','calibration','basis_max_delta_cm')},indent=2))
                return 0 if result['status']=='success' else 2
            elif args.command == "inspect-output":
                iterations, history = parse_output(Path(args.output).read_text(errors="replace"))
                if not iterations:
                    raise FitError("No full-precision fitted parameters found")
                print(json.dumps({"last_iteration": iterations[-1], "last_summary": history[-1] if history else None}, indent=2))
            else:
                pipeline = Pipeline(args.config, executable=args.duo, output=args.output, resume=args.resume)
                try:
                    result = pipeline.execute()
                except (FitError, OSError, KeyboardInterrupt) as exc:
                    write_json(pipeline.root/"failure.json", {"status": "interrupted" if isinstance(exc,KeyboardInterrupt) else "failed", "reason": str(exc)})
                    raise
                (pipeline.root/"failure.json").unlink(missing_ok=True)
                print(f"{result['status']}: {pipeline.root / 'report.md'}")
                return 0 if result["status"] == "success" else 2
        return 0
    except KeyboardInterrupt:
        print("Interrupted; rerun with --resume.", file=sys.stderr)
        return 130
    except (FitError, OSError, ValueError, KeyError) as exc:
        print(f"pecfit: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
