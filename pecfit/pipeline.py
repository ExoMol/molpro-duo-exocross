"""Bounded fit/restart/validate workflow. All acceptance uses fresh calculations."""
from __future__ import annotations

from dataclasses import replace
import hashlib
import json
import math
import os
from pathlib import Path

from .duo import Runner, assess_energies, digest, parse_en, write_csv, write_json
from .model import FitError, emo, import_template, read_pec, render, shape_check
from .states import read_states
from .forms import potential, kind, asymptote, bob, bob_shape, bob_keys, centrifugal_shape


DEFAULTS = {
    "schema_version": 1,
    "template": "sample_outputs/CO_duo_fit_01.inp",
    "duo_executable": None,
    "output_directory": "pecfit_runs",
    "pec_file": None,
    "states_file": None,
    "states_j_column": 4,
    "states_v_column": 5,
    "states_weight_scale": 100.0,
    "refine_only": False,
    "pec_units": "hartree",
    "pec_shift_minimum": True,
    "experimental_de_cm": None,
    "jmax_schedule": [10, 40, "all"],
    "vmax": 80,
    "timeout_seconds": 1800,
    "threads": 1,
    "abinitio": {"iterations": 400, "max_restarts": 4, "fit_scale": 0.1,
                  "fit_factor": 1e7, "energy_fit_factor": 1e-10,
                  "equilibrium_window_angstrom": 0.15, "equilibrium_rms_limit_cm": 100,
                  "outer_rms_limit_cm": 2000, "inner_max_abs_limit_cm": 100000},
    "experiment": {"iterations": 200, "max_restarts": 8, "fit_scale": 0.1,
                   "bob_prefit": False,
                   "line_search_steps": 0,
                   "initial_fit_scale": 0.01, "robust": 0.0001, "polish_robust": 0.0,
                   "fit_factor": 1e6, "abinitio_fit_factor": 1e-7,
                   "outlier_threshold_cm": 10, "max_retries": 3,
                   "target_rms_cm": 0.1, "target_metric": "rms_cm",
                   "minimum_iterations": 20,
                   "stability_tolerance": 1e-7, "max_right_order": 10,
                   "order_step": 2, "minimum_relative_improvement": 1e-5},
    "validation": {"vmax_increment": 20, "npoints_increment": 100,
                   "range_extension_angstrom": [0.0, 0.0],
                   "max_basis_refinements": 2, "energy_tolerance_cm": 0.001,
                   "shape_tolerance_cm": 0.001,
                   "bob_max_abs": 0.1, "bob_max_slope": 0.5, "bob_max_turns": 4,"bob_turning_tolerance":0.0,
                   "shape_range_angstrom": None,"potential_ceiling_cm":None, "check_centrifugal_wells": False}
}


def merge_config(data):
    config = json.loads(json.dumps(DEFAULTS))
    unknown = set(data)-set(config)
    if unknown:
        raise FitError(f"Unknown configuration keys: {sorted(unknown)}")
    for k, v in data.items():
        if isinstance(config[k], dict):
            if not isinstance(v, dict) or set(v)-set(config[k]):
                raise FitError(f"Unknown or invalid {k} options")
            config[k].update(v)
        else:
            config[k] = v
    if config["schema_version"] != 1:
        raise FitError("Unsupported schema_version")
    for section, defaults in ((config, DEFAULTS), (config["abinitio"], DEFAULTS["abinitio"]),
                              (config["experiment"], DEFAULTS["experiment"]),
                              (config["validation"], DEFAULTS["validation"])):
        for key, default in defaults.items():
            if isinstance(default, (int, float)) and not isinstance(default, bool):
                value = section[key]
                if not isinstance(value, (int, float)) or isinstance(value, bool):
                    raise FitError(f"{key} must be numeric")
    if not isinstance(config["pec_shift_minimum"], bool):
        raise FitError("pec_shift_minimum must be true or false")
    if not isinstance(config["refine_only"],bool):
        raise FitError("refine_only must be true or false")
    if not isinstance(config["validation"]["check_centrifugal_wells"],bool):
        raise FitError("check_centrifugal_wells must be true or false")
    for key in ("template", "output_directory", "pec_units"):
        if not isinstance(config[key], str) or not config[key]:
            raise FitError(f"{key} must be a nonempty string")
    for key in ("pec_file", "states_file", "duo_executable"):
        if config[key] is not None and (not isinstance(config[key], str) or not config[key]):
            raise FitError(f"{key} must be null or a nonempty string")
    de = config["experimental_de_cm"]
    if de is not None and (not isinstance(de, (int, float)) or isinstance(de,bool) or not math.isfinite(de) or de <= 0):
        raise FitError("experimental_de_cm must be null or a positive well depth")
    for section in (config, config["abinitio"], config["experiment"], config["validation"]):
        for key, value in section.items():
            if isinstance(value, (int, float)) and (not math.isfinite(value) or value < 0):
                raise FitError(f"{key} must be finite and nonnegative")
    for section, keys in ((config, ("vmax", "threads", "timeout_seconds")),
                          (config["abinitio"], ("iterations", "max_restarts")),
                          (config["experiment"], ("iterations", "max_restarts", "order_step", "minimum_iterations")),
                          (config["validation"], ("vmax_increment", "npoints_increment"))):
        for key in keys:
            if not isinstance(section[key], int) or section[key] < 1:
                raise FitError(f"{key} must be a positive integer")
    e = config["experiment"]
    if not isinstance(e["bob_prefit"],bool):
        raise FitError("bob_prefit must be true or false")
    if (any(not isinstance(config[k], int) or config[k] < 4 for k in ("states_j_column", "states_v_column"))
            or config["states_j_column"] == config["states_v_column"] or config["states_weight_scale"] <= 0):
        raise FitError("States columns must be distinct integers >=4 and weight scale positive")
    if e["target_metric"] not in ("rms_cm", "weighted_rms_cm"):
        raise FitError("target_metric must be rms_cm or weighted_rms_cm")
    if not 0 < e["initial_fit_scale"] <= 1 or not 0 < e["fit_scale"] <= 1 or not 0 < config["abinitio"]["fit_scale"] <= 1:
        raise FitError("Fit scales must be in (0,1]")
    if not isinstance(config["jmax_schedule"], list) or not config["jmax_schedule"]:
        raise FitError("jmax_schedule must be a nonempty list")
    for j in config["jmax_schedule"]:
        if j != "all" and (not isinstance(j, int) or j < 0):
            raise FitError("J schedule entries must be nonnegative integers or 'all'")
    ext = config["validation"]["range_extension_angstrom"]
    if not isinstance(ext, list) or len(ext) != 2 or any(not isinstance(x,(int,float)) or not math.isfinite(x) or x < 0 for x in ext):
        raise FitError("range_extension_angstrom requires two nonnegative numbers")
    shape_range=config["validation"]["shape_range_angstrom"]
    if shape_range is not None and (not isinstance(shape_range,list) or len(shape_range)!=2
            or any(not isinstance(x,(int,float)) or not math.isfinite(x) for x in shape_range)
            or not 0<shape_range[0]<shape_range[1]):
        raise FitError("shape_range_angstrom must be null or an increasing positive interval")
    if not isinstance(config["validation"]["bob_max_turns"],int):
        raise FitError("bob_max_turns must be an integer")
    ceiling=config['validation']['potential_ceiling_cm']
    if ceiling is not None and (not isinstance(ceiling,(int,float)) or isinstance(ceiling,bool) or not math.isfinite(ceiling) or ceiling<=0):
        raise FitError('potential_ceiling_cm must be null or finite and positive')
    for section, keys in ((e, ("max_retries", "max_right_order", "line_search_steps")),
                          (config["validation"], ("max_basis_refinements",))):
        for key in keys:
            if not isinstance(section[key], int):
                raise FitError(f"{key} must be an integer")
    return config


def load(path):
    path = Path(path).resolve()
    c = merge_config(json.loads(path.read_text(encoding="utf-8-sig")))
    template = (path.parent/c["template"]).resolve()
    m = import_template(template)
    sources = {str(template): digest(template)}
    if c["states_file"]:
        states_path = (path.parent/c["states_file"]).resolve()
        m.observations = read_states(states_path, j_column=c["states_j_column"], v_column=c["states_v_column"],
                                     weight_scale=c["states_weight_scale"], state=m.state)
        sources[str(states_path)] = digest(states_path)
    if c["pec_file"]:
        pec_path = (path.parent/c["pec_file"]).resolve()
        m.pec = read_pec(pec_path, c["pec_units"], c["pec_shift_minimum"])
        sources[str(pec_path)] = digest(pec_path)
    m.validate()
    if c["vmax"] <= max(o.v for o in m.observations):
        raise FitError("vmax must exceed the largest experimental v (Duo counts basis functions)")
    if c["vmax"] > m.grid[2]:
        raise FitError("vmax exceeds radial grid size")
    jall = max(o.j for o in m.observations)
    schedule = sorted({min(jall, j) if j != "all" else jall for j in c["jmax_schedule"]} | {jall})
    c["jmax_schedule"] = schedule
    return path, c, m, sources


def pec_diagnostics(model, p, settings):
    near, outer, inner, rows = [], [], [], []
    for r, observed in model.pec:
        calc = potential(r, p)
        if not math.isfinite(calc):
            raise FitError("Non-finite EMO at an ab initio point")
        residual = observed-calc
        if abs(r-p["RE"]) <= settings["equilibrium_window_angstrom"]:
            near.append(residual)
        elif r > p["RE"]:
            outer.append(residual)
        else:
            inner.append(residual)
        rows.append({"r_angstrom": r, "abinitio_cm": observed, "emo_cm": calc, "residual_cm": residual})
    def rms(seq):
        return math.sqrt(sum(x*x for x in seq)/len(seq)) if seq else None
    stats = {"equilibrium_rms_cm": rms(near), "outer_rms_cm": rms(outer),
             "inner_max_abs_cm": max(map(abs, inner), default=0)}
    stats["ok"] = (bool(near) and stats["equilibrium_rms_cm"] <= settings["equilibrium_rms_limit_cm"]
                   and (not outer or stats["outer_rms_cm"] <= settings["outer_rms_limit_cm"])
                   and stats["inner_max_abs_cm"] <= settings["inner_max_abs_limit_cm"])
    return stats, rows


class Pipeline:
    def __init__(self, config_path, *, executable=None, output=None, resume=False):
        self.path, self.config, self.model, sources = load(config_path)
        c = self.config
        exe = executable or os.environ.get("DUO_EXE") or c["duo_executable"]
        if not exe:
            import shutil
            exe = shutil.which("duo") or shutil.which("duo.exe")
        if not exe:
            raise FitError("Specify --duo /path/to/duo or set DUO_EXE (no binary is downloaded automatically)")
        if not executable and not os.environ.get("DUO_EXE") and c["duo_executable"] and ("/" in exe or "\\" in exe):
            exe = self.path.parent/exe
        self.root = Path(output).resolve() if output else (self.path.parent/c["output_directory"]).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.runner = Runner(exe, self.root/"attempts", timeout=c["timeout_seconds"], threads=c["threads"], resume=resume)
        from . import __version__
        code_hashes = {p.name: digest(p) for p in Path(__file__).parent.glob("*.py")}
        provenance = {"config": c, "sources": sources, "executable": str(self.runner.executable),
                      "executable_sha256": self.runner.executable_hash, "version": __version__, "code": code_hashes}
        signature = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
        manifest = self.root/"manifest.json"
        if manifest.exists():
            old = json.loads(manifest.read_text())
            if not resume:
                raise FitError("Output already contains a run; use --resume or a new --output directory")
            if old["signature"] != signature:
                raise FitError("Config, source, executable or pipeline code changed; use a new output directory")
        elif any(self.root.iterdir()):
            raise FitError("Output must be empty for a new run")
        write_json(manifest, {"signature": signature, **provenance})
        self.log = []
        self.p = dict(self.model.parameters)
        self.fitted = list(self.model.fitted)
        self.vmax = c["vmax"]
        self.grid = self.model.grid
        self.counter = 0

    def record(self, **data):
        self.log.append(data)
        write_json(self.root/"progress.json", {"events": self.log, "parameters": self.p,
                                               "vmax": self.vmax, "grid": self.grid})

    def run_duo(self, label, p, jmax, iterations, *, abinitio=False, scale=0.1, robust=0.0, threshold=None, grid=None, vmax=None):
        c = self.config
        a, e = c["abinitio"], c["experiment"]
        self.counter += 1
        text = render(self.model, p, jmax=jmax, vmax=vmax or (1 if abinitio else self.vmax),
                      iterations=iterations, scale=scale, robust=robust,
                      abi_factor=a["fit_factor"] if abinitio else e["abinitio_fit_factor"],
                      energy_factor=a["energy_fit_factor"] if abinitio else e["fit_factor"],
                      fitted=self.fitted, threshold=threshold, grid=grid or self.grid)
        return self.runner.run(f"{self.counter:04d}_{label}", text, fitting=iterations > 0,
                               stability_tolerance=e["stability_tolerance"] if iterations > 0 else None,
                               minimum_iterations=e["minimum_iterations"])

    def physical(self, p):
        replace(self.model, parameters=p, fitted=self.fitted, grid=self.grid).validate()
        shape_grid=self.config["validation"]["shape_range_angstrom"] or self.grid
        shape = shape_check(p, shape_grid, self.config["validation"]["shape_tolerance_cm"],max_potential=self.config['validation']['potential_ceiling_cm'])
        if not shape["ok"]:
            raise FitError(f"Rejected nonphysical PEC: {shape['reason']}")
        limits=self.config["validation"]
        br=bob_shape(p,shape_grid,max_abs=limits["bob_max_abs"],max_slope=limits["bob_max_slope"],max_turns=limits["bob_max_turns"],turning_tolerance=limits["bob_turning_tolerance"])
        if not br["ok"]:
            raise FitError(f"Rejected nonphysical BOB: {br['reason']}")
        shape["bob"]=br
        if limits["check_centrifugal_wells"]:
            centrifugal=centrifugal_shape(p,self.model.masses,max(o.j for o in self.model.observations),shape_grid)
            if not centrifugal["ok"]:
                raise FitError("Rejected PEC/BOB combination with extra centrifugal wells or barriers")
            shape["centrifugal"]=centrifugal
        return shape

    def evaluate(self, label, p, jmax, *, grid=None, vmax=None):
        self.physical(p)
        # No outlier rejection or robust reweighting in validation: all observations count.
        run = self.run_duo(label, p, jmax, 0, grid=grid, vmax=vmax)
        levels = parse_en(run.directory/"fit.en")
        stats, rows = assess_energies(levels, self.model.observations, jmax)
        write_csv(run.directory/"residuals.csv", rows)
        write_json(run.directory/"validation.json", stats)
        return stats, levels, run

    def fit_abinitio(self):
        a = self.config["abinitio"]
        scale = a["fit_scale"]
        best_score = math.inf
        for restart in range(a["max_restarts"]):
            try:
                run = self.run_duo(f"abinitio_{restart}", self.p, 0, a["iterations"], abinitio=True, scale=scale)
                p = run.parameters
                self.physical(p)
                # Re-evaluate parameters, since printed fit residuals can precede the update.
                check = self.run_duo(f"abinitio_check_{restart}", p, 0, 0, abinitio=True)
                stats, rows = pec_diagnostics(self.model, p, a)
                score = run.history[-1]["duo_weighted_rms"]
                if score > best_score*(1+1e-6):
                    raise FitError("Ab initio weighted objective increased")
                self.p, best_score = p, score
                write_csv(check.directory/"pec_residuals.csv", rows)
                self.record(stage="abinitio", accepted=True, diagnostics=stats, history=run.history[-1])
                print(f"  Ab initio: near-equilibrium RMS {stats['equilibrium_rms_cm']:.6g} cm-1", flush=True)
                if stats["ok"] and run.history[-1]["stability"] <= self.config["experiment"]["stability_tolerance"]:
                    write_json(self.root/"abinitio_parameters.json", self.p)
                    return
            except FitError as exc:
                self.record(stage="abinitio", accepted=False, reason=str(exc))
                scale *= 0.5
        raise FitError("Ab initio fit did not pass convergence and residual limits; see progress.json")

    def fit_experiment(self, jmax, label, initial=False, robust=None):
        e = self.config["experiment"]
        metric = e["target_metric"]
        robust = e["robust"] if robust is None else robust
        stats, _, _ = self.evaluate(label+"_before", self.p, jmax)
        scale = e["initial_fit_scale"] if initial else e["fit_scale"]
        for restart in range(e["max_restarts"]):
            accepted = False
            for retry in range(e["max_retries"]+1):
                # Activate only once every positive-weight level is already inside the cutoff.
                cutoff = e["outlier_threshold_cm"]
                threshold = cutoff if cutoff and stats["max_abs_cm"] < cutoff else None
                try:
                    run = self.run_duo(f"{label}_fit_{restart}_{retry}", self.p, jmax, e["iterations"],
                                       scale=scale, robust=robust, threshold=threshold)
                    proposed = run.parameters
                    candidate, new, check, fraction = self.validated_step(
                        proposed,stats,jmax,f"{label}_check_{restart}_{retry}")
                    # Acceptance is never based on a changing robust or outlier-rejected subset.
                    if new[metric] > stats[metric] + max(0.00005, stats[metric]*1e-6):
                        if run.history[-1]["stability"] <= e["stability_tolerance"]:
                            self.record(stage=label, accepted=False, diagnostics=new,
                                        reason="Weighted fit converged but worsened the acceptance metric; retained previous model")
                            return stats
                        raise FitError(f"Fixed-data {metric} worsened: {stats[metric]:.7g} -> {new[metric]:.7g}")
                    old_score = stats[metric]
                    self.p, stats, accepted = candidate, new, True
                    self.record(stage=label, accepted=True, fit_scale=scale, threshold=threshold,parameter_step_fraction=fraction,
                                diagnostics=new, validation_directory=str(check.directory.relative_to(self.root)),
                                history=run.history[-1])
                    print(f"  {label}: weighted RMS {new['weighted_rms_cm']:.6g}; RMS {new['rms_cm']:.6g} cm-1 ({new['count']} levels)", flush=True)
                    stable = fraction == 1 and run.history[-1]["stability"] <= e["stability_tolerance"]
                    improvement = (old_score-new[metric])/max(old_score, 1e-20)
                    if stable or improvement < e["minimum_relative_improvement"]:
                        return stats
                    scale = min(e["fit_scale"], scale*2)
                    break
                except FitError as exc:
                    self.record(stage=label, accepted=False, fit_scale=scale, reason=str(exc))
                    scale *= 0.5
            if not accepted:
                self.record(stage=label, reason="Retries exhausted; retained last validated parameters")
                return stats
        self.record(stage=label, reason="Restart limit reached; retained last validated parameters")
        return stats

    def validated_step(self,proposed,previous,jmax,label):
        """Backtrack an unsafe/worsening parameter update on fixed-data RMS."""
        e=self.config["experiment"];metric=e["target_metric"]
        last_error=None
        for step in range(e["line_search_steps"]+1):
            fraction=0.5**step
            p=proposed if step==0 else {k:self.p[k]+fraction*(value-self.p[k]) for k,value in proposed.items()}
            try:
                self.physical(p)
                new,_,check=self.evaluate(label+(f"_damp_{step}" if step else ""),p,jmax)
                if new[metric] <= previous[metric]+max(.00005,previous[metric]*1e-6):
                    return p,new,check,fraction
                # Preserve the stationary-objective handling when damping is disabled.
                if e["line_search_steps"]==0:
                    return p,new,check,fraction
                last_error=FitError(f"Fixed-data {metric} worsened: {previous[metric]:.7g} -> {new[metric]:.7g}")
            except FitError as exc:
                last_error=exc
            self.record(stage=label,accepted=False,parameter_step_fraction=fraction,reason=str(last_error))
        raise FitError(f"No acceptable damped step: {last_error}")

    def promote_order(self, jmax, stats):
        if kind(self.p)!="EMO":
            return stats
        e = self.config["experiment"]
        metric = e["target_metric"]
        while stats[metric] > e["target_rms_cm"] and self.p["NR"] < e["max_right_order"]:
            old_p, old_fitted = dict(self.p), list(self.fitted)
            old_stats = stats
            degree = min(int(self.p["NR"])+e["order_step"], e["max_right_order"])
            for k in range(int(max(self.p["NL"], self.p["NR"]))+1, max(int(self.p["NL"]), degree)+1):
                self.p[f"B{k}"] = 0.0
                self.fitted.append(f"B{k}")
            self.p["NR"] = float(degree)
            stats = self.fit_experiment(jmax, f"order_{degree}")
            if stats[metric] >= old_stats[metric]*(1-e["minimum_relative_improvement"]):
                self.p, self.fitted = old_p, old_fitted
                return old_stats
        return stats

    def execute(self):
        c, e, v = self.config, self.config["experiment"], self.config["validation"]
        if c["refine_only"]:
            self.physical(self.p)
            self.record(stage="refine_only",reason="Starting from the supplied validated refinement template")
        else:
            self.fit_abinitio()
        if c["experimental_de_cm"] is not None:
            # V0=0 enforced: configured experimental well depth equals Duo DE.
            de_key="AE" if kind(self.p)=="MLJ" else "DE"
            self.p[de_key] = c["experimental_de_cm"]
            self.fitted = [k for k in self.fitted if k != de_key]
            self.physical(self.p)
        for i, jmax in enumerate(c["jmax_schedule"]):
            if i == 0 and jmax>0 and e["bob_prefit"] and bob_keys(self.p):
                joint_fitted=list(self.fitted)
                bob_fitted=[k for k in self.fitted if k.startswith("BR_")]
                if bob_fitted:
                    try:
                        self.fitted=bob_fitted
                        self.fit_experiment(c["jmax_schedule"][-1],"bob_prefit")
                    finally:
                        self.fitted=joint_fitted
            stats = self.fit_experiment(jmax, f"experiment_J{jmax}", initial=i == 0)
            if i == 0 and jmax==0 and e["bob_prefit"] and bob_keys(self.p):
                joint_fitted=list(self.fitted)
                bob_fitted=[k for k in self.fitted if k.startswith("BR_")]
                if bob_fitted:
                    try:
                        self.fitted=bob_fitted
                        self.fit_experiment(c["jmax_schedule"][-1],"bob_prefit")
                    finally:
                        self.fitted=joint_fitted
        jmax = c["jmax_schedule"][-1]
        # The README also recommends disabling dynamic reweighting once close.
        # This polishes the objective defined by the original experimental weights.
        stats = self.fit_experiment(jmax, "final_polish", robust=e["polish_robust"])
        stats = self.promote_order(jmax, stats)
        converged = False
        delta = None
        for attempt in range(v["max_basis_refinements"]+1):
            stats, levels, final_run = self.evaluate(f"final_{attempt}", self.p, jmax)
            larger_vmax = self.vmax+v["vmax_increment"]
            lo, hi, n = self.grid
            larger_grid = (lo-v["range_extension_angstrom"][0], hi+v["range_extension_angstrom"][1], n+v["npoints_increment"])
            if larger_grid[0] <= 0 or larger_vmax > larger_grid[2]:
                raise FitError("Invalid expanded basis/grid")
            shape = shape_check(self.p, larger_grid, v["shape_tolerance_cm"],max_potential=v['potential_ceiling_cm'])
            if not shape["ok"]:
                raise FitError("Potential failed shape check on expanded grid")
            if not bob_shape(self.p,larger_grid,max_abs=v["bob_max_abs"],max_slope=v["bob_max_slope"],max_turns=v["bob_max_turns"],turning_tolerance=v["bob_turning_tolerance"])["ok"]:
                raise FitError("BOB failed shape check on expanded grid")
            _, larger, _ = self.evaluate(f"basis_check_{attempt}", self.p, jmax, grid=larger_grid, vmax=larger_vmax)
            differences = [{"J": o.j, "v": o.v, "energy_cm": levels[o.key]["calculated"],
                            "larger_basis_energy_cm": larger[o.key]["calculated"],
                            "difference_cm": larger[o.key]["calculated"]-levels[o.key]["calculated"]}
                           for o in self.model.observations]
            delta = max(abs(row["difference_cm"]) for row in differences)
            write_csv(self.root/"basis_convergence.csv", differences)
            self.record(stage="basis_check", max_delta_cm=delta, grid=self.grid, larger_grid=larger_grid,
                        vmax=self.vmax, larger_vmax=larger_vmax)
            if delta <= v["energy_tolerance_cm"]:
                converged = True
                break
            if attempt < v["max_basis_refinements"]:
                self.grid, self.vmax = larger_grid, larger_vmax
                stats = self.fit_experiment(jmax, f"basis_refit_{attempt}", robust=e["polish_robust"])
        quality = stats[e["target_metric"]] <= e["target_rms_cm"]
        status = "success" if quality and converged else "needs_review"
        final_text = (final_run.directory/"run.inp").read_text()
        (self.root/"final.inp").write_text(final_text, encoding="ascii")
        (self.root/"final_residuals.csv").write_bytes((final_run.directory/"residuals.csv").read_bytes())
        write_json(self.root/"final_parameters.json", self.p)
        shape = self.physical(self.p)
        lo, hi, _ = self.grid
        curve = [{"r_angstrom": lo+(hi-lo)*i/2000,
                  "potential_cm": potential(lo+(hi-lo)*i/2000, self.p)} for i in range(2001)]
        write_csv(self.root/"potential.csv", curve)
        if bob_keys(self.p):
            write_csv(self.root/"bob_rot.csv",[{"r_angstrom":row["r_angstrom"],"bob_dimensionless":bob(row["r_angstrom"],self.p)} for row in curve])
        result = {"status": status, "experiment": stats, "target_metric": e["target_metric"],
                  "target_cm": e["target_rms_cm"], "quality_pass": quality,
                  "basis_converged": converged, "basis_max_delta_cm": delta,
                  "basis_tolerance_cm": v["energy_tolerance_cm"], "shape": shape,
                  "Jmax": jmax, "vmax": self.vmax, "grid": self.grid,
                  "parameters": self.p, "validation_directory": str(final_run.directory.relative_to(self.root)),
                  "de_fixed_to_experiment": c["experimental_de_cm"] is not None}
        write_json(self.root/"result.json", result)
        self.report(result)
        return result

    def report(self, r):
        s = r["experiment"]
        text = f"""# Duo {kind(self.p)} fitting result: {r['status']}

Validated {s['count']} experimental levels through J={r['Jmax']}.

| Quantity | Result |
|---|---:|
| RMS over all positive-weight observations | {s['rms_cm']:.8g} cm-1 |
| RMS using original experimental weights | {s['weighted_rms_cm']:.8g} cm-1 |
| Largest absolute residual | {s['max_abs_cm']:.8g} cm-1 |
| Acceptance metric | {r['target_metric']} |
| Requested target | {r['target_cm']:.8g} cm-1 |
| Maximum change on increasing basis/grid | {r['basis_max_delta_cm']:.8g} cm-1 |
| Basis tolerance | {r['basis_tolerance_cm']:.8g} cm-1 |
| Single-well shape on specified radial interval | {r['shape']['ok']} |

`final.inp` is a self-contained Duo straight-through calculation (`itmax 0`).
`final_parameters.json`, `potential.csv`, `final_residuals.csv`, and
`basis_convergence.csv` contain the model and diagnostics. Every attempted fit,
including rejected steps, remains under `attempts/`; `progress.json` records decisions.
`manifest.json` records source, executable and pipeline SHA-256 hashes.

Validation retains every supplied observation; robust fitting weights and the
outlier cutoff do not define the reported acceptance metric. Zero input weights
are omitted from RMS statistics but retained in the residual table. Duo `.en`
energies are printed to 0.0001 cm-1; residuals and basis comparisons have that
resolution. A shape check on a finite interval does not validate extrapolation
to dissociation. This is a fitted PEC, not an independently validated line list.
"""
        (self.root/"report.md").write_text(text, encoding="utf-8")
