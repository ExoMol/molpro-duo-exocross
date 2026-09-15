"""Run Duo without a shell; read full-precision parameters and diagnostic files."""
from __future__ import annotations

from dataclasses import dataclass
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

from .model import FitError, number


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, data):
    path = Path(path)
    temp = path.with_suffix(path.suffix+".tmp")
    temp.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    temp.replace(path)


def parse_output(text):
    iterations = []
    # Only the Parameters section, never the echoed input or rounded table.
    # Robust-fit builds print Watson diagnostics between Iteration and Parameters.
    pattern = r"(?ms)^Iteration\s*=\s*(\d+)\s*\n(?:(?!^Iteration).)*?^Parameters:\s*\n(.*?)(?=^Fitted parameters \(rounded\):)"
    for m in re.finditer(pattern, text):
        fields = list(re.finditer(r"(?ims)^POTEN\s+[^\n]+\n.*?^Values\s*\n(.*?)^end\s*$", m[2]))
        if len(fields) != 1:
            raise FitError("Duo output does not contain exactly one fitted potential")
        p = {}
        for line in fields[0][1].splitlines():
            row = line.split()
            if row:
                p[row[0].upper()] = number(row[1])
        bob_fields = list(re.finditer(r"(?ims)^(?:BOBROT|BOB-ROT)\s+[^\n]+\n.*?^Values\s*\n(.*?)^end\s*$",m[2]))
        if len(bob_fields)>1:
            raise FitError("Only one rotational BOB field is supported")
        for field in bob_fields:
            for line in field[1].splitlines():
                row=line.split()
                if row:
                    p["BR_"+row[0].upper()]=number(row[1])
        iterations.append({"iteration": int(m[1]), "parameters": p})
    history = []
    for line in text.splitlines():
        if line.startswith("-->|"):
            row = [s.strip() for s in line.split("|")[1:-1]]
            if len(row) == 7:
                history.append(dict(zip(("iteration", "points", "nparameters", "duo_weighted_rms",
                                         "duo_energy_rms", "duo_pec_rms", "stability"), map(number, row))))
    return iterations, history


def converged_checkpoint(text, tolerance, minimum_iterations=3):
    """Select a complete parameter/summary pair after three stable iterations."""
    iterations, history = parse_output(text)
    if len(history) < 3 or any(h["stability"] > tolerance for h in history[-3:]):
        return None
    last = int(history[-1]["iteration"])
    if last < minimum_iterations:
        return None
    if not any(i["iteration"] == last for i in iterations):
        return None
    return {"iteration": last, "reason": "stability", "tolerance": tolerance,
            "minimum_iterations": minimum_iterations}


def parse_en(path):
    text = Path(path).read_text(errors="replace")
    sections = re.split(r"(?m)^\s*Iteration\s*=\s*\d+", text)
    result = {}
    for line in sections[-1].splitlines():
        if "(" not in line:
            continue
        prefix, *groups = re.split(r"[()]", line)
        row = prefix.split()
        if len(row) != 8:
            raise FitError(f"Malformed Duo energy row: {line}")
        calc = groups[0].split()
        obs = groups[2].split() if len(groups) > 2 and groups[2].strip() else None
        if len(calc) != 5 or (obs is not None and len(obs) != 5):
            raise FitError(f"Malformed Duo assignments: {line}")
        j = number(row[2])
        if int(j) != j:
            raise FitError("Half-integer J is outside single singlet-state scope")
        key = int(j), int(calc[1])
        if key in result:
            raise FitError(f"Duplicate calculated (J,v): {key}")
        result[key] = {"j": key[0], "v": key[1], "n": int(row[0]), "matched_n": int(row[1]),
                       "parity": row[3], "observed": number(row[4]), "calculated": number(row[5]),
                       "residual": number(row[6]), "duo_weight": number(row[7]),
                       "calculated_assignment": calc, "observed_assignment": obs,
                       "mark": groups[-1].strip()}
    if not result:
        raise FitError(f"No calculated levels in {path}")
    return result


def assess_energies(levels, observations, jmax):
    rows, errors = [], []
    for o in observations:
        if o.j > jmax:
            continue
        c = levels.get(o.key)
        if c is None:
            errors.append(f"Missing J={o.j}, v={o.v}")
            continue
        expected = [o.state, str(o.v), "0", "0.0", "0.0"]
        # Compare numerical quantum numbers, not their string formatting.
        def same(q):
            return q is not None and q[0] == expected[0] and all(number(a) == number(b) for a,b in zip(q[1:], expected[1:]))
        if (not same(c["calculated_assignment"]) or not same(c["observed_assignment"])
                or c["n"] != o.v+1 or c["matched_n"] != o.v+1
                or c["parity"] != ("+" if o.j % 2 == 0 else "-")
                or c["mark"] or abs(c["observed"]-o.energy) > 0.00011):
            errors.append(f"Assignment mismatch J={o.j}, v={o.v}")
        # Recompute using original observation precision and fixed input weights.
        rows.append({"J": o.j, "v": o.v, "observed_cm": o.energy, "calculated_cm": c["calculated"],
                     "residual_cm": o.energy-c["calculated"], "input_weight": o.weight,
                     "duo_weight": c["duo_weight"]})
    if errors:
        raise FitError("; ".join(errors[:8])+f" ({len(errors)} assignment/coverage errors)")
    active = [r for r in rows if r["input_weight"] > 0]
    if not active:
        raise FitError("No positive-weight observations in validation")
    sq = sum(r["residual_cm"]**2 for r in active)
    wsq = sum(r["input_weight"]*r["residual_cm"]**2 for r in active)
    return {"count": len(rows), "positive_weight_count": len(active),
            "rms_cm": math.sqrt(sq/len(active)),
            "weighted_rms_cm": math.sqrt(wsq/sum(r["input_weight"] for r in active)),
            "max_abs_cm": max(abs(r["residual_cm"]) for r in active),
            "excluded_by_duo_count": sum(r["duo_weight"] == 0 for r in active)}, rows


def write_csv(path, rows):
    if not rows:
        raise FitError(f"No rows for {path}")
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


@dataclass
class RunResult:
    directory: Path
    parameters: dict | None
    history: list
    controlled_stop: dict | None = None


class Runner:
    def __init__(self, executable, root, *, timeout=1800, threads=1, resume=False, poll_seconds=10):
        resolved = shutil.which(str(executable)) or str(Path(executable).resolve())
        self.executable = Path(resolved)
        if not self.executable.is_file():
            raise FitError(f"Duo executable not found: {executable}. Set --duo or DUO_EXE.")
        self.executable_hash = digest(self.executable)
        self.root = Path(root)
        self.timeout, self.threads, self.resume = timeout, threads, resume
        self.poll_seconds = poll_seconds

    def run(self, name, text, fitting=False, stability_tolerance=None, minimum_iterations=3):
        directory = self.root/name
        directory.mkdir(parents=True, exist_ok=True)
        signature = hashlib.sha256((text+self.executable_hash+str(self.threads)).encode()).hexdigest()
        stamp = directory/"completed.json"
        files = ("run.out", "fit.en", "fit.pot", "Potential_functions.dat")
        cached = False
        controlled_stop = None
        if self.resume and stamp.exists():
            data = json.loads(stamp.read_text())
            cached = (data.get("signature") == signature and
                      all((directory/f).is_file() and digest(directory/f) == data.get("artifacts", {}).get(f) for f in files))
            if cached:
                controlled_stop = data.get("controlled_stop")
                if controlled_stop and (stability_tolerance is None or controlled_stop["tolerance"] > stability_tolerance
                                        or controlled_stop["iteration"] < minimum_iterations):
                    cached = False
                    controlled_stop = None
        if not cached:
            # Remove only known products in this specific attempt directory.
            for f in (*files, "completed.json"):
                (directory/f).unlink(missing_ok=True)
            (directory/"run.inp").write_text(text, encoding="ascii")
            env = os.environ.copy()
            env.update(OMP_NUM_THREADS=str(self.threads), MKL_NUM_THREADS=str(self.threads), OPENBLAS_NUM_THREADS=str(self.threads))
            start = time.monotonic()
            last_notice = start
            print(f"Running {name}", flush=True)
            with (directory/"run.inp").open("rb") as inp, (directory/"run.out").open("wb") as out, (directory/"stderr.txt").open("wb") as err:
                process = subprocess.Popen([str(self.executable)], stdin=inp, stdout=out, stderr=err,
                                           cwd=directory, env=env, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
                try:
                    while True:
                        try:
                            code = process.wait(timeout=min(self.poll_seconds, self.timeout))
                            break
                        except subprocess.TimeoutExpired:
                            elapsed = time.monotonic()-start
                            if elapsed >= self.timeout:
                                raise FitError(f"Duo timed out after {elapsed:.0f}s in {directory}")
                            if fitting and stability_tolerance is not None:
                                # The tail includes several whole parameter/summary sections
                                # without repeatedly reading a large spectral output file.
                                with (directory/"run.out").open("rb") as snapshot:
                                    snapshot.seek(max(0, snapshot.seek(0, 2)-8*1024*1024))
                                    tail = snapshot.read().decode(errors="replace")
                                controlled_stop = converged_checkpoint(tail, stability_tolerance, minimum_iterations)
                                if controlled_stop:
                                    process.terminate()
                                    code = process.wait()
                                    print(f"  {name}: reached stability at iteration {controlled_stop['iteration']}; validating checkpoint", flush=True)
                                    break
                            if time.monotonic()-last_notice >= 30:
                                print(f"  {name}: {elapsed:.0f}s elapsed", flush=True)
                                last_notice = time.monotonic()
                finally:
                    if process.poll() is None:
                        process.kill()
                        process.wait()
            if code and not controlled_stop:
                hint = " Missing Windows DLL/runtime." if code & 0xffffffff == 0xc0000135 else ""
                raise FitError(f"Duo exited {code} in {directory}.{hint} See stderr.txt and run.out.")
        else:
            print(f"Reusing {name}", flush=True)
        output = (directory/"run.out").read_text(errors="replace")
        error = re.search(r"(?im)^.*(?:forrtl:|segmentation fault|error:|illegal number|NaN|Infinity).*$", output)
        if error:
            raise FitError(f"Duo reported an error in {directory}: {error[0][:240]}")
        if not controlled_stop and ("Timing data at" not in output[-7000:] or "Zero point energy" not in output):
            raise FitError(f"Duo output is incomplete in {directory}")
        iterations, history = parse_output(output)
        if controlled_stop:
            checkpoint = controlled_stop["iteration"]
            iterations = [i for i in iterations if i["iteration"] <= checkpoint]
            history = [h for h in history if h["iteration"] <= checkpoint]
            if (not iterations or not history or iterations[-1]["iteration"] != checkpoint
                    or history[-1]["iteration"] != checkpoint
                    or len(history) < 3 or any(h["stability"] > controlled_stop["tolerance"] for h in history[-3:])):
                raise FitError(f"Incomplete or invalid convergence checkpoint in {directory}")
        if fitting and (not iterations or not history or iterations[-1]["iteration"] != history[-1]["iteration"]):
            raise FitError(f"No complete fitted parameter iteration in {directory}")
        if not all((directory/f).is_file() and (directory/f).stat().st_size > 0 for f in files):
            raise FitError(f"Duo omitted required diagnostic files in {directory}")
        if not cached:
            write_json(stamp, {"signature": signature, "seconds": time.monotonic()-start,
                               "artifacts": {f: digest(directory/f) for f in files},
                               "controlled_stop": controlled_stop})
        return RunResult(directory, iterations[-1]["parameters"] if iterations else None, history, controlled_stop)
