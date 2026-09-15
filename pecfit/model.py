"""Validated input data and the EMO convention used by Duo/functions.f90."""
from __future__ import annotations

from dataclasses import dataclass, replace
import math
from pathlib import Path
import re

HARTREE_TO_CM = 219474.6313708


class FitError(ValueError):
    pass


def number(token: str) -> float:
    normalized = token.replace("D", "E").replace("d", "e")
    # Fortran may omit E to fit a three-digit exponent into a fixed field.
    normalized = re.sub(r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+))([+-]\d{3})$", r"\1E\2", normalized)
    try:
        value = float(normalized)
    except ValueError as exc:
        raise FitError(f"Invalid or overflowed numerical field: {token}") from exc
    if not math.isfinite(value):
        raise FitError(f"Non-finite number: {token}")
    return value


def clean(line: str) -> str:
    return re.split(r"[!#(]", line, maxsplit=1)[0].strip()


def block(text: str, header: str) -> re.Match:
    matches = list(re.finditer(r"(?im)^\s*" + header + r"[^\n]*\n.*?^\s*end\s*(?:[!#].*)?$", text, re.S))
    if len(matches) != 1:
        raise FitError(f"Expected exactly one {header!r} block; found {len(matches)}")
    return matches[0]


def values(text: str) -> list[str]:
    parts = re.split(r"(?im)^\s*values\s*\n", text, maxsplit=1)
    if len(parts) != 2:
        raise FitError("Missing VALUES")
    return [clean(s) for s in parts[1].splitlines() if clean(s) and clean(s).lower() != "end"]


@dataclass(frozen=True)
class Observation:
    j: int
    parity: str
    n: int
    energy: float
    state: str
    v: int
    weight: float

    @property
    def key(self):
        return self.j, self.v

    def line(self):
        return f"{self.j:5d} e {self.v+1:5d} {self.energy:.12f} {self.state} {self.v:4d} 0 0 0 {self.weight:.12g}"


def observations(lines: list[str]) -> list[Observation]:
    result = []
    seen = set()
    for line in lines:
        row = clean(line).split()
        if not row or row[0].lower() == "end":
            continue
        if len(row) != 10:
            raise FitError(f"Expected 10 energy columns, got: {line}")
        j, n, v = number(row[0]), number(row[2]), number(row[5])
        if any(x < 0 or int(x) != x for x in (j, v)) or n < 1 or int(n) != n:
            raise FitError(f"Invalid J, N or v: {line}")
        # For a single Sigma+ singlet all levels are e. +/- is total parity.
        parity = row[1].lower()
        allowed = "+" if int(j) % 2 == 0 else "-"
        if parity not in ("e", allowed) or any(number(x) != 0 for x in row[6:9]):
            raise FitError(f"Only X1Sigma+ levels (e, Lambda=Sigma=Omega=0) are supported: {line}")
        obs = Observation(int(j), "e", int(n), number(row[3]), row[4], int(v), number(row[9]))
        if obs.weight < 0 or obs.energy < 0 or obs.key in seen:
            raise FitError(f"Negative weight/energy or duplicate (J,v): {line}")
        if obs.n != obs.v + 1:
            raise FitError(f"Single-state N must equal v+1; review this assignment: {line}")
        seen.add(obs.key)
        result.append(obs)
    result.sort(key=lambda o: o.key)
    if not result or not any(o.weight > 0 and o.energy > 0 for o in result):
        raise FitError("No positive-weight nonzero experimental energies")
    return result


def read_pec(path: Path, units: str, shift: bool) -> list[tuple[float, float]]:
    rows = []
    for line in path.read_text().splitlines():
        tokens = clean(line).replace(",", " ").split()
        if not tokens:
            continue
        if len(tokens) != 2:
            raise FitError(f"PEC requires two columns r/Angstrom, energy: {line}")
        rows.append(tuple(number(x) for x in tokens))
    return convert_pec(rows, units, shift)


def convert_pec(rows, units, shift):
    if units.lower() not in ("eh", "hartree", "cm-1"):
        raise FitError("PEC units must be hartree, Eh or cm-1")
    if len(rows) < 3 or any(r <= 0 for r, _ in rows):
        raise FitError("PEC requires at least three positive radii")
    if any(b[0] <= a[0] for a, b in zip(rows, rows[1:])):
        raise FitError("PEC radii must be strictly increasing")
    offset = min(e for _, e in rows) if shift else 0.0
    factor = HARTREE_TO_CM if units.lower() in ("eh", "hartree") else 1.0
    return [(r, (e-offset)*factor) for r, e in rows]


def emo(r: float, p: dict[str, float]) -> float:
    """Duo EMO: DE is the asymptote, so depth = DE - V0 (cm-1)."""
    ref = p["RREF"] if p["RREF"] > 0 else p["RE"]
    side = "L" if r <= ref else "R"
    power, degree = int(p["P"+side]), int(p["N"+side])
    z = (r**power-ref**power)/(r**power+ref**power)
    beta = sum(p[f"B{k}"]*z**k for k in range(degree+1))
    try:
        return p["V0"] + (p["DE"]-p["V0"])*math.expm1(-beta*(r-p["RE"]))**2
    except OverflowError:
        return math.inf


@dataclass
class Model:
    atoms: str
    masses: str | None
    state: str
    name: str
    grid: tuple[float, float, int]
    parameters: dict[str, float]
    fitted: list[str]
    pec: list[tuple[float, float]]
    weighting: str
    observations: list[Observation]

    def validate(self):
        p = self.parameters
        from .forms import validate_parameters
        validate_parameters(p, self.fitted)
        lo, hi, n = self.grid
        if not 0 < lo < p["RE"] < hi or n < 20:
            raise FitError("Invalid radial grid or RE outside grid")
        if len(self.atoms.split()) != 2 or any(o.state != self.state for o in self.observations):
            raise FitError("Exactly two atoms and one experimental electronic state are required")
        if len(self.state) > 3:
            raise FitError("Use a state label of at most three characters (Duo prints three-character tags)")
        if self.masses:
            mass_values = clean(self.masses).split()
            if len(mass_values) != 2 or any(number(m) <= 0 for m in mass_values):
                raise FitError("Supply two positive finite isotope masses")
        if not any(o.j == 0 and o.v == 0 and abs(o.energy) < 1e-8 for o in self.observations):
            raise FitError("Include the J=0,v=0 zero-energy reference")


def import_template(path: Path) -> Model:
    text = path.read_text(encoding="utf-8-sig")
    unsupported = re.search(r"(?im)^\s*(?:bob-vib|bobvib|spin-orbit|spinorbit|spin-spin|spinspin|spin-rot|spinrot|diabatic|lambda-doubling|lambdaopq|lambdap2q|lambdaq|L2|LxLy|NAC)\b", text)
    if unsupported:
        raise FitError(f"Additional Hamiltonian fields are outside the one-PEC scope: {unsupported[0].strip()}")
    def keyword(key):
        m = re.search(r"(?im)^\s*" + key + r"\s+([^\n]+)", text)
        return m.group(1).strip() if m else None
    states = keyword("states")
    if not states or len(clean(states).split()) != 1:
        raise FitError("Template must contain exactly one electronic state")
    state = clean(states)
    analytic = block(text, r"poten\b").group()
    for key, expected in (("lambda", "0"), ("mult", "1"), ("symmetry", "+")):
        match = re.search(r"(?im)^\s*"+key+r"\s+(\S+)", analytic)
        if not match or (number(match[1]) != number(expected) if key in ("lambda","mult") else match[1] != expected):
            raise FitError(f"POTEN must have {key} {expected}")
    form = re.search(r"(?im)^\s*type\s+(\S+)",analytic)
    if not form or form[1].upper() not in ("EMO","MLJ"):
        raise FitError("POTEN type must be EMO or MLJ")
    p, fitted = {}, []
    for line in values(analytic):
        row = line.split()
        if len(row) not in (2, 3) or (len(row) == 3 and row[2].lower() != "fit"):
            raise FitError(f"Unsupported parameter line: {line}")
        key = row[0].upper()
        if key in p:
            raise FitError(f"Duplicate parameter: {key}")
        p[key] = number(row[1])
        if len(row) == 3:
            fitted.append(key)
    from .forms import kind
    if kind(p) != form[1].upper():
        raise FitError("Potential parameter names do not match declared type")
    if re.search(r"(?im)^\s*(?:bob-rot|bobrot)\s",text):
        br = block(text,r"(?:bob-rot|bobrot)\b").group()
        for key, expected in (("type","POLYNOM_DECAY_24"),("lambda","0"),("spin","0")):
            match = re.search(r"(?im)^\s*"+key+r"\s+(\S+)",br)
            if not match or match[1].upper()!=expected:
                raise FitError(f"BOB-ROT must have {key} {expected}")
        if re.search(r'(?im)^\s*bob-rot\s+(\S+)',br)[1] != state:
            raise FitError("BOB-ROT must belong to the single electronic state")
        for line in values(br):
            row=line.split()
            if len(row) not in (2,3) or (len(row)==3 and row[2].lower()!="fit"):
                raise FitError(f"Unsupported BOB parameter line: {line}")
            key="BR_"+row[0].upper()
            if key in p:
                raise FitError(f"Duplicate BOB parameter: {key}")
            p[key]=number(row[1])
            if len(row)==3:
                fitted.append(key)
    abi = block(text, r"abinitio\s+poten\b").group()
    units = re.search(r"(?im)^\s*units\s+(\S+)", abi)
    if not units:
        raise FitError("Explicit units are required in abinitio POTEN")
    rows = [tuple(number(x) for x in line.split()) for line in values(abi)]
    if any(len(row) != 2 for row in rows):
        raise FitError("Only two-column ab initio grids are supported")
    pec = convert_pec(rows, units.group(1), units.group(1).lower() in ("eh", "hartree"))
    weighting = re.search(r"(?im)^\s*weighting\s+([^\n]+)", abi)
    grid = block(text, r"grid\b").group()
    rg = re.search(r"(?im)^\s*range\s+(\S+)\s+(\S+)", grid)
    np = re.search(r"(?im)^\s*npoints\s+(\d+)", grid)
    if not rg or not np:
        raise FitError("GRID requires range and npoints")
    fit = block(text, r"fitting\b").group()
    parts = re.split(r"(?im)^\s*energies[^\n]*\n", fit, maxsplit=1)
    if len(parts) != 2:
        raise FitError("FITTING requires an ENERGIES table")
    name = re.search(r'(?im)^\s*name\s+([^\n]+)', analytic)
    model = Model(clean(keyword("atoms") or ""), keyword("masses"), state,
                  name.group(1).strip().strip('"') if name else "X1Sigma+",
                  (number(rg[1]), number(rg[2]), int(np[1])), p, fitted, pec,
                  clean(weighting[1]) if weighting else "PS1997 1e-5 80000.0",
                  observations(parts[1].splitlines()))
    model.validate()
    return model


def parameter_lines(p, fitted):
    from .forms import potential_keys
    keys = potential_keys(p)
    return "\n".join(f"{k:8s} {p[k]:.14E}" + (" fit" if k in fitted else "") for k in keys)


def render(model: Model, p, *, jmax, vmax, iterations, scale=0.1, robust=0.0,
           abi_factor=1e-7, energy_factor=1e6, threshold=None, fitted=None, grid=None):
    lo, hi, n = grid or model.grid
    from .forms import kind, bob_keys
    obs = [o for o in model.observations if o.j <= jmax and o.v < vmax]
    # Duo's vmax counts retained functions: vmax=1 retains v=0 only.
    if not obs:
        raise FitError("No observations in this stage")
    text = [f"atoms {model.atoms}"]
    if model.masses:
        text.append(f"masses {model.masses}")
    text += [f"states {model.state}", "print_pecs_and_couplings_to_file",
             "assign_V_by_count", f"jrot 0 - {jmax}", "symmetry Cs(M)",
             f"grid\n npoints {n}\n range {lo:.12g} {hi:.12g}\n type 0\nend",
             "diagonalizer\n SYEV\nend", f"contraction\n vib\n vmax {vmax}\nend",
             f'poten {model.state}\nname "{model.name}"\nlambda 0\nsymmetry +\nmult 1\ntype {kind(p)}\nvalues',
             parameter_lines(p, model.fitted if fitted is None else fitted), "end",
             f'abinitio poten {model.state}\nname "{model.name}"\nlambda 0\nsymmetry +\nmult 1\ntype grid\nunits cm-1',
             f"fit_factor {abi_factor:.12g}\nweighting {model.weighting}\nvalues",
             "\n".join(f"{r:.12g} {e:.14g}" for r, e in model.pec), "end",
             f"fitting\nJLIST 0 - {jmax}\nitmax {iterations}\nfit_factor {energy_factor:.12g}",
             f"fit_type DGELSS\nfit_scale {scale:.12g}\noutput fit\nlock -8000\nrobust {robust:.12g}"]
    if bob_keys(p):
        fit_keys = model.fitted if fitted is None else fitted
        br=[f'bob-rot {model.state}\nname "<X|BR|X>"\nlambda 0\nspin 0\ntype polynom_decay_24\nvalues']
        br += [f"{k[3:]:8s} {p[k]:.14E}"+(" fit" if k in fit_keys else "") for k in bob_keys(p)]
        br.append("end")
        text.insert(next(i for i,s in enumerate(text) if s.startswith("abinitio poten")),"\n".join(br))
    if threshold is not None:
        text.append(f"THRESH_OBS-CALC {threshold:.12g}")
    text += ["energies", "\n".join(o.line() for o in obs), "end", ""]
    return "\n".join(text)


def shape_check(p, grid, tolerance=0.01, points=2001, max_potential=None):
    """Finite, single-well, monotone branches on the declared radial interval."""
    lo, hi = grid[:2]
    from .forms import asymptote, potential
    if not lo < p["RE"] < hi or asymptote(p) <= 0:
        return {"ok": False, "reason": "invalid depth, RE or B0"}
    radii = sorted(set([lo+(hi-lo)*i/(points-1) for i in range(points)] + [p["RE"]]))
    curve = [(r, potential(r, p)) for r in radii]
    if not all(math.isfinite(v) for _, v in curve):
        return {"ok": False, "reason": "non-finite PEC"}
    highest=max(v for _,v in curve)
    if max_potential is not None and highest>max_potential:
        return {"ok":False,"reason":"repulsive wall exceeds configured ceiling","maximum_potential_cm":highest,
                "potential_ceiling_cm":max_potential,"range_angstrom":[lo,hi]}
    bad = [(r, s) for (r, v), (s, w) in zip(curve, curve[1:])
           if (s <= p["RE"] and w > v+tolerance) or (r >= p["RE"] and w < v-tolerance)]
    return {"ok": not bad, "reason": "monotone single well" if not bad else "extra turning points",
           "violating_intervals": len(bad), "range_angstrom": [lo, hi], "points": len(curve),"maximum_potential_cm":highest}
