"""Duo EMO/MLJ and dimensionless rotational BOB conventions.

BOB parameters are namespaced BR_ internally so their RE/B0 cannot overwrite
the potential parameters when reading a multi-object fitted checkpoint.
"""
import math
import re

from .model import FitError, emo


def kind(p):
    return "MLJ" if "TE" in p else "EMO"


def potential_keys(p):
    if kind(p) == "EMO":
        return ["V0", "RE", "DE", "RREF", "PL", "PR", "NL", "NR"] + [f"B{i}" for i in range(int(max(p["NL"], p["NR"]))+1)]
    coefficients = sorted((k for k in p if re.fullmatch(r"PHI\d+", k)),key=lambda k:int(k[3:]))
    return ["TE", "RE", "AE", "R12", "DELTA", "PHIINF", "N"] + coefficients


def bob_keys(p):
    if "BR_RE" not in p:
        return []
    coefficients = sorted((k for k in p if re.fullmatch(r"BR_B\d+", k)),key=lambda k:int(k[4:]))
    return ["BR_RE", "BR_BETA", "BR_GAMMA", "BR_P"] + coefficients + ["BR_BINF"]


def asymptote(p):
    return p["AE"] if kind(p) == "MLJ" else p["DE"]


def potential(r,p):
    if kind(p) == "EMO":
        return emo(r,p)
    z = 2*(r-p["RE"])/(r+p["RE"])
    polynomial = sum(p[k]*z**int(k[3:]) for k in p if re.fullmatch(r"PHI\d+",k))
    x = p["DELTA"]*(r-p["R12"])
    switch = math.exp(-x)/(1+math.exp(-x)) if x > 0 else 1/(1+math.exp(x))
    phi = switch*polynomial+(1-switch)*p["PHIINF"]
    try:
        y = -math.expm1(p["N"]*math.log(p["RE"]/r)-phi*z)
        return p["TE"]+(p["AE"]-p["TE"])*y*y
    except OverflowError:
        return math.inf


def bob(r,p):
    if "BR_RE" not in p:
        return 0.0
    dr = r-p["BR_RE"]
    z = dr*math.exp(-p["BR_BETA"]*dr*dr-p["BR_GAMMA"]*dr**4)
    power = p["BR_P"]
    y = (r**power-p["BR_RE"]**power)/(r**power+p["BR_RE"]**power)
    f = sum(p[k]*z**int(k[4:]) for k in p if re.fullmatch(r"BR_B\d+", k))
    return (1-y)*f+y*p["BR_BINF"]


def validate_parameters(p,fitted):
    if not p or not all(math.isfinite(v) for v in p.values()):
        raise FitError("All model parameters must be finite")
    if kind(p) == "EMO":
        required = ["V0", "RE", "DE", "RREF", "PL", "PR", "NL", "NR"]
        if any(k not in p for k in required):
            raise FitError("Missing EMO parameter")
        for k in ("PL","PR","NL","NR"):
            if p[k] != int(p[k]) or p[k] < (1 if k.startswith("P") else 0):
                raise FitError(f"Invalid EMO {k}")
        if p.get("B0",0) <= 0 or p["V0"] != 0:
            raise FitError("EMO requires positive B0 and fixed V0=0")
        if all(f"B{i}" in p for i in range(int(p["NR"])+1)) and sum(p[f"B{i}"] for i in range(int(p["NR"])+1))<=0:
            raise FitError("EMO asymptotic exponent must be positive to approach DE")
        allowed = {"RE","DE"} | {k for k in potential_keys(p) if re.fullmatch(r"B\d+",k)}
    else:
        required = ["TE","RE","AE","R12","DELTA","PHIINF","N","PHI0"]
        if any(k not in p for k in required):
            raise FitError("Missing MLJ parameter")
        coeffs = [k for k in potential_keys(p) if re.fullmatch(r"PHI\d+",k)]
        if coeffs != [f"PHI{i}" for i in range(len(coeffs))]:
            raise FitError("MLJ coefficients must be contiguous PHI0..PHIn")
        if p["TE"] != 0 or p["DELTA"] <= 0 or p["R12"] <= p["RE"] or p["N"] < 1 or p["N"] != int(p["N"]):
            raise FitError("MLJ requires fixed TE=0, positive DELTA/integer N and R12>RE")
        allowed = {"RE","AE"} | set(coeffs)
    if p["RE"] <= 0 or asymptote(p) <= 0:
        raise FitError("Positive equilibrium radius and well depth required")
    bk = bob_keys(p)
    if bk:
        if any(k not in p for k in bk) or not any(k == "BR_B0" for k in bk):
            raise FitError("Missing BOB parameter")
        coeffs = [k for k in bk if re.fullmatch(r"BR_B\d+",k)]
        if coeffs != [f"BR_B{i}" for i in range(len(coeffs))]:
            raise FitError("BOB coefficients must be contiguous B0..Bn")
        if p["BR_RE"] <= 0 or p["BR_BETA"] < 0 or p["BR_GAMMA"] < 0 or p["BR_P"] < 1 or p["BR_P"] != int(p["BR_P"]) or p["BR_BINF"] != 0:
            raise FitError("BOB requires positive RE/integer P, nonnegative damping, fixed BINF=0")
        allowed |= set(coeffs)
    if set(p) != set(potential_keys(p)+bk):
        raise FitError("Unexpected or missing potential/BOB coefficients")
    if not fitted or len(set(fitted)) != len(fitted) or not set(fitted) <= allowed:
        raise FitError("Only well depth, RE, potential coefficients and BOB B coefficients may be fitted")


def bob_shape(p,grid,*,max_abs=0.1,max_slope=0.5,max_turns=4,turning_tolerance=0.0,points=2001):
    if not bob_keys(p):
        return {"ok":True,"present":False}
    lo,hi=grid[:2]
    dr=(hi-lo)/(points-1)
    vals=[bob(lo+i*dr,p) for i in range(points)]
    if not all(math.isfinite(v) for v in vals):
        return {"ok":False,"present":True,"reason":"nonfinite BOB"}
    slopes=[(b-a)/dr for a,b in zip(vals,vals[1:])]
    # Ignore derivative roundoff in an effectively flat tail.
    signs=[1 if s>0 else -1 for s in slopes if abs(s)>1e-7]
    raw_turns=sum(a!=b for a,b in zip(signs,signs[1:]))
    turns=raw_turns
    if turning_tolerance>0:
        # Count reversals only after a specified dimensionless excursion. This
        # avoids treating negligible ripples in a decaying tail as extra lobes.
        turns=0;direction=0;extreme=vals[0]
        for value in vals[1:]:
            if direction==0:
                if abs(value-extreme)>turning_tolerance:
                    direction=1 if value>extreme else -1;extreme=value
            elif direction*(value-extreme)>0:
                extreme=value
            elif direction*(extreme-value)>turning_tolerance:
                turns+=1;direction=-direction;extreme=value
    amplitude=max(map(abs,vals)); slope=max(map(abs,slopes))
    ok=amplitude<=max_abs and slope<=max_slope and turns<=max_turns and min(vals)>-1
    return {"ok":ok,"present":True,"max_abs":amplitude,"max_slope_per_angstrom":slope,
            "turning_points":turns,"raw_turning_points":raw_turns,"minimum_rotational_factor":1+min(vals),
            "limits":{"max_abs":max_abs,"max_slope":max_slope,"max_turns":max_turns,"turning_tolerance":turning_tolerance},
            "reason":"bounded smooth correction" if ok else "BOB amplitude, slope or oscillation limit exceeded"}


def centrifugal_shape(p,masses,jmax,grid,points=4001):
    """One inner well and at most one centrifugal barrier for each integer J.

    A monotonic PEC may still have a shoulder creating extra rotational wells.
    The factor is h/(8*pi*pi*c*u*Angstrom**2), in cm-1 Angstrom**2 u.
    This shape diagnostic does not assert that above-threshold states are bound.
    """
    if not masses:
        raise FitError("Centrifugal shape checks require explicit isotope masses")
    m1,m2=map(float,masses.split());mu=m1*m2/(m1+m2)
    lo,hi=grid[:2];step=(hi-lo)/(points-1)
    rr=[lo+i*step for i in range(points)]
    vv=[potential(r,p) for r in rr]
    rot=[16.857629206/mu*(1+bob(r,p))/r**2 for r in rr]
    failures=[]
    for j in range(jmax+1):
        vals=[v+j*(j+1)*b for v,b in zip(vv,rot)]
        directions=[(i,1 if b>a else -1) for i,(a,b) in enumerate(zip(vals,vals[1:])) if abs(b-a)>1e-5]
        minima=[rr[k] for (i,a),(k,b) in zip(directions,directions[1:]) if a<0 and b>0]
        maxima=[rr[k] for (i,a),(k,b) in zip(directions,directions[1:]) if a>0 and b<0]
        if len(minima)!=1 or len(maxima)>1:
            failures.append({"J":j,"minima_angstrom":minima,"maxima_angstrom":maxima})
    return {"ok":not failures,"range_angstrom":[lo,hi],"failures":failures,
            "criterion":"one centrifugal well and at most one barrier for each observed J"}
