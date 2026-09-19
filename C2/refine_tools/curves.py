"""Curve diagnostics using the analytic definitions in Duo/functions.f90."""
import math
from coupled import fields


def value(kind, p, r):
    if kind == 'EMO':
        te, re, de, rref, pl, pr, nl, nr = p[:8]
        if rref <= 0:
            rref = re
        power, order = (pl, nl) if r <= rref else (pr, nr)
        z = (r**power-rref**power)/(r**power+rref**power)
        beta = sum(p[8+k]*z**k for k in range(int(order)+1))
        return te+(de-te)*(1-math.exp(-beta*(r-re)))**2
    if kind == 'TWO_COUPLED_EMOS':
        n1 = int(p[7])+9
        n2 = int(p[n1+7])+9
        a, b = value('EMO', p[:n1], r), value('EMO', p[n1:n1+n2], r)
        q = p[n1+n2:-1]
        x = r-q[2]
        coupling = q[0]+sum(c*x**k for k,c in enumerate(q[3:]))/math.cosh(q[1]*x)
        split = math.hypot(a-b, 2*coupling)
        return (a+b+(-1 if int(p[-1]) == 1 else 1)*split)/2
    if kind == 'POLYNOM_DECAY_24':
        re, beta, gamma, power = p[:4]
        d = r-re
        z = d*math.exp(-beta*d*d-gamma*d**4)
        y = (r**power-re**power)/(r**power+re**power)
        return (1-y)*sum(c*z**k for k,c in enumerate(p[4:-1]))+y*p[-1]
    if kind == 'BOBLEROY':
        re, _, power, _ = p[:4]
        z = (r**power-re**power)/(r**power+re**power)
        return (1-z)*sum(c*z**k for k,c in enumerate(p[4:-1]))+z*p[-1]
    raise ValueError(f'Unsupported analytic field: {kind}')


def turning_points(values, tolerance=1e-6):
    signs = [1 if b-a > tolerance else -1 for a,b in zip(values, values[1:]) if abs(b-a) > tolerance]
    return sum(a != b for a,b in zip(signs, signs[1:]))


def audit(candidate, baseline):
    current, original = fields(candidate), fields(baseline)
    radius = [.7+i*(6-.7)/2000 for i in range(2001)]
    result = {}
    for key, field in current.items():
        old = original[key]
        try:
            y = [value(field.kind, field.values, r) for r in radius]
            y0 = [value(old.kind, old.values, r) for r in radius]
            if not all(math.isfinite(x) for x in y):
                raise ValueError('Non-finite curve')
            delta = [abs(a-b) for a,b in zip(y,y0)]
            new_turns, old_turns = turning_points(y), turning_points(y0)
            if key[0] == 'POTEN':
                observed_region_delta = max(d for d,b in zip(delta,y0) if b <= 45000)
                ok = (max(y) < 1e8 and min(y) >= min(y0)-100 and new_turns <= old_turns
                      and observed_region_delta <= 500 and max(delta) <= 5000)
            else:
                observed_region_delta = None
                # Morphing functions are dimensionless factors; compare to the
                # original shape rather than confusing them with physical SOCs.
                limit = max(.25, .5*max(map(abs,y0)))
                ok = max(delta) <= limit and new_turns <= old_turns+2
                if key[0] == 'BOBROT':
                    ok = ok and max(map(abs,y)) < .05
                # Fine-structure functions carry cm-1, unlike SOC/L+ morphing
                # factors. Small/sign-changing corrections are not invalid just
                # because their relative change from a near-zero seed is large.
                fine_limits = {'SPINSPIN':5., 'SPINROT':.5, 'LAMBDAOPQ':5.,
                               'LAMBDAP2Q':1., 'LAMBDAQ':1.}
                if key[0] in fine_limits:
                    slope = max(abs(b-a)/(radius[1]-radius[0]) for a,b in zip(y,y[1:]))
                    ok = (max(map(abs,y)) <= fine_limits[key[0]] and slope <= 20*fine_limits[key[0]]
                          and new_turns <= old_turns+2)
            result[':'.join(map(str,key))] = dict(ok=ok, kind=field.kind, name=field.name,
                max_change=max(delta), well_region_max_change_cm=observed_region_delta,
                turns=new_turns, baseline_turns=old_turns, minimum=min(y), maximum=max(y))
        except (ValueError, OverflowError) as exc:
            result[':'.join(map(str,key))] = dict(ok=False, error=str(exc))
    return dict(ok=all(r['ok'] for r in result.values()), range_angstrom=[.7,6], points=2001, fields=result)
