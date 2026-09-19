"""Lossless parameter transfer and fixed-data diagnostics for coupled-state Duo fits.

Parameters are identified by object, states and position, never by name alone:
compound functions such as TWO_COUPLED_EMOS contain repeated parameter names.
The native Duo Hamiltonian and assignment algorithm remain authoritative.
"""
from dataclasses import dataclass
import csv
import json
import math
from pathlib import Path
import re


def number(s):
    x = float(s.replace('D', 'E').replace('d', 'e'))
    if not math.isfinite(x):
        raise ValueError(f'Non-finite number: {s}')
    return x


def canonical(s):
    s = s.upper().replace('-', '').replace('_', '')
    return {'SPINORBITX': 'SPINORBIT', 'LX': 'L+', 'BOBROT': 'BOBROT'}.get(s, s)


CLASSES = {'POTEN', 'SPINORBIT', 'L+', 'SPINSPIN', 'SPINROT', 'BOBROT',
           'LAMBDAOPQ', 'LAMBDAP2Q', 'LAMBDAQ', 'L2', 'DIABATIC'}


def uncomment(text):
    """Blank Duo's line-local comments, preserving offsets and quoted names."""
    chars = list(text)
    depth = 0
    quote = None
    bang = False
    for i, c in enumerate(text):
        if c == '\n':
            bang = False
            depth = 0
            quote = None
            continue
        if bang:
            chars[i] = ' '
        elif depth:
            if c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
            chars[i] = ' '
        elif quote:
            if c == quote:
                quote = None
        elif c in '\"\'':
            quote = c
        elif c == '(':
            depth = 1
            chars[i] = ' '
        elif c == '!':
            bang = True
            chars[i] = ' '
    return ''.join(chars)


@dataclass
class Field:
    key: tuple
    kind: str
    name: str
    labels: list
    values: list
    spans: list
    fitted: list


def fields(text):
    clean = uncomment(text)
    result = {}
    header = re.compile(r'(?im)^([\w+\-]+)[ \t]+(\d+)(?:[ \t]+(\d+))?[^\n]*$')
    for match in header.finditer(clean):
        cls = canonical(match[1])
        if cls not in CLASSES:
            continue
        end = re.search(r'(?im)^\s*end\s*$', clean[match.end():])
        if end is None:
            raise ValueError(f'Missing end: {match[0]}')
        stop = match.end() + end.start()
        block = clean[match.end():stop]
        val = re.search(r'(?im)^\s*values\s*$', block)
        if val is None:
            continue
        kinds = re.findall(r'(?im)^\s*type\s+(\S+)', block[:val.start()])
        name = re.search(r'(?im)^\s*name\s+([^\n]+)', block[:val.start()])
        if not kinds:
            raise ValueError(f'Missing type: {match[0]}')
        key = cls, int(match[2]), int(match[3] or match[2])
        if key in result:
            raise ValueError(f'Duplicate field: {key}')
        start = match.end() + val.end()
        labels, values, spans, fitted = [], [], [], []
        for row in re.finditer(r'(?m)^\s*(\w+)\s+([-+\d.eEdD]+)([^\n]*)$', clean[start:stop]):
            labels.append(row[1].upper())
            values.append(number(row[2]))
            spans.append((start + row.start(2), start + row.end(2)))
            fitted.append(bool(re.search(r'\bfit\b', row[3], re.I)))
        result[key] = Field(key, kinds[-1].upper(), name[1].strip().strip('\"') if name else '',
                            labels, values, spans, fitted)
    return result


def checkpoints(output):
    """Only complete full-precision parameter blocks before a summary table."""
    result = []
    for match in re.finditer(r'(?m)^Parameters:\s*$', output):
        tail = output[match.end():]
        end = re.search(r'(?m)^Fitted parameters \(rounded\):|^   [-]{20,}', tail)
        if end is None:
            continue
        f = fields(tail[:end.start()])
        if f:
            result.append(f)
    return result


def transfer(template, checkpoint, fraction=1.0):
    src = fields(template)
    if src.keys() != checkpoint.keys():
        raise ValueError(f'Field coverage mismatch: {src.keys() ^ checkpoint.keys()}')
    edits = []
    for key, field in src.items():
        target = checkpoint[key]
        if field.kind != target.kind or field.labels != target.labels:
            raise ValueError(f'Parameter layout changed: {key}')
        for old, new, span in zip(field.values, target.values, field.spans):
            if old != new:
                edits.append((*span, f'{old + fraction * (new-old):.14E}'))
    for start, end, value in sorted(edits, reverse=True):
        template = template[:start] + value + template[end:]
    return template


def observations(text, jmax=None):
    clean = uncomment(text)
    block = re.search(r'(?ims)^\s*FITTING\s*$.*?^\s*energies\s*$\n(.*?)^\s*end\s*$', clean)
    if block is None:
        raise ValueError('Expected one FITTING energies block')
    result = []
    for line in block[1].splitlines():
        p = line.split()
        if not p:
            continue
        if len(p) < 10:
            raise ValueError(f'Malformed observation: {line}')
        j = number(p[0])
        if jmax is not None and j > jmax:
            continue
        q = [number(x) for x in p[4:9]]
        if number(p[9]) < 0:
            raise ValueError('Negative observation weights are not supported')
        result.append(dict(J=j, parity=p[1], original_rank=int(p[2]), observed_cm=number(p[3]),
                           state=int(q[0]), v=int(q[1]), Lambda=q[2], Sigma=q[3], Omega=q[4],
                           input_weight=number(p[9])))
    return result


ROW = re.compile(r'^\s*(\d+)\s+(\d+)\s+([\d.]+)\s+([+-])\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([\d.Ee+-]+)\s*\(([^()]*)\)\(([^()]*)\)(.*)$')


def energy_tables(output):
    tables, rows = [], []
    for line in output.splitlines():
        m = ROW.match(line)
        if not m:
            continue
        idx = int(m[1])
        if idx == 1 and rows:
            tables.append(rows)
            rows = []
        rows.append(dict(index=idx, rank=int(m[2]), J=number(m[3]), parity=m[4],
                         observed=number(m[5]), calculated=number(m[6]), residual=number(m[7]),
                         duo_weight=number(m[8]), calc_q=[number(x) for x in m[9].split()],
                         obs_q=[number(x) for x in m[10].split()], mark=m[11].strip()))
    if rows:
        tables.append(rows)
    return tables


def metrics(rows):
    active = [r for r in rows if r['input_weight'] > 0]
    if not active:
        return {'count': len(rows), 'positive_weight_count': 0}
    residuals = [r['residual_cm'] for r in active]
    return dict(count=len(rows), positive_weight_count=len(active),
                rms_cm=math.sqrt(sum(x*x for x in residuals)/len(active)),
                weighted_rms_cm=math.sqrt(sum(r['input_weight']*r['residual_cm']**2 for r in active)/sum(r['input_weight'] for r in active)),
                max_abs_cm=max(map(abs, residuals)),
                assignment_mismatches=sum(not r['assignment_matches'] for r in rows),
                marked_count=sum(bool(r['mark']) for r in rows),
                unmatched_positive_weight_count=sum(bool(r['mark']) for r in active),
                state_v_mismatches=sum(r['state'] != r['calculated_state'] or r['v'] != r['calculated_v'] for r in active),
                zero_duo_weight_count=sum(r['duo_weight'] == 0 for r in active))


def assess(text, output, jmax):
    all_obs = observations(text, jmax)
    states = {key[1] for key in fields(text) if key[0] == 'POTEN'}
    omitted = [o for o in all_obs if o['state'] not in states]
    if any(o['input_weight'] > 0 for o in omitted):
        raise ValueError('Positive-weight observations reference states absent from the model')
    obs = [o for o in all_obs if o['state'] in states]
    tables = energy_tables(output)
    if not tables or len(tables[-1]) != len(obs):
        raise ValueError(f'Incomplete observation coverage: {len(tables[-1]) if tables else 0}/{len(obs)}')
    result = []
    for i, (o, c) in enumerate(zip(obs, tables[-1]), 1):
        q = [o[k] for k in ('state', 'v', 'Lambda', 'Sigma', 'Omega')]
        if c['index'] != i or c['J'] != o['J'] or c['parity'] != o['parity'] or c['obs_q'] != q or abs(c['observed']-o['observed_cm']) > 0.000051:
            raise ValueError(f'Observation identity changed at row {i}: {o}, {c}')
        cq = c['calc_q']
        same = len(cq) == 5 and cq[:2] == q[:2] and list(map(abs, cq[2:])) == list(map(abs, q[2:]))
        result.append(dict(o, calculated_cm=c['calculated'], residual_cm=o['observed_cm']-c['calculated'],
                           matched_rank=c['rank'], calculated_state=int(cq[0]), calculated_v=int(cq[1]),
                           calculated_Lambda=cq[2], calculated_Sigma=cq[3], calculated_Omega=cq[4],
                           assignment_matches=same, mark=c['mark'], duo_weight=c['duo_weight']))
    stats = metrics(result)
    stats['inactive_unknown_state_rows'] = omitted
    stats['by_state'] = {str(s): metrics([r for r in result if r['state'] == s]) for s in sorted({r['state'] for r in result})}
    return stats, result


def save_assessment(directory, jmax):
    directory = Path(directory)
    stats, rows = assess((directory/'input.inp').read_text(), (directory/'duo.out').read_text(errors='replace'), jmax)
    precise = directory/'rovibronic_energies.dat'
    if precise.exists():
        energies = {}
        for line in precise.read_text().splitlines():
            if '||' not in line:
                continue
            p = line.split('||')[0].split()
            if len(p) != 10:
                raise ValueError('Malformed precise energy row')
            key = number(p[0]), p[9], int(p[1])
            if key in energies:
                raise ValueError(f'Duplicate precise energy: {key}')
            energies[key] = number(p[2])
        for row in rows:
            value = energies[(row['J'],row['parity'],row['matched_rank'])]
            if abs(value-row['calculated_cm']) > 0.000051:
                raise ValueError('Precise energies disagree with native fitting output')
            row['calculated_cm'] = value
            row['residual_cm'] = row['observed_cm']-value
        inactive = stats['inactive_unknown_state_rows']
        stats = metrics(rows)
        stats['inactive_unknown_state_rows'] = inactive
        stats['by_state'] = {str(s): metrics([r for r in rows if r['state'] == s]) for s in sorted({r['state'] for r in rows})}
        stats['energy_precision'] = 'native 12 decimal places'
    cutoff = re.search(r'(?im)^thresh_obs-calc\s+([\d.eEdD+-]+)', (directory/'input.inp').read_text())
    if cutoff:
        stats['input_residual_cutoff_cm'] = number(cutoff[1])
        stats['positive_weight_rows_above_cutoff'] = sum(abs(r['residual_cm']) > number(cutoff[1]) for r in rows if r['input_weight'] > 0)
    (directory/'assessment.json').write_text(json.dumps(stats, indent=2)+'\n')
    with (directory/'residuals.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return stats


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('directory', type=Path)
    ap.add_argument('--jmax', type=int, required=True)
    a = ap.parse_args()
    print(json.dumps(save_assessment(a.directory, a.jmax), indent=2))
