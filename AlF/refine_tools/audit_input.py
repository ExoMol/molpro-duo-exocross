"""Prepare reversible AlF input repairs; retain every source observation in an audit."""
from pathlib import Path
import collections
import hashlib
import json
import re
import argparse
import coupled

ROOT = Path(__file__).resolve().parent
SOURCE = Path('C:/sergei/programs/duo-inputs/AlF/27AlF_X1E+_A1Pi_EMO_02.inp')

def main():
    global ROOT, SOURCE
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input',type=Path,nargs='?',default=SOURCE)
    ap.add_argument('--output',type=Path,default=ROOT)
    args=ap.parse_args(); SOURCE=args.input.resolve(); ROOT=args.output.resolve()
    ROOT.mkdir(parents=True,exist_ok=True)
    original = SOURCE.read_text(encoding='cp1252')
    text = original.replace('\xa0', ' ')
    # No ab initio Lx grid exists in the supplied file. Interpret its existing
    # constant-one analytic function directly, and record this model assumption.
    text = re.sub(r'(?im)^morphing[ \t]*\n', '', text)
    excluded = []
    rows = coupled.observations(text)
    lines = []
    for n, line in enumerate(text.splitlines(), 1):
        p = line.split()
        if len(p) == 10 and re.fullmatch(r'-\d+', p[0]) and p[1] in ('+', '-'):
            excluded.append(dict(line=n, text=line, reason='negative J and energy; already skipped by native JLIST'))
            continue
        lines.append(line)
    text = '\n'.join(lines) + '\n'
    (ROOT/'runnable_original.inp').write_text(text)
    # A singlet Pi level has |Omega| = |Lambda| = 1. Preserve signed Lambda.
    clean = coupled.uncomment(text)
    match = re.search(r'(?ims)^energies[^\n]*\n(.*?)^end\s*$', clean)
    edits = []
    count = 0
    for m in re.finditer(r'(?m)^([^\n]+)$', match[1]):
        tokens = list(re.finditer(r'\S+', m[1]))
        if len(tokens) != 10 or tokens[4][0] != '2':
            continue
        omega = tokens[8]
        if omega[0] != '0':
            raise ValueError('Unexpected source Pi Omega')
        start = match.start(1) + m.start() + omega.start()
        edits.append((start, start + len(omega[0]), '-1'))
        count += 1
    for start, stop, value in reversed(edits):
        text = text[:start] + value + text[stop:]
    text = re.sub(r'(?im)^Jlist[^\n]*', 'Jlist 0 - 98', text)
    (ROOT/'prepared.inp').write_text(text)
    active = [r for r in rows if r['J'] >= 0]
    groups = collections.defaultdict(list)
    for row in active:
        key = tuple(row[k] for k in ('state', 'v', 'J', 'parity'))
        groups[key].append(row)
    audit = dict(source=str(SOURCE), source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                 nbsp_replaced=original.count('\xa0'), excluded_source_markers=excluded,
                 source_rows=len(rows), retained_rows=len(active),
                 positive_weight_rows=sum(r['input_weight'] > 0 for r in active),
                 omega_repairs=count, jlist_extended_from=92, jlist_extended_to=98,
                 newly_included_J_gt_92=sum(r['J'] > 92 for r in active),
                 Lx_assumption='No ab initio grid supplied: remove MORPHING; retain fixed analytic Lx=1',
                 duplicate_physical_levels=[dict(key=k, rows=v) for k,v in groups.items() if len(v)>1])
    match = re.search(r'(?ims)^energies[^\n]*\n(.*?)^end\s*$', text)
    merged, indices, duplicate_audit = [], {}, []
    for line in text[match.start(1):match.end(1)].splitlines():
        p = line.split()
        if not p:
            continue
        key = tuple(p[k] for k in (0, 1, 4, 5, 6, 7, 8))
        if key not in indices:
            indices[key] = len(merged)
            merged.append([p, [line]])
        else:
            old, sources = merged[indices[key]]
            w1, w2 = float(old[9]), float(p[9])
            old[3] = f'{(w1*float(old[3])+w2*float(p[3]))/(w1+w2):.10f}'
            old[9] = f'{w1+w2:g}'
            sources.append(line)
    for p, sources in merged:
        if len(sources)>1:
            duplicate_audit.append(dict(source_rows=sources, working_row=' '.join(p)))
    text = text[:match.start(1)] + '\n'.join(' '.join(p) for p,_ in merged) + '\n' + text[match.end(1):]
    (ROOT/'unique_levels.inp').write_text(text)
    audit['duplicate_weighted_means'] = duplicate_audit
    audit['duplicate_policy'] = 'Same physical QNs: retain summed weight and weighted mean; all raw records retained. Ordinary squared loss differs only by a fixed within-group constant. Parities are not guessed.'
    (ROOT/'input_audit.json').write_text(json.dumps(audit, indent=2)+'\n')
    print(json.dumps(audit, indent=2))

if __name__ == '__main__':
    main()
