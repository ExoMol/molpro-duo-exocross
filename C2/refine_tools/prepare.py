"""Prepare duplicate observations explicitly and record every altered row."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import coupled


def prepare(text, *, merge_duplicates=False, tolerance=.001, lock=None, jmax=None):
    clean = coupled.uncomment(text)
    match = re.search(r'(?ims)^\s*FITTING\s*$.*?^\s*energies\s*$\n(.*?)^\s*end\s*$', clean)
    if match is None:
        raise ValueError('Expected a FITTING energies block')
    start, stop = match.span(1)
    groups, records, merged = {}, [], []
    for line in text[start:stop].splitlines():
        p = coupled.uncomment(line).split()
        if not p:
            records.append(dict(source=[line], p=None))
            continue
        key = tuple(p[k] for k in (0,1,2,4,5,6,7,8))
        energy, weight = coupled.number(p[3]), coupled.number(p[9])
        old = groups.get(key) if merge_duplicates else None
        if old:
            if abs(float(old['p'][3])-energy) > tolerance:
                raise ValueError(f'Inconsistent duplicate energy for {key}')
            old['source'].append(line)
            old['weighted_energy'] += energy*weight
            old['weight'] += weight
        else:
            item = dict(p=p, source=[line], weighted_energy=energy*weight, weight=weight)
            groups[key] = item
            records.append(item)
    output = []
    for item in records:
        if len(item['source']) > 1:
            if item['weight'] <= 0:
                raise ValueError('Cannot merge zero-weight duplicate observations')
            p = item['p'][:]
            p[3] = f"{item['weighted_energy']/item['weight']:.10f}"
            p[9] = f"{item['weight']:.10g}"
            merged.append(dict(original_rows=item['source'], working_row=' '.join(p)))
            output.append(' '.join(p))
        else:
            output.append(item['source'][0])
    result = text[:start]+'\n'.join(output)+'\n'+text[stop:]
    if lock is not None:
        result = re.sub(r'(?im)^lock\s+[^\n]*', f'lock {lock}', result, count=1)
    if jmax is not None:
        result = re.sub(r'(?im)^JLIST\s+[^\n]*', f'JLIST 0 - {jmax}', result, count=1)
    return result, dict(duplicate_groups=merged, lock=lock, jmax=jmax,
        note='Duplicate weighted averaging preserves ordinary weighted least squares up to a constant, but does not exactly preserve a robust loss. Original records are retained here. Negative lock relaxes approximate Lambda/Sigma/Omega labels while retaining state/v within the given energy window.')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input', type=Path)
    ap.add_argument('output', type=Path)
    ap.add_argument('--merge-duplicates', action='store_true')
    ap.add_argument('--lock', type=float)
    ap.add_argument('--jmax', type=int)
    a = ap.parse_args()
    text, audit = prepare(a.input.read_text(), merge_duplicates=a.merge_duplicates, lock=a.lock, jmax=a.jmax)
    audit['source'] = str(a.input.resolve())
    audit['source_sha256'] = hashlib.sha256(a.input.read_bytes()).hexdigest()
    with a.output.open('x') as stream:
        stream.write(text)
    with a.output.with_suffix('.preparation.json').open('x') as stream:
        json.dump(audit, stream, indent=2)
    print(f"Prepared {a.output}; merged {len(audit['duplicate_groups'])} duplicate groups.")


if __name__ == '__main__':
    main()
