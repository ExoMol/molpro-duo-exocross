"""Audited outer loop for coupled-state fits performed by native Duo.

Every proposed parameter set is evaluated in a fresh zero-iteration calculation.
Acceptance uses the original, fixed observation weights, not the robust weights
or Duo's threshold-filtered summary. This does not certify QN assignments.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import shutil
import os
import math
import csv
import re
deps = Path(__file__).resolve().parent.parent/'.pecfit_stage/plot_deps'
if deps.exists(): sys.path.insert(0,str(deps))
import coupled
import barrier

def objective(directory, stats, a):
    if a.criterion=='weighted_rms': return stats['weighted_rms_cm']
    with (directory/'residuals.csv').open() as stream: rows=list(csv.DictReader(stream))
    total=sum(float(r['input_weight']) for r in rows)
    return sum(math.log1p(.001*float(r['residual_cm'])**2*float(r['input_weight'])/(total*a.robust**2)) for r in rows)

def irls_weights(rows, robust):
    if robust<=0: raise ValueError('Positive robust scale required')
    total=sum(float(r['input_weight']) for r in rows)
    weights=[1/(robust**2*total/float(r['input_weight'])+.001*float(r['residual_cm'])**2) if float(r['input_weight'])>0 else 0. for r in rows]
    norm=sum(weights)
    return [w/norm for w in weights]


def native(source, output, a, fitting=False):
    robust_override=a.robust
    if fitting and a.robust is not None and a.robust>0:
        # Duo updates robust weights AFTER its first parameter step. Restarting
        # it for one checkpoint would otherwise repeat ordinary least squares.
        # Build the exact Watson IRLS weights from the last fresh full-J result.
        text=source.read_text()
        with (source.parent/'residuals.csv').open() as stream: rows=list(csv.DictReader(stream))
        weights=irls_weights(rows,a.robust)
        match=re.search(r'(?ims)^energies[^\n]*\n(.*?)^end\s*$',text)
        records=[line.split() for line in match[1].splitlines() if line.strip()]
        if len(records)!=len(rows): raise ValueError('IRLS row coverage mismatch')
        for p,row,w in zip(records,rows,weights):
            if float(p[0])!=float(row['J']) or p[1]!=row['parity'] or abs(float(p[3])-float(row['observed_cm']))>1e-8:
                raise ValueError('IRLS row identity mismatch')
            p[9]=f'{w:.16E}'
        weighted=output.with_suffix('.weighted.inp')
        weighted.write_text(text[:match.start(1)]+'\n'.join(' '.join(p) for p in records)+'\n'+text[match.end(1):])
        output.with_suffix('.weights.json').write_text(json.dumps(dict(robust=a.robust,alpha=.001,source=str(source),weights=weights),indent=2))
        source=weighted
        robust_override=0.
    command = [sys.executable, '-B', str(Path(__file__).with_name('run_native.py')),
               str(source), str(output), '--duo', str(a.duo), '--threads', str(a.threads),
               '--timeout', str(a.timeout)]
    if fitting:
        command += ['--iterations', str(a.proposal_steps), '--scale', str(a.scale), '--stop-after-checkpoint', str(a.proposal_steps)]
        command += ['--jlist', *map(str,a.fit_jlist)] if a.fit_jlist else ['--jmax', str(a.jmax)]
    else:
        command += ['--jmax', str(a.jmax), '--precise']
    if robust_override is not None:
        command += ['--robust', str(robust_override)]
    if a.lock is not None:
        command += ['--lock', str(a.lock)]
    result=subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        print(result.stdout,flush=True)
        raise RuntimeError(f'Native calculation failed: {output}')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input', type=Path)
    ap.add_argument('--duo', required=True, type=Path)
    ap.add_argument('--output', required=True, type=Path)
    ap.add_argument('--jmax', required=True, type=int)
    ap.add_argument('--fit-jlist', nargs='+', type=int)
    ap.add_argument('--rounds', type=int, default=3)
    ap.add_argument('--proposal-steps',type=int,default=1)
    ap.add_argument('--scale', type=float, default=.1)
    ap.add_argument('--robust', type=float)
    ap.add_argument('--criterion',choices=['weighted_rms','cauchy'],default='weighted_rms')
    ap.add_argument('--lock', type=float)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--timeout', type=float, default=7200)
    a = ap.parse_args()
    if a.criterion=='cauchy' and (a.robust is None or a.robust<=0):
        ap.error('Cauchy acceptance requires positive --robust')
    if a.rounds < 1 or a.proposal_steps<1 or not 0 < a.scale <= 1 or a.jmax < 0:
        ap.error('Require positive rounds, 0 < scale <= 1 and nonnegative Jmax')
    if a.fit_jlist and (min(a.fit_jlist) != 0 or max(a.fit_jlist) > a.jmax):
        ap.error('Fit J list must include zero and stay within the validation range')
    a.output.mkdir(parents=True, exist_ok=False)
    baseline_text = a.input.read_text()
    if not barrier.audit(baseline_text)['ok']:
        raise ValueError('Starting curves fail the AlF shape constraints')
    manifest = {k: str(v) if isinstance(v, Path) else v for k,v in vars(a).items()}
    manifest['source_sha256'] = hashlib.sha256(a.input.read_bytes()).hexdigest()
    manifest['code_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in Path(__file__).parent.glob('*.py')}
    (a.output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    snapshot = a.output/'code_snapshot'
    snapshot.mkdir()
    for path in Path(__file__).parent.glob('*.py'):
        shutil.copy2(path, snapshot/path.name)
    native(a.input, a.output/'baseline', a)
    baseline = coupled.save_assessment(a.output/'baseline', a.jmax)
    _, baseline_rows = coupled.assess((a.output/'baseline/input.inp').read_text(),
                                     (a.output/'baseline/duo.out').read_text(), a.jmax)
    baseline_unmatched = {i for i,r in enumerate(baseline_rows) if r['mark']}
    current = a.output/'baseline/input.inp'
    best_score = objective(a.output/'baseline',baseline,a)
    best_rms = baseline['weighted_rms_cm']
    history = []
    for iteration in range(1,a.rounds+1):
        proposal_dir = a.output/f'round_{iteration:02d}_proposal'
        native(current, proposal_dir, a, fitting=True)
        checkpoint = coupled.checkpoints((proposal_dir/'duo.out').read_text())[-1]
        before = current.read_text()
        accepted = None
        for fraction in (1., .5, .25, .125, .0625):
            text = coupled.transfer(before, checkpoint, fraction)
            shape = barrier.audit(text)
            label = f'round_{iteration:02d}_fraction_{fraction:g}'
            candidate = a.output/f'{label}.inp'
            candidate.write_text(text)
            (a.output/f'{label}.shape.json').write_text(json.dumps(shape, indent=2)+'\n')
            entry = dict(round=iteration, fraction=fraction, shape_ok=shape['ok'])
            if shape['ok']:
                native(candidate, a.output/label, a)
                score = coupled.save_assessment(a.output/label, a.jmax)
                _, candidate_rows = coupled.assess((a.output/label/'input.inp').read_text(),
                                                   (a.output/label/'duo.out').read_text(), a.jmax)
                unmatched = {i for i,r in enumerate(candidate_rows) if r['mark']}
                entry.update(weighted_rms_cm=score['weighted_rms_cm'], rms_cm=score['rms_cm'],
                             marked_count=score['marked_count'])
                # Retain the same full observation set and reject new unmatched
                # observations. Approximate Omega labels are audited separately.
                candidate_score=objective(a.output/label,score,a)
                entry['objective']=candidate_score
                if (candidate_score < best_score*(1-1e-7)
                        and unmatched <= baseline_unmatched):
                    best_score = candidate_score
                    best_rms = score['weighted_rms_cm']
                    current = a.output/label/'input.inp'
                    accepted = label
                    entry['accepted'] = True
            history.append(entry)
            print(json.dumps(entry),flush=True)
            (a.output/'history.json').write_text(json.dumps(history, indent=2)+'\n')
            if accepted:
                break
        if not accepted:
            break
    (a.output/'best.inp').write_text(current.read_text())
    result = dict(status='needs_review', baseline_weighted_rms_cm=baseline['weighted_rms_cm'],
                  best_weighted_rms_cm=best_rms, criterion=a.criterion, best_objective=best_score,
                  selected_calculation=str(current.parent),
                  validation='Fresh full-J evaluation completed; basis convergence and mixed-state assignments require separate review.')
    (a.output/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
