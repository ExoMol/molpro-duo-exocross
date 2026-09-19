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
import coupled
import curves


def native(source, output, a, fitting=False):
    command = [sys.executable, '-B', str(Path(__file__).with_name('run_native.py')),
               str(source), str(output), '--duo', str(a.duo), '--threads', str(a.threads),
               '--timeout', str(a.timeout)]
    if fitting:
        command += ['--iterations', '1', '--scale', str(a.scale), '--stop-after-checkpoint', '1']
        command += ['--jlist', *map(str,a.fit_jlist)] if a.fit_jlist else ['--jmax', str(a.jmax)]
    else:
        command += ['--jmax', str(a.jmax)]
    if a.robust is not None:
        command += ['--robust', str(a.robust)]
    if a.lock is not None:
        command += ['--lock', str(a.lock)]
    subprocess.run(command, check=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input', type=Path)
    ap.add_argument('--duo', required=True, type=Path)
    ap.add_argument('--output', required=True, type=Path)
    ap.add_argument('--jmax', required=True, type=int)
    ap.add_argument('--fit-jlist', nargs='+', type=int)
    ap.add_argument('--rounds', type=int, default=3)
    ap.add_argument('--scale', type=float, default=.1)
    ap.add_argument('--robust', type=float)
    ap.add_argument('--lock', type=float)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--timeout', type=float, default=7200)
    a = ap.parse_args()
    if a.rounds < 1 or not 0 < a.scale <= 1 or a.jmax < 0:
        ap.error('Require positive rounds, 0 < scale <= 1 and nonnegative Jmax')
    if a.fit_jlist and (min(a.fit_jlist) != 0 or max(a.fit_jlist) > a.jmax):
        ap.error('Fit J list must include zero and stay within the validation range')
    a.output.mkdir(parents=True, exist_ok=False)
    baseline_text = a.input.read_text()
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
    best_score = baseline['weighted_rms_cm']
    history = []
    for iteration in range(1,a.rounds+1):
        proposal_dir = a.output/f'round_{iteration:02d}_proposal'
        native(current, proposal_dir, a, fitting=True)
        checkpoint = coupled.checkpoints((proposal_dir/'duo.out').read_text())[-1]
        before = current.read_text()
        accepted = None
        for fraction in (1., .5, .25):
            text = coupled.transfer(before, checkpoint, fraction)
            shape = curves.audit(text, baseline_text)
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
                if (score['weighted_rms_cm'] < best_score*(1-1e-5)
                        and unmatched <= baseline_unmatched):
                    best_score = score['weighted_rms_cm']
                    current = a.output/label/'input.inp'
                    accepted = label
                    entry['accepted'] = True
            history.append(entry)
            (a.output/'history.json').write_text(json.dumps(history, indent=2)+'\n')
            if accepted:
                break
        if not accepted:
            break
    (a.output/'best.inp').write_text(current.read_text())
    result = dict(status='needs_review', baseline_weighted_rms_cm=baseline['weighted_rms_cm'],
                  best_weighted_rms_cm=best_score, selected_calculation=str(current.parent),
                  validation='Fresh full-J evaluation completed; basis convergence and mixed-state assignments require separate review.')
    (a.output/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
