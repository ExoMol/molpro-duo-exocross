"""Run a C2 Duo input in an isolated calculation directory."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('input', type=Path)
    ap.add_argument('output', type=Path)
    ap.add_argument('--duo', type=Path, default=Path(os.environ.get('DUO_EXE', 'duo')))
    ap.add_argument('--iterations', type=int, default=0)
    ap.add_argument('--jmax', type=int)
    ap.add_argument('--jlist', nargs='+', type=int)
    ap.add_argument('--scale', type=float)
    ap.add_argument('--robust', type=float)
    ap.add_argument('--lock', type=float)
    ap.add_argument('--stop-after-checkpoint', type=int)
    ap.add_argument('--precise', action='store_true')
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--timeout', type=float, default=7200)
    a = ap.parse_args()
    a.duo = Path(shutil.which(str(a.duo)) or a.duo)
    if not a.duo.is_file():
        ap.error('Native Duo executable not found; supply --duo or DUO_EXE')
    source = a.input.resolve()
    text = source.read_text()
    start = re.search(r'(?im)^\s*FITTING\s*$', text).start()
    before, fit = text[:start], text[start:]
    fit = re.sub(r'(?im)^\s*itmax\s+[^\n]*', f'itmax {a.iterations}', fit, count=1)
    fit = re.sub(r'(?im)^\s*output\s+[^\n]*', 'output fit', fit, count=1)
    if a.jmax is not None:
        fit = re.sub(r'(?im)^\s*JLIST\s+[^\n]*', f'JLIST 0 - {a.jmax}', fit, count=1)
    if a.jlist is not None:
        fit = re.sub(r'(?im)^\s*JLIST\s+[^\n]*', 'JLIST ' + ' '.join(map(str, a.jlist)), fit, count=1)
    if a.robust is not None:
        fit = re.sub(r'(?im)^\s*robust\s+[^\n]*', f'robust {a.robust}', fit, count=1)
    if a.lock is not None:
        fit = re.sub(r'(?im)^\s*lock\s+[^\n]*', f'lock {a.lock}', fit, count=1)
    if a.scale is not None:
        if re.search(r'(?im)^\s*fit_scale\b', fit):
            fit = re.sub(r'(?im)^\s*fit_scale\s+[^\n]*', f'fit_scale {a.scale}', fit, count=1)
        else:
            fit = fit.replace('output fit', f'fit_scale {a.scale}\noutput fit', 1)
    if a.precise:
        if a.iterations != 0 or a.jmax is None:
            ap.error('--precise requires a zero-iteration calculation and --jmax')
        if not re.search(r'(?im)^PRINT_ROVIBRONIC_ENERGIES_TO_FILE\b', before):
            before = 'PRINT_ROVIBRONIC_ENERGIES_TO_FILE\n' + before
        before = re.sub(r'(?im)^jrot\s+[^\n]*', f'jrot 0 - {a.jmax}', before, count=1)
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    inp = out / 'input.inp'
    inp.write_text(before + fit)
    env = os.environ.copy()
    # Optional local Intel runtime paths. Other builds use the caller's PATH.
    runtimes = [Path('C:/Program Files (x86)/Intel/oneAPI/compiler/2025.3/bin'),
                Path('C:/Program Files (x86)/Intel/oneAPI/mkl/2025.3/bin')]
    env['PATH'] = os.pathsep.join([str(p) for p in runtimes if p.is_dir()] + [env['PATH']])
    for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
        env[key] = str(a.threads)
    manifest = {'source': str(source), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                'input_sha256': hashlib.sha256(inp.read_bytes()).hexdigest(), 'duo': str(a.duo.resolve()),
                'duo_sha256': hashlib.sha256(a.duo.read_bytes()).hexdigest(), 'threads': a.threads}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    begin = time.monotonic()
    with inp.open('rb') as stdin, (out / 'duo.out').open('wb') as stdout:
        proc = subprocess.Popen([str(a.duo.resolve())], stdin=stdin, stdout=stdout, stderr=subprocess.STDOUT,
                                cwd=out, env=env)
        manifest['pid'] = proc.pid
        (out / 'manifest.json').write_text(json.dumps(manifest, indent=2))
        stopped = False
        while proc.poll() is None:
            time.sleep(1)
            if time.monotonic()-begin > a.timeout:
                proc.kill()
                proc.wait()
                manifest.update(returncode=proc.returncode, elapsed_seconds=time.monotonic()-begin,
                                status='timeout', complete_checkpoint=False,
                                controlled_checkpoint_stop=False)
                (out / 'manifest.json').write_text(json.dumps(manifest, indent=2))
                raise TimeoutError(f'Duo exceeded {a.timeout} seconds; outputs retained in {out}')
            if a.stop_after_checkpoint:
                log = (out/'duo.out').read_text(errors='replace')
                count = len(re.findall(r'(?m)^-->\|', log))
                if count >= a.stop_after_checkpoint:
                    proc.kill()
                    proc.wait()
                    stopped = True
                    break
    manifest.update(returncode=proc.returncode, elapsed_seconds=time.monotonic()-begin,
                    controlled_checkpoint_stop=stopped)
    log = (out / 'duo.out').read_text(errors='replace')
    manifest['complete_checkpoint'] = bool(re.search(r'(?m)^-->\|', log))
    manifest['native_diagnostics_present'] = all((out / name).is_file() for name in ('fit.en', 'fit.pot'))
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest, indent=2), flush=True)
    print(log[-1500:], flush=True)
    if not manifest['complete_checkpoint'] or not manifest['native_diagnostics_present']:
        raise SystemExit('Duo did not finish a valid calculation; Fortran STOP may return status zero.')
    raise SystemExit(0 if stopped else proc.returncode)

if __name__ == '__main__':
    main()
