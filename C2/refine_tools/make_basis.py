import argparse
from pathlib import Path
import re

ap = argparse.ArgumentParser()
ap.add_argument('input', type=Path)
ap.add_argument('output', type=Path)
ap.add_argument('--npoints', type=int, default=501)
ap.add_argument('--rmin', type=float, default=.8)
ap.add_argument('--rmax', type=float, default=4.5)
ap.add_argument('--vadd', type=int, default=10)
a = ap.parse_args()
text = a.input.read_text()
text = re.sub(r'(?im)^\s*npoints\s+[^\n]*', f'  npoints {a.npoints}', text, count=1)
text = re.sub(r'(?im)^\s*range\s+[^\n]*', f'  range {a.rmin}, {a.rmax}', text, count=1)
m = re.search(r'(?im)^[ \t]*vmax[ \t]+([\d \t]+)', text)
text = text[:m.start()] + '  vmax ' + ' '.join(str(int(x)+a.vadd) for x in m[1].split()) + text[m.end():]
with a.output.open('x') as out:
    out.write(text)
