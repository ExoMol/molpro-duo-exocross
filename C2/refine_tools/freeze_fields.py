import argparse
from pathlib import Path
import re
import coupled

ap = argparse.ArgumentParser()
ap.add_argument('input', type=Path)
ap.add_argument('output', type=Path)
ap.add_argument('--classes', nargs='+', default=['SPINORBIT','L+'])
ap.add_argument('--add-coefficient', nargs=2, action='append', metavar=('STATE','LABEL'))
a = ap.parse_args()
text = a.input.read_text()
edits = []
for key, field in coupled.fields(text).items():
    for label, (_, stop), fitted in zip(field.labels, field.spans, field.fitted):
        end = text.find('\n', stop)
        suffix = text[stop:end]
        if key[0] in a.classes and fitted:
            suffix = re.sub(r'\bfit\b', '   ', suffix, count=1, flags=re.I)
            edits.append((stop, end, suffix))
        if key[0] == 'POTEN' and [str(key[1]),label] in (a.add_coefficient or []) and not fitted:
            edits.append((stop, stop, ' fit'))
for start, end, value in sorted(edits, reverse=True):
    text = text[:start]+value+text[end:]
with a.output.open('x') as stream:
    stream.write(text)
print('Active parameters:', sum(sum(f.fitted) for f in coupled.fields(text).values()))
