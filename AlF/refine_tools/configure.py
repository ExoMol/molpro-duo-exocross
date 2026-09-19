"""Edit parameter positions explicitly; compound blocks may repeat their names."""
from pathlib import Path
import argparse
import re
import coupled

def configure(text, fit=(), freeze=(), values=()):
    fields=coupled.fields(text); edits=[]
    requests={}
    for action,items in [('fit',fit),('freeze',freeze),('value',values)]:
        for item in items:
            k,*val=item.split('='); cls,s1,s2,index=k.split(':')
            key=(coupled.canonical(cls),int(s1),int(s2)); idx=int(index)
            f=fields[key]; start,stop=f.spans[idx]
            request=requests.setdefault((key,idx),dict(start=start,stop=stop,value=f.values[idx],fit=f.fitted[idx]))
            if action=='value': request['value']=float(val[0])
            else: request['fit']=action=='fit'
    for (key,i),req in requests.items():
        start,stop=req['start'],req['stop']
        end=text.find('\n',stop)
        suffix=text[stop:end]
        suffix=re.sub(r'\bfit\b','',suffix,flags=re.I)
        edits.append((start,end,f"{req['value']:.14E}"+(' fit' if req['fit'] else '')+suffix))
    for a,b,value in sorted(edits,reverse=True): text=text[:a]+value+text[b:]
    return text

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input',type=Path); ap.add_argument('output',type=Path)
    ap.add_argument('--fit',action='append',default=[])
    ap.add_argument('--freeze',action='append',default=[])
    ap.add_argument('--value',action='append',default=[])
    a=ap.parse_args()
    a.output.write_text(configure(a.input.read_text(),a.fit,a.freeze,a.value))
