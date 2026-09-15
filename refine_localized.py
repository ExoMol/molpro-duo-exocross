"""Optional refinement for single-state bound and localized quasibound levels."""
import argparse
import json
import os


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config')
    parser.add_argument('--duo',required=True)
    parser.add_argument('--output')
    parser.add_argument('--max-iterations',type=int,default=100)
    parser.add_argument('--threads',type=int,default=1,help='BLAS threads for the Python eigensolver')
    parser.add_argument('--evaluate-only',action='store_true')
    parser.add_argument('--robust',type=float,help='Duo/Watson robust scale; defaults to the configuration value; 0 disables reweighting')
    parser.add_argument('--fit-below',type=float,help='Optional explicit calibration cutoff in cm-1; every source level remains in final diagnostics')
    parser.add_argument('--optimizer',choices=('guarded','penalized'),default='guarded')
    parser.add_argument('--preserve-c6',action='store_true',help='Preserve the MLJ seed C6 when varying De or Re')
    args=parser.parse_args()
    if args.threads<1 or args.max_iterations<1:parser.error('Threads and iterations must be positive')
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]=str(args.threads)
    try:
        from pecfit.localized import execute
        from pecfit.model import FitError
    except ImportError as exc:
        parser.error(f'Install the optional numerical dependencies: pip install ".[localized]" ({exc})')
    try:
        result=execute(args.config,args.duo,output=args.output,max_iterations=args.max_iterations,
                       evaluate_only=args.evaluate_only,preserve_c6=args.preserve_c6,robust=args.robust,fit_below=args.fit_below,optimizer=args.optimizer)
    except (FitError,OSError,ValueError) as exc:
        print(f'Localized refinement failed: {exc}')
        return 1
    print(json.dumps({k:result[k] for k in ('status','experiment','bound_experiment','quasibound_count','basis_max_delta_cm')},indent=2))
    return 0 if result['status']=='success' else 2


if __name__=='__main__':raise SystemExit(main())
