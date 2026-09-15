# CO automatic fitting validation — 2026-09-11

The complete Python workflow was exercised against the supplied CO data using
the local Windows IFX Duo build. All 27 tests pass, including real Duo
execution, monitored stopping, timeout handling, parameter extraction,
assignment coverage, recovery and cache integrity. A CLI `--resume` replay used
only verified cached calculations and left the final model, residuals and basis
comparison byte-for-byte unchanged. It returned exit status 2.

| Check | Result |
|---|---:|
| Experimental levels | 2293 |
| J range | 0–123 |
| Unweighted RMS | 0.14429389 cm-1 |
| Original-weight RMS | 0.000334646 cm-1 |
| Maximum absolute residual | 0.575390 cm-1 |
| Required unweighted RMS | 0.1 cm-1 |
| Final right-hand EMO order NR | 8 |
| Maximum basis/grid energy change | 0.0001000 cm-1 |
| Required basis tolerance | 0.001 cm-1 |
| Shape check on 0.8–3.0 Angstrom | passed |
| Outcome | needs_review |

The higher-order trial worsened unweighted RMS and was rejected. The retained
model passes assignment, shape and basis checks, but does **not** reach the
requested 0.1 cm-1 unweighted accuracy. No levels were hidden by a residual
cutoff in final validation. The pipeline explicitly reports this limitation;
its weighted RMS is not substituted for the requested unweighted target.
No experimental De was supplied, so DE remained fitted.

The final artifacts are in [pecfit_runs](pecfit_runs/):
[Duo input](pecfit_runs/final.inp), [report](pecfit_runs/report.md),
[residuals](pecfit_runs/final_residuals.csv),
[basis comparison](pecfit_runs/basis_convergence.csv), and
[resume log](pecfit_runs/cli_resume.log). Generated run directories are ignored
by Git. The reusable code and instructions are in [PEC_FITTING.md](../PEC_FITTING.md).

Largest residuals (observed minus calculated):

| J | v | Residual / cm-1 | Original weight |
|---:|---:|---:|---:|
| 115 | 7 | -0.575390 | 0.00983284 |
| 110 | 8 | -0.556560 | 0.0205529 |
| 118 | 6 | -0.544300 | 0.0180018 |
| 108 | 8 | -0.537140 | 0.0316522 |
| 111 | 7 | -0.534890 | 0.0209996 |
| 116 | 6 | -0.523590 | 0.0200481 |
| 115 | 6 | -0.520550 | 0.0196464 |
| 109 | 7 | -0.517250 | 0.0209864 |

Executable: `C:\sergei\programs\Duo\build\ifx\duo.exe`.
SHA-256: `3d45ab6e21673fd13a04d25f9ea9b20145a81a45d21d9601d2c2db49787c1e85`.
Full configuration and source/code fingerprints are in `pecfit_runs/manifest.json`.
During development, completed identical Duo calculations were migrated between
test directories; each reused calculation's input, executable and output hashes
were checked and Python acceptance decisions were replayed. The final polishing,
higher-order and basis checks were run with the delivered implementation.
