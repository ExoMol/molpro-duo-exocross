# Coupled-state refinement with native Duo

These tools preserve the complete supplied Hamiltonian: multiple potentials,
spin-orbit and other couplings, morphing functions, linked parameters, symmetry,
ab initio constraints and experimental observations. Native Duo performs the
coupled rovibronic calculations and least-squares derivatives. The original
single-state `pecfit` commands remain unchanged.

Python 3.10+ is sufficient; supply a working native Duo executable. Run from
the repository root. Output directories must be new, so existing results cannot
be overwritten accidentally.

```powershell
python C2/refine_tools/prepare.py C2/source_input.inp C2/prepared_new.inp --merge-duplicates --lock -7 --jmax 87
python C2/refine_tools/refine.py C2/prepared_new.inp --duo C:/path/to/duo.exe --output C2/pecfit_runs/new_trial --jmax 87 --fit-jlist 0 1 2 5 10 15 20 25 30 35 40 50 60 70 80 87 --rounds 3 --scale 0.1
```

Preparation is explicit. Duplicate consolidation requires matching J, parity,
level number and all supplied quantum numbers, and energies within 0.001 cm-1.
It sums weights and uses their weighted mean energy, recording every original
row. This preserves ordinary weighted least squares up to a constant; robust
losses are not exactly invariant. Negative `lock` relaxes approximate
Lambda/Sigma/Omega labels while requiring the same electronic state and v within
the stated energy window. It is an assignment hypothesis to review, not proof
of an experimental reassignment. Do not apply it blindly to another molecule.

The outer loop obtains a full-precision native parameter checkpoint, then
evaluates each trial in a fresh calculation on **all** J values through `--jmax`.
Sparse `--fit-jlist` sampling reduces derivative cost; no sampled fit is accepted
without the full-range check. Repeated parameter names inside compound
functions are transferred by object, state pair and parameter position.
Input metadata and the observation table survive parameter transfer.

Acceptance requires improved RMS using the fixed supplied weights, curve checks,
complete observation coverage and no newly unmatched observations. Backtracking
tries fractions 1, 0.5 and 0.25. Robust weights and the native residual cutoff are
never used to hide observations in the reported RMS. Controlled checkpoint stops
are recorded and followed by fresh evaluations. Native Fortran `STOP` can return
exit status zero; missing checkpoints or diagnostics still count as failure.

`--robust` is optional and otherwise retains the input setting. The C2 experiments
also tested a separately labelled `--robust 0` polish, froze the main spin-orbit
and L+ morphing parameters, and activated c-state B3. See the refinement report
for the exact selected input rather than assuming the example command recreates
all exploratory stages.

```powershell
python C2/refine_tools/run_native.py C2/refined_model/C2_refined.inp C2/pecfit_runs/recheck --duo C:/path/to/duo.exe --jmax 87 --precise
python C2/refine_tools/coupled.py C2/pecfit_runs/recheck --jmax 87
python -m unittest discover -s C2/refine_tools/tests -v
```

`make_basis.py` creates a separate input with changed radial grid and larger
contractions for convergence checks. `check_native_curves.py` compares the
analytic PEC diagnostic functions with native Duo's rounded ab initio tables.
`freeze_fields.py` changes fit flags explicitly; it does not change curve values.

Current diagnostic scope is numeric electronic-state labels and energy fits
using EMO, TWO_COUPLED_EMOS, POLYNOM_DECAY_24 and BOBLEROY functions. Unknown forms
fail the curve audit. The shape limits were chosen for this C2 trial and must be
reviewed for another molecule: PEC topology cannot add extrema relative to the
seed; PEC changes are limited in the observed energy region and outer grid;
SOC/L+ morphing factors have conservative relative-change bounds; BOB magnitude
is below 0.05; fine-structure functions have dimensioned magnitude and slope
limits. Existing avoided crossings/double-well structure are retained.

The final status remains `needs_review`: the outer loop does not certify mixed
state assignments, complete basis convergence or a global optimum. All input,
executable and code hashes, proposals, native outputs and rejection diagnostics
are retained. New runs snapshot the Python tools. There is no automatic resume;
start a new output directory from a selected input.
