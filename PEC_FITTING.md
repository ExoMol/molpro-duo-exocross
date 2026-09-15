# Automatic single-state potential fitting with Duo

`pecfit` implements the two-step procedure in `CO/readme.txt`: fit an extended
Morse oscillator (EMO) to the ab initio potential, then refine its parameters
against experimental or MARVEL term values. It uses Duo's own least-squares
fitter and nuclear-motion solver. Python handles input generation, parameter
transfer, staged J ranges, restarts, rejected steps, and final validation.

Scope: one electronic X1Sigma+ state, one isotopologue, one EMO or MLJ potential,
an optional rotational BOB `polynom_decay_24` function, and
experimental **term energies**, not transition frequencies. MOLPRO calculations,
MARVEL inversion, intensity fitting, and line-list generation are separate tasks.
The generated fitting inputs contain the potential, optional BOB, and energy data; the
original sample inputs and their dipole/intensity blocks are not edited.

## Run CO

Use Python 3.10 or newer. The core Duo workflow has no Python runtime dependencies. Run commands
from the repository root; installing the package is optional.

```powershell
python -m pecfit plan CO/pecfit.json
python -m pecfit run CO/pecfit.json --duo "C:/path/to/duo.exe"
```

On Linux/macOS, supply the native Duo executable instead:

```sh
python -m pecfit run CO/pecfit.json --duo /path/to/duo
```

Alternatively, put `duo` on PATH, set `DUO_EXE`, or set `duo_executable` in the
JSON configuration. Relative configuration paths are relative to that JSON
file. Relative command-line paths are relative to the current directory.

The Windows executable `C:/sergei/programs/Duo/executables/duo-win-220626.exe`
was used for development tests. The repository's `CO/duo-project.exe` returned
Windows status `0xC0000135` (missing runtime/DLL) on the development machine.
Use a working native build; executable binaries are not downloaded by the
pipeline. A DOS command prompt can launch Windows Duo executables, but a Linux
binary such as `j-duo.x` cannot run directly on Windows.
The full CO validation also uses the locally built `Duo/build/ifx/duo.exe`.
Intel-linked builds need the Intel oneAPI runtime environment initialised.

The supplied CO configuration uses the **301-point, already shifted cm-1 grid**
in `CO/sample_outputs/CO_duo_fit_01.inp`. It does not replace it with the sparse
26-point Hartree grid in `pec.txt`. This reproduces the input used in the
refinement instructions, including the prepared grid's extrapolated region.

Output goes to `CO/pecfit_runs/`. To restart after interruption:

```powershell
python -m pecfit run CO/pecfit.json --duo "C:/path/to/duo.exe" --resume
```

Resume deterministically replays stage decisions and reuses complete attempts
only when their input, executable, thread count, and diagnostic file hashes
match. Configuration, source files, and pipeline code must be unchanged.
Interrupted or damaged attempts are rerun. Use a new directory for changed
settings:

```powershell
python -m pecfit run CO/pecfit.json --duo "C:/path/to/duo.exe" --output CO/fit_trial_02
```

## Run HCl with ExoMol states energies

The supplied `HCl/pecfit.json` reads `HCl/pec.txt` (Hartree) and
`HCl/1H-35Cl__HITRAN-HCl.states` directly. All 710 levels are included: J=0..41,
v=0..17, with the available combinations in the states file. The states file
is an energy target, without an assumption that every listed energy was
independently measured. No MARVEL conversion is needed.

```powershell
python -m pecfit plan HCl/pecfit.json
python -m pecfit run HCl/pecfit.json --duo "C:/path/to/duo.exe"
```

For this five-column file, the mapping is ID, energy/cm-1, statistical weight,
J, v. The fitting records use state X, e parity, Lambda=Sigma=Omega=0,
`N=v+1`, and weight `100/sqrt(v+1)/sqrt(J+1)`. The states ID and statistical
weight are **not** the Duo level number or fitting weight. Energies retain their
original zero and precision. `HCl/pecfit.energies.txt` is a readable conversion;
the pipeline rereads the original states file and hashes it for resume checks.

Both `experiment.robust` and `experiment.polish_robust` are 0.0001 in the HCl
configuration. Validation uses no robust reweighting and includes every level.
The isotope masses are 1.00782503223 and 34.968852682 u, matching the existing
HCl Duo output. The initial radial grid extends 0.7..5.0 Angstrom, with 501
points and 80 vibrational functions. Convergence checks also extend the radial
boundaries by 0.05 Angstrom inward and 0.5 Angstrom outward. De remains fitted;
no experimental dissociation energy was supplied.

Create the same input setup in a **new** configuration with:

```powershell
python -m pecfit init --config HCl/new_fit.json --atoms H Cl --masses 1.00782503223 34.968852682 --pec HCl/pec.txt --states HCl/1H-35Cl__HITRAN-HCl.states --grid 0.7 5.0 501
```

Then set `experiment.polish_robust` to 0.0001 and the boundary extensions as in
the supplied HCl configuration. The general default polishing robust remains 0
for the CO procedure. To export only the Duo energy table:

```powershell
python -m pecfit convert-states HCl/1H-35Cl__HITRAN-HCl.states HCl/new_energies.txt
```

The states reader supports plain text and `.bz2`. Its default J/v columns are
4/5; override with `--j-column`/`--v-column` (one based), or `states_j_column` /
`states_v_column` in JSON, for another explicitly identified single-state file.
`--weight-scale` / `states_weight_scale` replaces the numerator 100. Column
positions are not guessed, and extra states, parity branches or non-singlet
quantum numbers must not be passed through this single-state mapping. Duplicate
IDs, duplicate (J,v), invalid numbers and missing columns stop conversion.

## What happens automatically

1. Validate units, EMO structure, a single electronic state, the energy origin,
   quantum numbers, duplicate observations, and basis size. For this state,
   `N=v+1`, Lambda=Sigma=Omega=0, and e parity are required. Total +/- parity is
   also accepted if it agrees with J. Invalid assignments stop the run with an
   explicit diagnostic instead of being guessed from nearby energies.
2. Fit the ab initio grid with `vmax 1`, J=0, ab initio `fit_factor 1e7`, energy
   `fit_factor 1e-10`, and `fit_scale 0.1`. Only the computable v=0 reference
   is included in this inexpensive stage. Full-precision parameters are read
   from Duo's `Parameters:` sections, never from the rounded uncertainty table.
3. Require stability <= 1e-7 and configurable ab initio residual limits. Defaults
   are RMS <= 100 cm-1 within 0.15 Angstrom of Re, outer RMS <= 2000 cm-1, and
   inner maximum residual <= 100000 cm-1. These are permissive starting-model
   checks; they are not the final spectroscopic accuracy target. The supplied
   CO starting fit has about 55 cm-1 RMS in that relatively broad Re window.
4. Refine at Jmax = 10, 40, then the full observed Jmax (123 for CO). Starting
   with some rotational data helps constrain Re; a J=0-only stage can be added
   to the configuration if desired. The
   final stage always includes all supplied levels, even if `all` is omitted
   from the configured schedule. CO uses 80 retained vibrational functions.
   **Duo's `vmax` is a count:** 80 retains v=0..79; it must exceed the largest
   observed v (41 in the CO data).
5. Use the README's experimental weights: ab initio factor 1e-7, energy factor
   1e6, robust 0.0001, initially scale 0.01, increasing toward 0.1 after accepted
   steps. Runs are bounded chunks with automatic continuation. Divergent,
   non-finite, non-monotonic or worsening candidates are rejected; smaller-step
   retries start from the last accepted parameters.
   After at least 20 iterations, a monitor stops a fit once three complete iteration summaries meet the
   configured stability tolerance, instead of waiting for Duo's tighter
   internal tolerance. It retains the corresponding full-precision parameter
   checkpoint and validates it with a new normal Duo calculation. Raw diagnostic
   files from a monitored fit can end mid-iteration; use the following `check`
   directory's diagnostics or the final residual table. `completed.json` records
   the controlled stop and the selected iteration explicitly.
   The minimum avoids a false early plateau while newly added zero-valued EMO
   coefficients begin moving (visible in the supplied `CO_duo_fit_05.out`).
   If Duo converges on its weighted objective but worsens the selected acceptance
   metric, the previous model is retained; smaller steps are reserved for
   nonconverged overshoots and failed candidates.
6. Enable `THRESH_OBS-CALC 10` only when **all positive-weight residuals are
   already below 10 cm-1**. This avoids excluding most data while the starting
   model is still inaccurate. Every candidate is checked using a separate
   `itmax 0`, `robust 0` calculation with no cutoff. Missing or mismatched
   assignments fail that check. Both unweighted and original-weight RMS are
   recomputed on the same observations, regardless of Duo's changing robust weights.
7. If the requested target is still missed, try increasing NR and adding
   zero-initialized coefficients, up to `max_right_order`. For CO the default
   maximum is 10, corresponding to the higher-order examples. Higher orders
   must improve the selected metric and pass the same shape test.
   Before considering higher orders, a full-range polishing stage switches
   robust to `experiment.polish_robust` (default 0) as suggested in the README.
   Basis refits use this polishing setting as well.
8. Recalculate all observations with larger vibrational and radial bases.
   Defaults add 20 vibrational functions and 100 radial points; require a
   maximum energy change <= 0.001 cm-1. If needed, promote the larger basis,
   refit, and repeat up to the configured limit. Radial boundary extensions
   can also be configured.

Fit chunks, retries, restarts, polynomial order, and basis refinements all have
finite limits. A numerical error produces exit code 1 and `failure.json`.
Ctrl-C kills the active Duo process and returns 130. A completed calculation
that misses the requested residual or basis target returns **2**, with status
`needs_review` and the best validated result. Exit code **0** means the selected
accuracy target and the basis/shape/coverage checks passed. Automation does not
guarantee that a chosen model can reproduce every molecule's data.

## Metrics and physical conventions

Both RMS measures are reported in cm-1, over the original positive-weight data:

```text
rms_cm          = sqrt(sum((Eobs-Ecalc)^2) / N)
weighted_rms_cm = sqrt(sum(w*(Eobs-Ecalc)^2) / sum(w))
```

The default target is **unweighted `rms_cm <= 0.1`**, following the energy RMS
criterion in the CO instructions. Duo still optimises the supplied weighted
least-squares objective. Set `experiment.target_metric` to `weighted_rms_cm`
if the acceptance target should use those weights instead. Weighted RMS passing
does not imply unweighted RMS passing. The target metric also determines whether
a new candidate is accepted. A fixed metric is
essential when comparing fits with changing robust weights. Duo's own normalized
fit statistics and stability are retained in `progress.json`.

Weights in the input table are consumed exactly as supplied; no assumption that
they are inverse variances is imposed. Convert uncertainties to your intended
weights before creating the table. Zero-weight observations remain in the
residual and coverage reports but are excluded from RMS statistics.

For Hartree data, the converter applies:

```text
(E - min(E)) * 219474.6313708
```

All radii must be in Angstrom. The pipeline fixes V0=0. Duo's EMO convention is
`V(r) = V0 + (DE-V0) * (1-exp(-beta(r)*(r-Re)))^2`: DE is the asymptote relative
to the chosen origin. The Surkus power and polynomial order use the left or
right branch at RREF; nonpositive RREF means Re. These conventions match the
local Duo `functions.f90` implementation.

If a trusted experimental **De** is available, set `experimental_de_cm`; it is
applied after the ab initio stage and held fixed throughout refinement. Supply
the well depth De, not the dissociation energy D0 from v=0. The CO default leaves
DE fitted because the supplied instructions do not prescribe a numerical
experimental value. The pipeline does not infer one from the sample outputs.

## Reuse for another molecule

If you already have a matching single-state Duo fitting input:

```powershell
python -m pecfit init --template HCl/my_fit.inp --config HCl/pecfit.json
python -m pecfit plan HCl/pecfit.json
python -m pecfit run HCl/pecfit.json --duo "C:/path/to/duo.exe"
```

Or start from a two-column ab initio PEC and a ten-column experimental energy
table. This generates an initial EMO template and JSON configuration:

```powershell
python -m pecfit init --config HCl/pecfit.json --atoms H Cl --masses 1.00782503223 34.968852682 --pec HCl/pec.txt --units hartree --energies HCl/energies.txt --grid 0.8 4.0 401 --de 37000
```

`--de` here is only an **initial guess**. To fix a measured De during refinement,
use `experimental_de_cm` in JSON. When masses are omitted, Duo uses its atomic
defaults. Supply isotope masses explicitly for a specific isotopologue.

The energy table has no header, allows `#`/`!` comments, and follows the README:

```text
# J parity N energy/cm-1 state v Lambda Sigma Omega weight
0 e 1 0.0      X 0 0 0 0 1.0
0 e 2 2885.9   X 1 0 0 0 1.0
1 e 1 20.8     X 0 0 0 0 1.0
```

These numbers illustrate the format only. Use your actual measured term values.
Include the J=0,v=0 zero reference. The experimental table must already describe
the intended isotopologue and common energy origin. For ExoMol single-state
files use the explicit `--states` mapping above. Raw MARVEL formats and arbitrary
multi-state ExoMol files require separate preparation.

To use a new PEC with an existing fitting template, pass both `--template` and
`--pec` to `init`, or set `pec_file`, `pec_units`, and `pec_shift_minimum` in JSON.
Already shifted cm-1 data can use `pec_shift_minimum: false`. The template's
prepared ab initio cm-1 grid is preserved without rescaling by default.

## Settings and outputs

`init` writes all available settings; `CO/pecfit.json` specifies only overrides.
Defaults are documented in `pecfit/pipeline.py::DEFAULTS`. Unknown keys fail
validation to catch spelling mistakes. Useful settings include:

| Setting | Meaning |
|---|---|
| `abinitio.iterations`, `experiment.iterations` | Iterations per Duo chunk |
| `max_restarts` in each stage | Maximum continuation chunks |
| `experiment.max_retries` | Smaller-step retries after a rejected chunk |
| `experiment.robust` | Robust parameter; 0 disables robust fitting |
| `experiment.polish_robust` | Robust parameter for final full-range polishing and basis refits; default 0 |
| `experiment.outlier_threshold_cm` | Outlier cutoff; 0 disables it |
| `experiment.minimum_iterations` | Minimum iteration before monitored stability stopping; default 20 |
| `experiment.max_right_order` | Maximum NR; set to initial NR to disable growth |
| `refine_only` | Start from a prepared analytic potential instead of repeating the ab initio stage |
| `experiment.bob_prefit` | Fit BOB coefficients separately before the joint refinement |
| `experiment.line_search_steps` | Backtrack unsafe or worsening parameter updates; 0 disables backtracking |
| `validation.shape_range_angstrom` | Explicit radial interval for PEC/BOB checks, independent of the calculation box |
| `validation.check_centrifugal_wells` | Reject extra wells/barriers in the effective potential for each observed J; requires explicit masses |
| `validation.bob_max_abs`, `bob_max_slope`, `bob_max_turns` | Bounds on the dimensionless BOB curve and its variation |
| `validation.bob_turning_tolerance` | Ignore reversals smaller than this dimensionless excursion when counting BOB lobes; default 0 |
| `validation.range_extension_angstrom` | Amount to extend inner and outer grid boundaries during convergence checks |
| `threads` | OMP, MKL and OpenBLAS thread count; default 1 |
| `states_file` | Optional ExoMol states source, relative to config; replaces the template energy table |
| `states_j_column`, `states_v_column` | One-based J/v columns; defaults 4/5 |
| `states_weight_scale` | Numerator of scale/sqrt(v+1)/sqrt(J+1); default 100 |
| `timeout_seconds` | Maximum time per Duo process |

`final.inp` is ready for a Duo straight-through calculation. The result folder
also contains `final_parameters.json`, `potential.csv`, `final_residuals.csv`,
`basis_convergence.csv`, a machine-readable `result.json`, and `report.md`.
All attempts retain their input, output, stderr, `.en`, `.pot`, and
`Potential_functions.dat`. Separate directories avoid filename collisions.
`progress.json` records accepted and rejected attempts; `manifest.json` stores
configuration and SHA-256 hashes of sources, code, and the executable.

The `.en` format used by the supplied Duo builds prints energies to 0.0001 cm-1;
the validation metrics share that resolution. Shape checks sample 2001 radial
points plus Re on the declared interval. They detect additional turning points
there, but do not establish physical extrapolation beyond the tested interval,
parameter uncertainty, or predictive accuracy on independent data. The radial
boundary extension defaults to zero: enlarge it explicitly to test boundaries.
Changing basis size alone does not test the suitability of the radial domain.

## Tests

```powershell
python -m unittest discover -s tests -v
$env:PEC_FIT_DUO = "C:/path/to/duo.exe"
python -m unittest discover -s tests -v
```

The optional real-Duo test verifies straight-through energy parsing, assignment
coverage, cache reuse, and rerunning a corrupt diagnostic. Unit tests cover
Hartree conversion, EMO branches/asymptote, model shape, full-precision parameter
extraction, observation coverage, fixed-weight statistics, and safe resume.
The full CO run is the scientific integration test; its output is deliberately
not used as a silently updated numerical reference.

## HCl model comparisons and localized levels

`HCl/models/` contains the requested NR=10 EMO trial (missing B10 initialized
to zero), MLJ trials, and models with rotational BOB. MLJ parameters and BOB
parameters have separate namespaces internally; repeated names such as RE and
B0 cannot overwrite each other. Higher EMO order can be promoted automatically;
MLJ degree and structural parameters are chosen explicitly in each template.

The HCl test revealed that `assign_v_by_count` can assign a continuum box state
to a supplied high-J vibrational label. Increasing the box can insert additional
continuum states below a localized quasibound level. A low fit RMS alone therefore
does not validate those assignments. The original rank-based workflow continues
to reject failed basis checks.

The optional localized refiner uses the same single-state sinc-DVR Hamiltonian,
selects levels by inner-well node count and probability, and verifies all final
energies against native Duo's 12-decimal rovibronic output. It automatically
enlarges Duo's contracted vibrational basis when needed to match the primitive
calculation, then checks two larger radial boxes. Its finite-box quasibound
energies are stabilization estimates, not resonance poles. A failed full-data
box check remains `needs_review`; stable bound-level results are reported
separately.

Install the optional numerical and plotting dependencies:

```sh
python -m pip install ".[localized,plots]"
```

Refine all source levels, retaining the configured robust parameter:

```sh
python -m pecfit refine-localized HCl/models/physical_emo14_polish.json --duo /path/to/duo --optimizer penalized --robust 0.0001 --output HCl/pecfit_runs/new_all_levels
```

For the explicitly labelled conservative HCl calibration used in this study:

```sh
python -m pecfit refine-localized HCl/models/physical_emo14_polish.json --duo /path/to/duo --optimizer penalized --robust 0 --fit-below 35000 --output HCl/pecfit_runs/new_bound_comparison
```

This cutoff selects 649 calibration levels. All 710 original levels remain in
the final diagnostic tables; calibration, bound-state, and full-data errors have
different denominators and are labelled separately. The number classified as
bound is model dependent. `--robust 0` runs a separate fixed-weight comparison;
the default uses the configuration's `Robust=0.0001`, with the Watson alpha
0.001 and normalized input weights used by Duo. No level is silently discarded.

`--evaluate-only` validates a supplied template without fitting. Each localized
run requires a new output directory. Unlike the core Duo runner it does not
support cache replay/resume; `progress.json` retains every accepted checkpoint.
`final.inp` uses the audited energy ranks, and `final_residuals.csv` preserves
their mapping to the original physical v labels. Regenerate this input after
changing the potential or radial grid. The source states file is never relabelled.

For HCl, the localized trials fix Re to the published isotopologue value
1.27454677 Angstrom to reduce its correlation with the BOB constant. The BOB
function is an effective correction for this isotopologue; this does not infer
a mass-dependent correction for other isotopologues. `--preserve-c6` adjusts
MLJ PHIINF as De changes to retain the seed's C6. It is not applicable to EMO.

Use `compare_pec_models.py` to collect completed core or localized runs on the
same source energies and weights. Its automatic full-data selection excludes
models with failed shape or basis checks. `--run` invokes the core Duo workflow;
run localized refinements with the separate command shown above.

See `HCl/MODEL_REFINEMENT.md` for the actual results, literature provenance,
shape plots, and numerical limitations. The published HCl MLR3 model is not the
same function as Duo MLJ, and the KH coefficients supplied as an example were
not used as HCl spectroscopic parameters.

The final physically screened configurations are `physical_mlj13_bob6.json`,
`physical_emo14_bob6.json`, and their separately labelled `_polish` variants.
Use `--optimizer penalized` to guide trial steps with soft shape constraints;
every returned candidate must still pass the strict shape checks. The guarded
optimizer instead backtracks every parameter step. Both record the best
physical checkpoint and retain the original observation weights for reporting.
The penalized optimizer includes a weak prior on changes to the seed tail;
its sampled shape penalties guide optimization but never replace validation.

`validation.potential_ceiling_cm` optionally rejects a monotone but excessive
repulsive wall. The HCl physical trials use 1e8 cm-1 on 0.5--12 Angstrom. This is
a broad guard against pathological extrapolation and numerical ill-conditioning,
not an experimental constraint on the wall. It must be considered again when
reusing the workflow for another molecule. The high-order HCl trials which
violate this limit were not retained as physical results.

Compare completed trials without rerunning them:

```sh
python compare_pec_models.py HCl/models/*.json --output HCl/pecfit_runs/comparison --calibration-max-energy 35000 --potential-ceiling 1e8 --plots
```

The common-subset selector cannot certify the full dataset. The comparison
reports `best` for the full dataset and `best_calibration` for the explicitly
specified common subset, with independent basis/shape eligibility checks.
