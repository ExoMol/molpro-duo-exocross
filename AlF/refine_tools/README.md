# AlF refinement tools

This driver builds a two-state X1Sigma+/A1Pi model. The A potential uses `coupled-pec`, with sub-types `EMO repulsive_exp EMO` and the lower adiabatic component. The AlF-specific long-range estimates, barrier target, and their limitations are in `../LITERATURE.md`.

Install the scientific dependencies in your Python environment:

```powershell
python -m pip install -r AlF/refine_tools/requirements.txt
```

Run the complete workflow from the repository root, using a new output directory:

```powershell
python AlF/refine_tools/pipeline.py AlF/reference_inputs/27AlF_X1E+_A1Pi_EMO_02.inp --duo C:/sergei/programs/Duo/build/ifx/duo.exe --output AlF/pecfit_runs/new_fit
```

The workflow audits the source, constructs the barrier seed, fits to experiment, checks radial-grid/basis/box convergence, evaluates every retained original observation, and writes `final.inp`, plots, and JSON diagnostics. A full run takes several minutes or longer depending on hardware and the number of refinement rounds. `--ordinary-rounds` and `--robust-rounds` control the outer loops. Publication-quality residuals and final parameters must be reviewed together with the flagged experimental levels.

For a fresh evaluation of the delivered model:

```powershell
python AlF/refine_tools/run_native.py AlF/refined_model/27AlF_X_A_coupled_pec.inp AlF/pecfit_runs/check --duo C:/sergei/programs/Duo/build/ifx/duo.exe --jmax 98 --precise
python AlF/refine_tools/coupled.py AlF/pecfit_runs/check --jmax 98
python AlF/refine_tools/check_curves.py AlF/pecfit_runs/check
```

To refine further from an existing model:

```powershell
python AlF/refine_tools/refine.py AlF/refined_model/27AlF_X_A_coupled_pec.inp --duo C:/sergei/programs/Duo/build/ifx/duo.exe --output AlF/pecfit_runs/continued_fit --jmax 98 --rounds 4 --proposal-steps 3 --scale 1 --robust 0.00001 --criterion cauchy
```

The delivered grid is larger than the initial fitting grid, so this command may run more slowly. Every proposed parameter set is checked for sensible PEC/BOB/coupling shapes and then recalculated with native Duo before acceptance. `weighted_rms` is the alternative acceptance criterion for ordinary least squares. Cauchy acceptance may increase ordinary RMS while improving the robust objective; both are retained in the history.

## Data handling

The original file is retained unchanged. `input_audit.json` records character normalization, the missing Lx morphing reference, the singlet-Pi Omega correction, the two source records already disabled by negative J, the extension of JLIST to 98, and six conflicting duplicate assignments. Duplicate levels are represented by their weighted mean and summed weight during fitting. Final residuals expand these means back to all original records. The original twelve zero-weight levels remain zero-weight. No additional observation is rejected because its residual is large.

The duplicates may contain parity transcription errors. The preparation deliberately does not guess different parities. A robust loss on the duplicate means is not mathematically identical to a robust loss on separate duplicate measurements; the ordinary weighted squared loss differs only by a fixed within-group constant.

## Safeguards and reproducibility

- Compound parameter labels can repeat. Parameter transfer uses object, state indices, and position; it never identifies a coefficient by `B6` alone.
- Native parameter tables contain *proposed* values. Their preceding residual table is not evidence for those values. Each accepted proposal receives a fresh zero-iteration evaluation.
- Restarted robust fits carry Watson IRLS weights explicitly. Native Duo updates robust weights after its first step, so restarting one-step fits with only `ROBUST` would otherwise repeat ordinary least squares.
- Every run records input/executable hashes and preserves stdout. The refinement also snapshots its Python source.
- The selected PECs are checked from 0.7 to 1000 Angstrom. X must have one well; A must have one inner well and one barrier, with at most one shallow outer well. Equal PL/PR = 6 or 8 is enforced. Small rotational and electronic angular-momentum corrections are bounded.
- `convergence.py` compares 601/801 radial points, vmax 60/90, and a final 1001-point, 8-Angstrom, vmax-100 calculation on the same levels. Its pass tolerance is 0.0001 cm-1 per positive-weight level.
- Scientific tests: `python -m unittest discover -s AlF/refine_tools -p test_refinement.py -v`.

The prepared inputs use integer J, numeric state labels, and the supplied AlF level conventions. The driver, literature constants, and shape bounds are specific to AlF; revise them explicitly before applying this barrier model to a different molecule. NumPy, SciPy, Matplotlib, and a Duo executable that supports `coupled-pec`/`repulsive_exp` are required. The input has no intensity calculation, so its inherited dipole functions are not validated by this energy fit.
