# molpro-duo-exocross

A pipeline for compouting molecular spectra from first principles

## Automatic potential fitting

The Python pipeline fits an EMO or MLJ potential, with optional rotational BOB, to ab initio data and then refines it
against experimental energies with Duo. See [PEC_FITTING.md](PEC_FITTING.md)
for CO commands, validation criteria, and reuse with another X1Sigma+ molecule.

Quick start: `python -m pecfit run CO/pecfit.json --duo /path/to/duo`.
`HCl/pecfit.json` uses ExoMol states energies directly, with the supplied J/v weights and robust 0.0001. Run `python -m pecfit run HCl/pecfit.json --duo /path/to/duo`.

For the HCl model comparison, literature references, shape checks and quasibound-state limitations, see [HCl/MODEL_REFINEMENT.md](HCl/MODEL_REFINEMENT.md). The optional `refine-localized` command tracks localized levels and verifies them with native Duo.

The [C2 eight-state refinement](C2/REFINEMENT_REPORT.md) includes coupled-state fitting tools, audited assignment handling, and numerical validation.

The [AlF two-state refinement](AlF/REFINEMENT_REPORT.md) uses a coupled EMO/repulsive barrier potential, literature-based long-range coefficients, and an automatic fitting and validation pipeline.
