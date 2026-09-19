# AlF physical assumptions and literature

The model has exactly two electronic states: X1Sigma+ (Lambda = 0, multiplicity 1) and A1Pi (Lambda = 1, multiplicity 1). The CH example supplies the **functional form**, not AlF numerical parameters or symmetry.

## Fixed long-range information

[Qin, Bai & Liu, MNRAS 510, 3011 (2022), Table 1 and equations 7–9](https://doi.org/10.1093/mnras/stab3598) assign X and A to the same ground-atom Al(2P) + F(2P) limit. They estimate C6 = 75.23 with the London formula and fit C5 = 8.07 for A. These are approximate theoretical estimates, not precision spectroscopic determinations. The table does not explicitly label coefficient units: the atomic-unit interpretation of C6 follows their London formula and stated atomic-unit polarizabilities; C5 is interpreted in the consistent atomic-unit convention.

Using Eh = 219474.63136320 cm-1 and a0 = 0.529177210903 Angstrom gives:

| Coefficient in the AlF input | Value | Units |
|---|---:|---|
| B5 = -C5 | -73495.9702828864 | cm-1 Angstrom^5 |
| B6 = -C6 | -362561.9198027829 | cm-1 Angstrom^6 |

The supplied absolute dissociation energy, 55564.02471 cm-1 relative to the X potential minimum, is retained. It agrees with the fixed AlF value in [Yousefi & Bernath, ApJS 237, 8 (2018)](https://doi.org/10.3847/1538-4365/aacc6a). It is not refitted from the low vibrational levels.

## Barrier and coupling

[Andreazza & de Almeida, MNRAS 437, 2932 (2014)](https://academic.oup.com/mnras/article/437/3/2932/1034811) report an A-state barrier of 0.763 eV at 4.8 bohr, using earlier electronic-structure curves. These correspond to approximately 6154.01 cm-1 above dissociation and 2.54005 Angstrom. They are initialization targets, not additional measured levels.

The native Duo implementation used here is:

```
Vrep(r) = D + A exp(-delta/r)/r^gamma + sum_k Bk/r^k
VA(r) = (VEMO(r) + Vrep(r) - sqrt((VEMO(r)-Vrep(r))^2 + 4 W(r)^2))/2
```

W is a second EMO with zero asymptote. `COMPON 1` selects the lower branch. `Nparameters` describes sub-functions within **one** electronic PEC; it does not add electronic states to the Hamiltonian.

Gamma is fixed at 8, so the positive repulsive term decays faster than the retained r^-5 and r^-6 attractions. Copying CH's gamma = 5.5 would make that term dominate C6 at sufficiently large distance. The diabatic EMO asymptote (80000 cm-1), repulsive amplitude, delta, and W parameters are construction choices constrained by the barrier target and smoothness. They are not independently measured AlF constants. The shallow outer van der Waals well follows from a barrier followed by an attractive tail.

## Scope of the empirical fit

The observed levels sample X v = 0–9 and A v = 0–6; they cannot uniquely determine the barrier or long-range tail. PL and PR are equal to 6 in all EMO sub-functions. X remains an EMO: its limiting energy is fixed, but its exponential tail is not a dispersion expansion.

The interaction of A1Pi with b3Sigma+ is documented by [Doppelbauer et al., Molecular Physics (2020)](https://doi.org/10.1080/00268976.2020.1810351). A two-state X/A model cannot explicitly represent those triplet perturbations. Large or localized residuals are flagged for data/perturbation review; their cause is not assigned automatically, and they remain in reported all-data statistics.
