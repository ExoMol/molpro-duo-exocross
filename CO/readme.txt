1. Ab initio line list
The simple master pipeline file is run_all.sh. This is to generate an ab initio line list without any refinement to experiment.


2. Refining the ab initio model by fitting to experimental/MARVEL energies.
See CO_duo_fit_01.inp

2.1 Convert PEC from Hartree to cm-1 (x219474.6313708) and shift it relative to the minimum: (E_i-min(E_i,i=1..Nmax))*219474.6313708 

2.2 Form two "poten" sections, one with analytic representation using some initial parameters with the keyword "fit" where these parameters will need to be optimased, e.g.:

poten X
name "X1Sigma+"
lambda 0
symmetry +
mult   1
type   EMO
values
V0           0.00000000000000E+00
RE           1.5                   fit
DE           50000                 fit
RREF        -1.00000000000000E+00
PL           6.00000000000000E+00
PR           6.00000000000000E+00
NL           4.00000000000000E+00
NR           8.00000000000000E+00
B0           2.000000000000000000  fit
B1           0.000000000000000000  fit
B2           0.000000000000000000  fit
B3           0.000000000000000000  fit
B4           0.0  fit
B5           0.0  fit
B6           0.0  fit
B7           0.0  fit
B8           0.0  fit
end
The second "poten" should contain the ab initio energies (cm-1) on a grid as follows:

abinitio poten X
name "X1Sigma+"
symmetry +
lambda 0
mult 1
type grid
units cm-1
fit_factor 1e7
Weighting PS1997 1e-5 80000.0
values
        0.80000000        1.25483301E+05
        0.80733333        1.16718737E+05
        0.81466667        1.08550612E+05
        0.82200000        1.00945150E+05
        0.82933333        9.38685749E+04
        0.83666667        8.72871099E+04
....
end

Here 
Weighting PS1997 1e-5 80000.0
is the weighting scheme (see manal) 
and 
fit_factor 1e7  
is to tell Duo to increase the fitting weights of these ab initio data by a huge factor. 

2.3 Add the fitting section containing experimental energies and fitting instructions:
FITTING
JLIST    0 - 0
itmax 5000
fit_factor  1e-10
fit_type  DGELSS
fit_scale  0.1
output   f01
lock     -8000
robust  0
energies  (J Parity N Energy State v Lambda Sigma Omega Weight)
       0  e        1              0   X        0     0       0       0    320000.00000000
       0  e        2    2143.271107   X        1     0       0       0       442.22264949
       0  e        3     4260.06217   X        2     0       0       0       861.74597227
       0  e        4    6350.439055   X        3     0       0       0     10679.19692439
       0  e        5    8414.469295   X        4     0       0       0       285.71428571
       0  e        6    10452.22217   X        5     0       0       0         6.48148148
       0  e        7    12463.76852   X        6     0       0       0         2.61894369
....
end

itmax 5000 
is the maximal number of iterations
fit_factor  1e-10
is to set the fitting weights of the energies to a tiny number (~0) and to effectively exclude them from the fit. 

output   f02 
is the name of the auxiliary files .pot and .en containing residuals of the PECs and energies, respectively. 

2.4 In CONTRACTION set vmax to 1 to reduce the computational costs of energies:
CONTRACTION
  vib
  vmax  1
END

2.5 Run Duo with a larger number of iterations (here 5000) until the convergence ("stability") conv = [rms(i)-rms(i-1)/rms(i) ~ 1e-7-1e-8 

2.6 Check the .pot file (f01.pot) if the residuals are acceptable. Something of the order of 10 cm-1 (around the eqauilibrium) to 1000 cm-1 at the dissociation and up to 10^5 at the very small geometries. These parameters will be used as initial, so it is not too critical to have them very accurately. 

2.6 If the fit looks good, go to the last iteration in the Duo .out file and copy the final values of the parameters 

V0           0.00000000000000E+00
RE           1.12944523974061E+00  fit
DE           7.07666762772840E+04  fit
RREF        -1.00000000000000E+00
PL           6.00000000000000E+00
PR           6.00000000000000E+00
NL           4.00000000000000E+00
NR           8.00000000000000E+00
B0           2.58948420238538E+00  fit
B1          -1.79041687477648E-01  fit
B2          -3.28179835019414E-01  fit
B3           4.10056644064876E-02  fit
B4           7.27135282308566E-02  fit
B5           7.04536428226530E+01  fit
B6          -2.38633968016442E+02  fit
B7           2.80353764398463E+02  fit
B8          -1.11386729145406E+02  fit


and paste into the Duo input file (give it a new name). 

2.7 Fitting to the experimental energies: CO_duo_fit_02.inp

2.7.1 In CONTRACTION set vmax to a large number, e.g. 8 to increase the accuracy of Duo's energies:
CONTRACTION
  vib
  vmax  80
END

2.7.2 Set 
fit_factor 1e-7 
in  abinitio poten X
abinitio poten X
name "X1Sigma+"
symmetry +
lambda 0
mult 1
type grid
units cm-1
fit_factor 1e-7
Weighting PS1997 1e-5 80000.0
values
        0.80000000        1.25483301E+05
        0.80733333        1.16718737E+05
....

2.7.3 
Set 
fit_factor  1e6
in the Fitting section:
FITTING
JLIST    0 - 10
itmax 1000
fit_factor  1e6
fit_type  DGELSS
fit_scale  0.01
output   f04
lock     -8000
robust  0.0001
energies  (J Parity N Energy State v Lambda Sigma Omega Weight)
       0  e        1              0   X        0     0       0       0    320000.00000000
       0  e        2    2143.271107   X        1     0       0       0       442.22264949
....

2.7.4 Set 
fit_scale  0.01
to a value 0.01 - 0.1 
to make the Newton-fit less aggressive. The change of the parameters estimated using the Newton optimisation will be scaled by this value (0.01). The default is 1. 

2.7.5 Ideally, change the parameter De to an experimental value and exclude it from the fit by deleting "FIT"

2.7.6 Decide on the  range of the J values in JLIST. Usually we start with a small Jmax, 0,1-10 and then extend to the full experimental range. 

2.7.7. Make sure the experimental energies are in the right form and are matched to the theoretical states correctly before doing the refinement. To this end it is useful to inspect the .en file in order to check for any wrongly assigned and therefore wrongly placed states. 

2.8 Run the fit and hope it does not blow because of some parameter values getting crazy. The target accuracy is an rms = 0.01 - 0.1 cm-1 or better. 

2.8.1 If the run is successful (the residuals are significantly smaller), copy paste the fitted values from the latest iteration back into the Duo input file (change the name of the file as well as the name of "output" in "FITTING"), increase Jmax or run it again if the convergence has not been reached within the itmax limit. 
See CO_duo_fit_03.inp

2.8.2 Experiment with the basis set sizes, e.g. increase vmax to see if it improves the rms or just produces better convergence of the energies. Kill the job if there no improvement or small, copy/pase the finale parameter if they are sensible, change other inputs and restart. 

2.8.3 Keep an eye on accidental mismatches of the Duo energies and experimental values caused  by accidental misassignments of the Duo states - correct them if required after each run. In order to prevent any harms during the calculations from outliers caused by such mismatches, use the threshold keyword 
THRESH_OBS-CALC 10
to exclude outliers (>10 cm-1 here) from the fit. 

2.8.4 Try increasing the number of the fitting parameters 
CO_duo_fit_05.inp
but be careful not to overfit. An overfit could be recognised by an unphysical (not Morse-like)shape if the potential.
To check, run Duo-fit with itmax set to zero (straight-through) with the keyword PRINT_PECS_AND_COUPLINGS_TO_FILE included anywhere in the input file 
(see CO_duo_fit_06.inp) 
and plot Potential_functions.dat produced. 



2.8.5 Try setting the  keyword robust to zero (to switch off the dynamic re-weighting procedure Robust) or changing to a smaller/large value
robust 0 

















