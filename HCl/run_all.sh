#!/bin/bash

# 1. MOLPRO calculations to produce a CO PEC and DMC.
# 1.1 run-molpro.sh contains linux instruction for MOLPRO
# 1.2 CO_loop6.com is a MOLPRO input file for CCSD(T) calculations
# 1.3 The output CO_loop6.out contains the PEC and DMC data which are extracted (grepped) and saved to 
#     pec.txt and dmc.txt
# 2. Duo calculations to produce .trans and .states files 
# 2.1 Script run_duo_template_intensity.sh merge CO_duo_01.template CO_dmc_01.template (basic Duo templates) with  
#     pec.txt and dmc.txt to create a Duo input file CO_duo_01.inp
# 2.2 The Duo calculations are run as 
#     ./j-duo.x < CO_duo_01.inp > CO_duo_01.out
# 2.3 Finally, the line list consists of two files: CO_duo_01.trans and CO_duo_01.states
# 3.  ExoCross calculations 
# 3.1 h_g1_T.sh is used as an example of testing the line list computed by calculating CO cross sections at T= 1000 K
#     using the Gaussian line profile of HWHM of 1 cm-1 on a grid of 1 cm-1 covering the range 0-20000 cm-1

# 1.MOLPRO
#./run-molpro.sh molpro_loop6.com
# 2. Duo
./run_duo_template_intensity.sh duo_01.template duo_dmc_01.template
# 3. ExoCross
./h_g1_T.sh  duo_01 1000
