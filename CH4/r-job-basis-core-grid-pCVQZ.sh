#!/bin/csh
#
# Generation of the input file for MOLPRO
#

#module load mkl/10.2.5/035
#module load compilers/intel/11.1/072
#module load mpi/openmpi/1.4.1/intel
#module load molpro/2009.1/openmpi/intel

cd $2
set pwd = `pwd`
echo $pwd

set basis = $3
set core = $4
set grid = $5


set irun = $1
set ffile = "./inp/ch4_grid{$grid}.dat"


#./read_geom.sh 1

set fname = ch4.${irun}.$basis.$core.$grid.kroll
#set xyzfile = ch4.${irun}.xyz3.dat
#set datfile = ch4.${irun}.d3.dat

echo $fname
#
cat<<endb> $fname.com

gthresh,twoint=1.d-10,prefac=1.d-18,zero=1.d-10,energy=1.d-10;
memory,512,m

PROC CCT-f12
 {hf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,18;wf,orbital}
 {ccsd(t)-f12b;$core; thresh,energy=1.d-10,coeff=1.d-9,thrint=1.d-10,zero=1.d-10;maxit,100}
 !{OPTG,grad=5,energy=8,saveact='xxx.opt';print,history;inactive,alpha}
 !{readvar,'xxx.opt'}
ENDPROC

PROC calc-mp2
 {hf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,18;wf,orbital}
 {mp2}
ENDPROC

!BASIS=$basis


basis={
!
! CARBON       (16s,10p,5d,3f,2g) -> [8s,8p,5d,3f,2g]
! CARBON       (16s,10p,5d,3f,2g) -> [8s,8p,5d,3f,2g]
s, C , 96770.0000000, 14500.0000000, 3300.0000000, 935.8000000, 306.2000000, 111.3000000, 43.9000000, 18.4000000, 8.0540000, 3.6370000, 1.6560000, 0.6333000, 0.2545000, 0.1019000, 0.0394000, 4.2560220
c, 1.10, 0.0000250, 0.0001900, 0.0010000, 0.0041830, 0.0148590, 0.0453010, 0.1165040, 0.2402490, 0.3587990, 0.2939410
c, 1.10, -0.0000050, -0.0000410, -0.0002130, -0.0008970, -0.0031870, -0.0099610, -0.0263750, -0.0600010, -0.1068250, -0.1441660
c, 11.11, 1
c, 12.12, 1
c, 13.13, 1
c, 14.14, 1
c, 15.15, 1
c, 16.16, 1
p, C , 101.8000000, 24.0400000, 7.5710000, 2.7320000, 1.0850000, 0.4496000, 0.1876000, 0.0760600, 0.0272000, 7.9358560
c, 1.3, 0.0008910, 0.0069760, 0.0316690
c, 4.4, 1
c, 5.5, 1
c, 6.6, 1
c, 7.7, 1
c, 8.8, 1
c, 9.9, 1
c, 10.10, 1
d, C , 3.7056000, 1.4212000, 0.5451000, 0.2091000, 15.8727840
c, 1.1, 1
c, 2.2, 1
c, 3.3, 1
c, 4.4, 1
c, 5.5, 1
f, C , 1.4438000, 0.5931000, 0.2436000
c, 1.1, 1
c, 2.2, 1
c, 3.3, 1
g, C , 1.1825000, 0.4685000
c, 1.1, 1
c, 2.2, 1
h=$basis
ANSATZ=3c,fix=1,canonical=1
}





dkroll=1

symmetry,nosym;
orient,noorient;

geometry={angstrom;
 C;
 H, 1,xR1;
 H, 1,xR2, 2,xA1;
 H, 1,xR3, 2,xA2, 3,xA4, 1;
 H, 1,xR4, 3,xA3, 2,xA5, 1;
}

! C;
! H, 1,xR1;
! H, 1,xR2, 2,xA12;
! H, 1,xR3, 2,xA13, 3,xA23, 1;
! H, 1,xR4, 3,xA24, 2,xA14, 1;

i1 = 1
field=[-0.005, 0.005]

`awk -f print_geom.awk irow=$irun $ffile`

 R1CH(i1) = xR1
 R2CH(i1) = xR2
 R3CH(i1) = xR3
 R4CH(i1) = xR4
 A12(i1)  = xA1
 A13(i1)  = xA2
 A24(i1)  = xA3
 A23(i1)  = xA4
 A14(i1)  = xA5
 !
 CCT-f12
 !
 eeee(i1)=energy
 !
 NNZ(i1) = 'XX'
 IZ(i1) = ${irun}
 !
 table,NNZ,IZ,R1CH,R2CH,R3CH,R4CH,A12,A13,A24,A23,A14,eeee
 DIGITS,0,  0,   8,   8,   8,   8,  8,  8,  8,  8,  8,  12
 !
 !title PES of CH4
 !
---

endb

setenv OMP_NUM_THREADS 1

limit
#echo $TMP

setenv wdir $TMPDIR


# Run Molpro
echo "System TMPDIR = $TMPDIR"
echo "wdir = $wdir"

cd $wdir

setenv molproexe /shared/ucl/apps/Molpro/2009.1/bin/
#setenv PATH /shared/ucl/apps/Molpro/2009.1/bin/:"${PATH}"

echo "Running molpro -n 1 -d $TMPDIR -W $wdir < $pwd/$fname.com > $pwd/$fname.out"

$molproexe/molpro -n 1 -I $TMPDIR -d $TMPDIR -W $wdir < $pwd/$fname.com > $pwd/$fname.out


#if (-e $xyzfile) then
#   /bin/cp  *.dat  $pwd
#endif


#/bin/rm -rf $HOME/${fname}.codine
#/bin/rm -rf /usr/scratch/$USER/${fname}


