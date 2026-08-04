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
set grid = $4
set core = $5


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

BASIS=$basis
ANSATZ=3c,fix=1,canonical=1


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


