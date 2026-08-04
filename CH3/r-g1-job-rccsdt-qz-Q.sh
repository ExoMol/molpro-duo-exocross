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

set grid = 1


set basis = aug-cc-pVQZ
set core = ""

set outdir = qz-Q-g$grid


set irun = $1
set ffile = "./grid{$grid}.dat"


if ( ! -e $pwd/${outdir}-dat) then
   mkdir $pwd/${outdir}-dat
endif

if ( ! -e $pwd/${outdir}-out) then
   mkdir $pwd/${outdir}-out
endif


#./read_geom.sh 1

set fname = CH3-rccsdt-qz-Q.${irun}.$outdir
#set xyzfile = ch4.${irun}.xyz3.dat
#set datfile = ch4.${irun}.d3.dat

echo $fname
#
cat<<endb> $fname.com

memory,1000,m
gthresh,energy=1.d-10,zero=1.d-14,thrint=1.d-14,oneint=1.d-14,twoint=1.d-14,prefac=1.d-20


PROC CCT
{hf,ipnit=1,maxdis=30;wf,9,1,1;maxit,100;accu,16}
{rccsd(t),thrden=1.0d-08}
ENDPROC

PROC calc-mp2
 {hf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,18;wf,orbital}
 {mp2}
ENDPROC



noorient;
symmetry,nosym
geometry={
  C;
  H, 1, a1;
  H, 1, a2, 2, b3;
  H, 1, a3, 2, b2, 3, b1, 1
}

basis=aug-cc-pVQZ

i1 = 1
field=[-0.005, 0.005]

`awk -f print_geom.awk irow=$irun $ffile`

a1 =xR1
a2 =xR2
a3 =xR3
b1 =xA1
b2 =xA2
b3 =xA3

NNZ(i1) = 'ppp'
IZ(i1) = ${irun}
!
r1(i1) = a1
r2(i1) = a2
r3(i1) = a3
r4(i1) = b1
r5(i1) = b2
r6(i1) = b3
!


  gexpec, dmx,dmy,dmz,QM


  
  !
  CCT
  !
  !calc-mp2
  !
  eeee(i1)=energy
  !
  quadxx(ii)=qmxx !save quadrupole momemts
  quadyy(ii)=qmyy
  quadzz(ii)=qmzz
  !
  NNZ(i1) = 'ppp'
  IZ(i1) = ${irun}
  !
  table,NNZ,IZ,r1,r2,r3,r4,r5,r6,eeee,quadxx,quadyy,quadzz
  DIGITS,0,  0, 8, 8, 8, 8, 8, 8,  12,    12,    12,    12
  !
  save,$fname.dat,new
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



#module load mkl/10.2.5/035
#module load compilers/intel/11.1/072
#module load mpi/openmpi/1.4.1/intel

module load mkl/10.2.5/035
module load compilers/intel/11.1/072
module load mpi/openmpi/1.4.1/intel
setenv molproexe /shared/ucl/apps/Molpro/2012.1-bindist/bin/

#setenv molproexe /shared/ucl/apps/Molpro/2009.1/bin/
#setenv PATH /shared/ucl/apps/Molpro/2009.1/bin/:"${PATH}"


echo "Running molpro -n 1 -d $TMPDIR -W $wdir < $pwd/$fname.com > $pwd/$fname.out"

time $molproexe/molpro -n 1 -I $TMPDIR -d $TMPDIR -W $wdir < $pwd/${fname}.com > & $pwd/${fname}.out


#$molproexe/molpro -n 1 -I $TMPDIR -d $TMPDIR -W $wdir < $pwd/${fname}.com > $wdir/${fname}.out

pwd
ls

if (-e $fname.dat) then
   /bin/cp  $fname.dat  $pwd/${outdir}-dat
endif

if (-e $fname.xyz.dat) then
   /bin/cp  $fname.xyz.dat  $pwd/${outdir}-dat
endif

if (-e $pwd/$fname.out) then
    gzip $pwd/$fname.out
   /bin/mv  $pwd/$fname.out.gz  $pwd/${outdir}-out
endif

if (-e $pwd/$fname.com) then
   /bin/rm  $pwd/$fname.com
endif

