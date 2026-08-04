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


set basis = TZ
set core = ""

set outdir = pol-tz-g$grid


set irun = $1
set ffile = "./grid{$grid}.dat"


if ( ! -e $pwd/${outdir}-dat) then
   mkdir $pwd/${outdir}-dat
endif

if ( ! -e $pwd/${outdir}-out) then
   mkdir $pwd/${outdir}-out
endif


#./read_geom.sh 1

set fname = CH3-uccsdt-tz.${irun}.$outdir
#set xyzfile = ch4.${irun}.xyz3.dat
#set datfile = ch4.${irun}.d3.dat

echo $fname
#
cat<<endb> $fname.com

memory,1000,m
gthresh,energy=1.d-10,zero=1.d-14,thrint=1.d-14,oneint=1.d-14,twoint=1.d-14,prefac=1.d-20


PROC CCT
{hf,ipnit=1,maxdis=30;wf,9,1,1;maxit,100;accu,16}
{uccsd(t),thrden=1.0d-08}
ENDPROC

h= 0.002d0

PROC DCCT
 field=nn*h
 dip,field(1),field(2),field(3)
 {hf,ipnit=1,maxdis=30;wf,9,1,1;maxit,100;accu,17}
 {uccsd(t),thrden=1.0d-09}
ENDPROC

basis=aug-cc-pVTZ

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

i1 = 1

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
  
  !
  nn = [ 0, 0, 0];DCCT;e000  = energy
  !
  nn = [ 2, 0, 0];DCCT;ep200 = energy
  nn = [-2, 0, 0];DCCT;em200 = energy
  !
  nn = [ 0, 2, 0];DCCT;e0p20 = energy
  nn = [ 0,-2, 0];DCCT;e0m20 = energy
  !
  nn = [ 0, 0, 2];DCCT;e00p2 = energy
  nn = [ 0, 0,-2];DCCT;e00m2 = energy
  !
  nn = [ 1, 1, 0];DCCT;ep1p10 = energy
  nn = [-1, 1, 0];DCCT;em1p10 = energy
  nn = [ 1,-1, 0];DCCT;ep1m10 = energy
  nn = [-1,-1, 0];DCCT;em1m10 = energy
  !
  nn = [ 0, 1, 1];DCCT;e0p1p1 = energy
  nn = [ 0,-1, 1];DCCT;e0m1p1 = energy
  nn = [ 0, 1,-1];DCCT;e0p1m1 = energy
  nn = [ 0,-1,-1];DCCT;e0m1m1 = energy
  !
  nn = [ 1, 0, 1];DCCT;ep10p1 = energy
  nn = [-1, 0, 1];DCCT;em10p1 = energy
  nn = [ 1, 0,-1];DCCT;ep10m1 = energy
  nn = [-1, 0,-1];DCCT;em10m1 = energy
  !
  alpha_xx(i1) = 0.25d0*( ep200+em200-2.d0*e000 )/( h*h )
  alpha_yy(i1) = 0.25d0*( e0p20+e0m20-2.d0*e000 )/( h*h )
  alpha_zz(i1) = 0.25d0*( e00p2+e00m2-2.d0*e000 )/( h*h )
  !
  alpha_xy(i1) = 0.25d0*( ep1p10-em1p10-ep1m10+em1m10 )/(h*h)
  alpha_xz(i1) = 0.25d0*( ep10p1-em10p1-ep10m1+em10m1 )/(h*h)
  alpha_yz(i1) = 0.25d0*( e0p1p1-e0m1p1-e0p1m1+e0m1m1 )/(h*h)
  !
  dipx(i1) = 0.25d0*( ep200-em200 )/h
  dipy(i1) = 0.25d0*( e0p20-e0m20 )/h
  dipz(i1) = 0.25d0*( e00p2-e00m2 )/h

  eeee(i1)=energy
  !
  NNZ(i1) = 'ppp'
  IZ(i1) = ${irun}
  !
  text ### CH3 AVTZ
  table,NNZ,IZ,r1,r2,r3,r4,r5,r6,e000,dipx,dipy,dipz,alpha_xx,alpha_yy,alpha_zz,alpha_xy,alpha_xz,alpha_yz
  DIGITS,0,  0, 8, 8, 8, 8, 8, 8,  12,  12,  12,  12,      12,      12,      12,      12,      12,      12
  !
  save,${irun}.dat,new
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

if (-e $irun.dat) then
   /bin/cp  $irun.dat  $pwd/${outdir}-dat
   /bin/cat  $irun.dat >>  $pwd/${outdir}-dat/all.dat
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

