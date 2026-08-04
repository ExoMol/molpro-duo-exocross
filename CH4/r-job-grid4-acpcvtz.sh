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

set outdir = f12b-g4-avctz

set basis = aug-cc-pCVTZ
set grid = 4
set core = ""


set irun = $1
set ffile = "./inp/ch4_grid{$grid}.dat"


if ( ! -e $pwd/${outdir}-dat) then
   mkdir $pwd/${outdir}-dat
endif

if ( ! -e $pwd/${outdir}-out) then
   mkdir $pwd/${outdir}-out
endif


#./read_geom.sh 1

set fname = ch4.${irun}.$outdir
#set xyzfile = ch4.${irun}.xyz3.dat
#set datfile = ch4.${irun}.d3.dat

echo $fname
#
cat<<endb> $fname.com

gthresh,twoint=1.d-10,prefac=1.d-18,zero=1.d-10,energy=1.d-10;
memory,512,m

PROC CCT-f12
 {hf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,18;wf,orbital}
 !{ccsd(t)-f12c;core;thresh,energy=1.d-10,coeff=1.d-9,thrint=1.d-10,zero=1.d-10;maxit,100}
 !{ccsd(t)-f12c,thrden=1.0d-12,ri_basis=ri_basis_set,df_basis=AV5Z/MP2FIT,df_basis_exch=V5Z/JKFIT,gem_beta=1.1; thresh,energy=1.d-9,coeff=1.d-9,thrint=1.d-10,zero=1.d-10;maxit,100}
 {ccsd(t)-f12b,df_basis=aug-cc-pwcv5z/mp2fit,df_basis_exch=aug-cc-pv5z/jkfit,ri_basis=cc-pCVTZ-f12/optri; thresh,energy=1.d-9,coeff=1.d-9,thrint=1.d-10,zero=1.d-10;maxit,100}
 !{ccsd(t)-f12b,thrden=1.0d-10,ri_basis=ri_basis_set,df_basis=AWCV5Z/MP2FIT,df_basis_exch=V5Z/JKFIT,gem_beta=1.1;thresh,energy=1.d-9,coeff=1.d-9,thrint=1.d-10,zero=1.d-10;maxit,100}
 !!{OPTG,grad=5,energy=8,saveact='xxx.opt';print,history;inactive,alpha}
 !{readvar,'xxx.opt'}
ENDPROC

PROC calc-mp2
 {hf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,18;wf,orbital}
 {mp2}
ENDPROC


ANSATZ=3c,fix=1,canonical=1

! cc-pVQZ-F12

basis={
!
C=aug-cc-pCVTZ
H=aug-cc-pCVTZ
!
}



!dkroll=1

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

 xxen = xEn
 !
 irunmode = mod($irun,5)
 !
 if (xEn.lt.8000.or.mod($irun,5).eq.0) then
  !
  !
  gexpec,rel,darwin,massv
  !
  CCT-f12
  !
  eeee(i1)=energy
  !
  relen(i1) =erel
  !
  NNZ(i1) = 'ppp'
  IZ(i1) = ${irun}
  !
  table,NNZ,IZ,R1CH,R2CH,R3CH,R4CH,A12,A13,A24,A23,A14,eeee,relen
  DIGITS,0,  0,   8,   8,   8,   8,  8,  8,  8,  8,  8,  12,   12
  !
  save,$fname.dat,new
  !
  !title PES of CH4
  !
 endif
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



module load mkl/10.2.5/035
module load compilers/intel/11.1/072
module load mpi/openmpi/1.4.1/intel

setenv molproexe /shared/ucl/apps/Molpro/2012.1-bindist/bin/

#setenv molproexe /shared/ucl/apps/Molpro/2009.1/bin/
#setenv PATH /shared/ucl/apps/Molpro/2009.1/bin/:"${PATH}"


echo "Running molpro -n 1 -d $TMPDIR -W $wdir < $pwd/$fname.com > $pwd/$fname.out"

time $molproexe/molpro -n 1 -I $TMPDIR -d $TMPDIR -W $wdir < $pwd/${fname}.com > $pwd/${fname}.out

#$molproexe/molpro -n 1 -I $TMPDIR -d $TMPDIR -W $wdir < $pwd/${fname}.com > $wdir/${fname}.out

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

