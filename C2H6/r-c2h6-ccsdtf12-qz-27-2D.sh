#!/bin/csh
#
# Generation of the input file for MOLPRO
#


cd $2
set pwd = `pwd`
echo $pwd


#grid file number
set grid = 3
set job = 27

echo "2"

set outdir = dqz-$job-$grid

set irun = $1

set ffile = "./grids/grid-2d-v3.dat"

echo "3"




if ( ! -e $pwd/${outdir}-dat) then
   mkdir $pwd/${outdir}-dat
endif

if ( ! -e $pwd/${outdir}-out) then
   mkdir $pwd/${outdir}-out
endif

echo "4"


#./read_geom.sh 1

set fname = c2h6.${irun}.$outdir.$job
#set xyzfile = ch4.${irun}.xyz3.dat
#set datfile = ch4.${irun}.d3.dat



echo "5"


echo $fname
#
cat<<endb> $fname.com

memory,5000,m
gthresh,energy=1.d-10,zero=1.d-14,thrint=1.d-14,oneint=1.d-14,twoint=1.d-14,prefac=1.d-20

PROC CCT-f12
!h= 0.002d0
!field=nn*h
!dip,field(1),field(2),field(3)
{rhf,ipnit=1,maxdis=30;wf,18,1,0;maxit,150;accu,18}
!
{ccsd(t)-f12b,thrden=1.0d-10,ri_basis=optri,df_basis=awcv5z/mp2fit,df_basis_exch=av5z/jkfit,gem_beta=1.0,MAXIT=300}
!
!{mp2}
!put,molden,$fname.molden
ENDPROC

PROC calc-mp2
 {hf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,18;wf,orbital}
 {mp2}
ENDPROC



symmetry,nosym
geometry={angstrom
c
c, 1, rcc
h, 1, rch1, 2, acch1
h, 2, rch4, 1, acch4, 3, ahh1
h, 1, rch2, 2, acch2, 4, ahh2
h, 2, rch5, 1, acch5, 5, ahh3
h, 1, rch3, 2, acch3, 6, ahh4
h, 2, rch6, 1, acch6, 7, ahh5}


!basis = Vqz
basis= vqz-f12
ansatz=3c,fix=1,canonical=1


i1 = 1
field=[-0.005, 0.005]

`awk -f print_geom.awk irow=$irun $ffile`

ZPE = -79.70940816
! Optimized variables
!
emp20    = -79.660149520022
r1e      =    1.52521635
r2e      =    1.09050281
alpha1e  =    111.20495396


rcc    =xR1+  r1e     
rch1  = xR2+  r2e     
rch2  = xR3+  r2e     
rch3  = xR4+  r2e     
rch4  = xR5+  r2e     
rch5  = xR6+  r2e     
rch6  = xR7+  r2e     

acch1 = xA1+  alpha1e 
acch2 = xA2+  alpha1e 
acch3 = xA3+  alpha1e 

acch4 = xA4+  alpha1e 
acch5 = xA5+  alpha1e 
acch6 = xA6+  alpha1e 

s14   = xA7+  0.0
s15   = xA8+  0.0
s16   = xA9+  0.0
s17   = xA10+ 0.0
s18   = xA11+ 0.0

pi = 180.0

tau14 = sqrt(2.0)*S17/3.0-sqrt(2.0)*S15/3.0+S18
tau24 = sqrt(2.0)*S17/3.0+sqrt(2.0)*S15/6.0+sqrt(6.0)*S14/6.0+S18-2.0/3.0*pi
tau25 = -sqrt(6.0)*S16/6.0+sqrt(6.0)*S14/6.0-sqrt(2.0)*S17/6.0+sqrt(2.0)*S15/6.0+S18
tau35 = -sqrt(6.0)*S16/6.0-sqrt(6.0)*S14/6.0-sqrt(2.0)*S17/6.0+sqrt(2.0)*S15/6.0+S18-2.0/3.0*pi
tau36 = -sqrt(2.0)*S17/6.0+sqrt(2.0)*S15/6.0+sqrt(6.0)*S16/6.0-sqrt(6.0)*S14/6.0+S18

ahh1 =tau14
ahh2 =tau24
ahh3 =tau25
ahh4 =tau35
ahh5 =tau36


NNZ(i1) = 'ppp'
IZ(i1) = ${irun}
!
r1(i1) =  rcc  
r2(i1) =  rch1 
r3(i1) =  rch2 
r4(i1) =  rch3 
r5(i1) =  rch4 
r6(i1) =  rch5 
r7(i1) =  rch6 
r8(i1) =  acch1
r9(i1) =  acch2
r10(i1) = acch3
r11(i1) = acch4
r12(i1) = acch5
r13(i1) = acch6
r14(i1) = ahh1 
r15(i1) = ahh2 
r16(i1) = ahh3 
r17(i1) = ahh4 
r18(i1) = ahh5 


q1(i1) =  xR1 
q2(i1) =  xR2 
q3(i1) =  xR3 
q4(i1) =  xR4 
q5(i1) =  xR5 
q6(i1) =  xR6 
q7(i1) =  xR7 
q8(i1) =  xA1 
q9(i1) =  xA2 
q10(i1) = xA3 
q11(i1) = xA4 
q12(i1) = xA5 
q13(i1) = xA6 
q14(i1) = xA7 
q15(i1) = xA8 
q16(i1) = xA9 
q17(i1) = xA10
q18(i1) = xA11





  !
  !CCT-f12
  !
  calc-mp2
  !
  !eeee(i1)=energy
  !
  NNZ(i1) = 'ppp'
  IZ(i1) = ${irun}


  !calc-MP2
  !
  emp2cm = (energy - emp20)*tocm
  !
  !eeee(i1) = energy
  !
  !gexpec,rel,darwin,massv;
  !
  !rele(i1) = erel
  !
  if (emp2cm.lt.60000.0) then
       !  
       nn = [ 0, 0, 0];CCT-f12;e000  = energy
       eeee(i1)=(energy-ZPE)*tocm
       !
       !nn = [ 1, 0, 0];CCT-f12;ep200 = energy
       !nn = [-1, 0, 0];CCT-f12;em200 = energy
       !
       !nn = [ 0, 1, 0];CCT-f12;e0p20 = energy
       !nn = [ 0,-1, 0];CCT-f12;e0m20 = energy
       !
       !nn = [ 0, 0, 1];CCT-f12;e00p2 = energy
       !nn = [ 0, 0,-1];CCT-f12;e00m2 = energy
       !
       !dipx(i1) = 0.5d0*( ep200-em200 )/h*todebye;
       !dipy(i1) = 0.5d0*( e0p20-e0m20 )/h*todebye;
       !dipz(i1) = 0.5d0*( e00p2-e00m2 )/h*todebye;
     
       !
       NNZ(i1) = 'ppp'
       IZ(i1) = ${irun}
       !
       text ### C2H6
       !table,NNZ,IZ,r1,r2,r3,r4,r5,r6,eeee,dipx,dipy,dipz
       !DIGITS,0,  0, 8, 8, 8, 8, 8, 8,  12,  12,  12,  12
       !
       !table,NNZ,IZ,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16,q17,q18,eeee,dipx,dipy,dipz
       !DIGITS,0,  0, 8, 8, 8, 8, 8, 8  8, 8, 8,  8  ,8,  8,  8,  8,  8,  8,  8,  8, 8, 8, 8, 8, 8, 8, 8, 8, 8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  12,  12,  12,  12
       
       ww(i1) = 1.0       


       !table ,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,eeee, ww,NNZ,IZ,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16,q17,q18 !,dipx,dipy,dipz
       !DIGITS, 8, 8, 8, 8, 8, 8  8, 8, 8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  12,  2,  4, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8,  8,  8,  8,  8,  8,  8,  8,  8,  8 !,  12,  12,  12


       table ,NNZ,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,eeee, ww,IZ,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16,q17,q18
       DIGITS,  0, 8, 8, 8, 8, 8, 8  8, 8, 8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  12,  2, 0, 8, 8, 8, 8, 8, 8, 8, 8, 8,  8,  8,  8,  8,  8,  8,  8,  8,  8


       !
       save,$fname.dat,new
       !
       SHOW,energy;
       !
       !
       put, xyz, $fname.x.dat,new
       !
       !goto,END:
       !
  endif




  !
  !
---

endb


setenv OMP_NUM_THREADS 1

limit
#echo $TMP

setenv wdir $TMPDIR

cp $pwd/$fname.com $wdir

# Run Molpro
echo "System TMPDIR = $TMPDIR"
echo "wdir = $wdir"

cd $wdir


module load molpro/2020.1/openmp


echo "Running molpro -n 1 -d $wdir -W $wdir < $pwd/$fname.com > $pwd/$fname.out"

time molpro -n 1 -I $wdir -d $wdir -W $wdir < ${fname}.com > ${fname}.out

#cp $wdir/$fname.out $pwd


#$molproexe/molpro -n 1 -I $TMPDIR -d $TMPDIR -W $wdir < $pwd/${fname}.com > $wdir/${fname}.out

if (-e $wdir/$fname.dat) then
   /bin/cp  $wdir/$fname.dat  $pwd/${outdir}-dat
endif

if (-e $wdir/$fname.x.dat) then
   /bin/cp  $wdir/$fname.x.dat  $pwd/${outdir}-dat
endif

if (-e $wdir/$fname.out) then
    gzip $wdir/$fname.out
   /bin/mv  $wdir/$fname.out.gz  $pwd/${outdir}-out
endif


cd $pwd


if (-e $pwd/$fname.com) then
   /bin/rm  $pwd/$fname.com
endif




