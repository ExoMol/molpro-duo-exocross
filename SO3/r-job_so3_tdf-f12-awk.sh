#!/bin/csh 
#
# Generation of the input file for MOLPRO 
#

cd $2
set pwd = `pwd`
echo $pwd

set irun = $1



set fname = so3.${irun}.tdf
#
cat<<endb> $fname.com

gthresh,energy=1.d-10,zero=1.d-14,thrint=1.d-14,oneint=1.d-14,twoint=1.d-14,prefac=1.d-20



memory,1600,m
PROC CCT_opt
 {hf;orbprint,-1;maxit,100;accu,18;wf,charge=1}
 {ccsd(t); thresh,energy=1.d-9,coeff=1.d-9,thrint=1.d-10,zero=1.d-10;maxit,100}
 {STATUS,nocheck}
 {OPTG,grad=5,energy=8,saveact='xxx.opt';print,history;inactive,alpha}
 {readvar,'xxx.opt'}
!{frequencies;}
ENDPROC

PROC CCT
 {hf;orbprint,-1;maxit,100;accu,18;wf,charge=1}
 {ccsd(t); thresh,energy=1.d-9,coeff=1.d-9,thrint=1.d-10,zero=1.d-10;maxit,100}
 {STATUS,nocheck}
ENDPROC


PROC HF-proc
 {rhf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,17;wf,40,1,orbital,IGNORE_ERROR}
ENDPROC


PROC CCT-f12
 !{rhf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,17;wf,40,1,orbital,IGNORE_ERROR}
 {ccsd(t)-f12b,thrden=1.0d-10,ri_basis=ri_basis_set,df_basis=AWCV5Z/MP2FIT,df_basis_exch=V5Z/JKFIT,gem_beta=1.2;maxit,100}
 {STATUS,nocheck}
 !{OPTG,grad=5,energy=8,saveact='xxx.opt';print,history;inactive,alpha}
 !{readvar,'xxx.opt'}
ENDPROC

ANSATZ=3C,fix=1,hybrid=1,canonical=1



BASIS
!
O=aug-cc-pVTZ,S=aug-cc-pV(T+d)Z
!
set,ri_basis_set
!
!
! K.E. Yousaf and K.A. Peterson, Chem. Phys. Lett. 476, 303 (2009).
! Auxiliary RI (OptRI) matched to the aug-cc-pV(T+d)Z Basis Set for O
!
s,O,1.135258E+01,3.375893E+00,4.533015E-01,1.104366E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;
!
p,O,5.456091E+00,2.352584E+00,1.067952E+00,3.301448E-01,8.936806E-02;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;
c,5.5,1.000000E+00;
!
d,O,1.206985E+01,4.278700E+00,1.374781E+00,3.193272E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;
!
f,O,6.915251E+00,3.242493E+00,7.499346E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
!
g,O,3.141319E+00,1.075471E+00;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
!
! K.E. Yousaf and K.A. Peterson, Chem. Phys. Lett. 476, 303 (2009).
! Auxiliary RI (OptRI) matched to the aug-cc-pV(T+d)Z Basis Set for S
!
s,S,7.154448E+00,1.705679E+00,3.879646E-01,7.994920E-02;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;
!
p,S,1.234307E+01,3.496947E+00,1.630397E+00,5.394127E-01,5.646997E-02;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;
c,5.5,1.000000E+00;
!
d,S,1.929014E+01,7.371754E+00,1.876286E+00,4.365223E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;
!
f,S,5.439897E+00,8.572979E-01,3.723605E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
!
g,S,1.051297E+00,4.871735E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
END



symmetry,nosym;
orient,noorient;

geometry={angstrom;
 S;
 O, 1, R1SO;
 O, 1, R2SO, 2, AA3;
 O, 1, R3SO, 2, AA2, 3, AA1, 1;
}

!BASIS=VDZ-F12
!ANSATZ=3c,fix=1,canonical=1
!orbital,2102.2,IGNORE_ERROR


hartree = 219474.63
emp20 = -623.083175267015

!stretching
jt0 = 1
gRR=[1.43, 1.44, 1.42, 1.46, 1.40, 1.35, 1.60, 1.30, 1.70, 1.10, 2.00]
!bending
gAA=[120., 119., 121.0, 117.0, 123., 115., 110.0, 100.0, 90.0, 70.0 ]


Vsmax= #gRR
Vbmax= #gAA


imin =     ${irun}
imax =     ${irun} 

i=0
i1 = 0
n1 = 0
do jt1=jt0,Vsmax
 do jt2=jt1,Vsmax
 do jt3=jt2,Vsmax
   do k1=1,Vbmax 
   do k2=1,k1
   do k3=1,k2
      !
      i=i+1
      !
       if (jt1.eq.jt2.and.jt1.eq.jt3.and.k1+k2+k3.le.24) then
        irunjob = 1;
       else if (k1.eq.k2.and.k1.eq.k3.and.jt1+jt2+jt3.le.24) then
        irunjob = 2;
       else if ((jt1+jt2+jt3).le.10) then
        irunjob = 3;
       else if ((k1+k2+k3).le.10) then
        irunjob = 4;
       else if (abs(jt1+jt2+jt3+k1+k2+k3).le.15) then
        irunjob = 5;
       else if (mod(i,8).eq.0) then
        irunjob = 10;
       else
        irunjob = 0;
       end if;
      !
      SHOW,irunjob;
      !
      if (irunjob.ne.0) then
         !
         i1 = i1+1
         if (imin.le.i1.and.i1.le.imax) then
            !
            n1 = n1 + 1 
            !
            SHOW,n1;
            !
            R1SO = gRR(jt1)
            R2SO = gRR(jt2)
            R3SO = gRR(jt3)
            AA1  = gAA(k1)
            AA2  = gAA(k2)
            AA3  = gAA(k3)
            !
            R1SOt(n1) =gRR(jt1)
            R2SOt(n1) =gRR(jt2)
            R3SOt(n1) =gRR(jt3)
            AAt1(n1)  = gAA(k1)
            AAt2(n1)  = gAA(k2)
            AAt3(n1)  = gAA(k3)
            !
            gexpec,rel,darwin,massv;
            !
            HF-proc
            !
            MP2
            !
            emp2cm = (energy - emp20)*tocm
            !
            if (emp2cm.lt.50000.0) then
               !
               CCT-f12
               eeee(n1)=energy
               !
               show massv*,darw,erel
               massvccsdt(n1)=massv
               darwccsdt(n1)=darwin
               erelen(n1)=erel
               !
               SHOW,energy;
               !
               xxx(n1) = 'XXX'
               yyy(n1) = '111'
               !
               text ### SO3 energies 
               table,xxx,R1SOt,R2SOt,R3SOt,AAt1,AAt2,AAt3,eeee,massvccsdt,darwccsdt,erelen
               DIGITS, 3,    5,    5,    5,   5,   5,   5,  12,        12,       12,    12 
               !
            endif
            !
         end if
      end if
   end do
   end do
   end do
 end do
 end do
end do

    
text ### SO3 nergies 
table,yyy,R1SOt,R2SOt,R3SOt,AAt1,AAt2,AAt3,eeee,massvccsdt,darwccsdt,erelen
DIGITS, 3,    5,    5,    5,   5,   5,   5,  12,        12,       12,    12 

save,so3.tdf.${irun}.dat,new

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

$molproexe/molpro -n 1 -I $TMPDIR -d $TMPDIR -W $wdir < $pwd/${fname}.com > $pwd/${fname}.out

if (-e so3*.dat) then
   /bin/cp  *.dat  $pwd
endif

