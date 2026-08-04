#!/bin/bash -l 
#
# Generation of the input file for MOLPRO 
#

cd $2
export pwd=`pwd`
echo $pwd

export irun=$1
export ijob=$3

export outdir=q.c-12.${ijob}

echo $outdir

#5z-7220_noSO_THRNRM_G2_${ijob}


if [ ! -e $pwd/${outdir}-dat ]; then
   mkdir $pwd/${outdir}-dat
fi

if [ ! -e $pwd/${outdir}-out ]; then
   mkdir $pwd/${outdir}-out
fi


export fname=P2H2.${irun}.${outdir}

echo $fname
#
cat<<endb> $fname.com

emp20 =  -679.1914

gthresh,twoint=1.d-10,prefac=1.d-18,zero=1.d-10,energy=1.d-10;
memory,512,m

PROC calc-mp2
 {hf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,18;wf,orbital}
 {mp2}
ENDPROC


PROC CCT_opt
 {hf;orbprint,-1;maxit,100;accu,18}
 {ccsd(t); thresh,energy=1.d-9,coeff=1.d-9,thrint=1.d-10,zero=1.d-10;maxit,100}
 {STATUS,nocheck}
 {OPTG,grad=5,energy=8,saveact='xxx.opt';print,history;inactive,alpha}
 {readvar,'xxx.opt'}
!{frequencies;}
ENDPROC



PROC CCT-f12
 {rhf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,18;}
 {ccsd(t)-f12b,thrden=1.0d-12,ri_basis=ri_basis_set,df_basis=AWCV5Z/MP2FIT,df_basis_exch=V5Z/JKFIT,gem_beta=1.2;
 thresh,energy=1.d-10,coeff=1.d-10,thrint=1.d-10,zero=1.d-11;maxit,100}
 !{STATUS,nocheck}
 !{OPTG,grad=5,energy=8,saveact='xxx.opt';print,history;inactive,alpha}
 !{readvar,'xxx.opt'}
ENDPROC

ANSATZ=3C,fix=1,hybrid=1,canonical=1

! cc-pCVQZ-F12_MP2FIT

basis={

s,P,6.355220E+02,3.923920E+02,2.932070E+02,2.114600E+02,1.453450E+02,8.420470E+01,3.478690E+01,2.152780E+01,1.324260E+01,7.093080E+00,4.176580E+00,2.341360E+00,1.222360E+00,7.405470E-01,4.947150E-01,2.827430E-01,1.653590E-01,1.141550E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;
c,5.5,1.000000E+00;
c,6.6,1.000000E+00;
c,7.7,1.000000E+00;
c,8.8,1.000000E+00;
c,9.9,1.000000E+00;
c,10.10,1.000000E+00;
c,11.11,1.000000E+00;
c,12.12,1.000000E+00;
c,13.13,1.000000E+00;
c,14.14,1.000000E+00;
c,15.15,1.000000E+00;
c,16.16,1.000000E+00;
c,17.17,1.000000E+00;
c,18.18,1.000000E+00;

p,P,6.634800E+02,4.137660E+02,2.585700E+02,1.542730E+02,9.265490E+01,3.742510E+01,2.215460E+01,1.075400E+01,5.008000E+00,2.571690E+00,1.357630E+00,8.164480E-01,5.440010E-01,3.062160E-01,1.702200E-01,1.192660E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;
c,5.5,1.000000E+00;
c,6.6,1.000000E+00;
c,7.7,1.000000E+00;
c,8.8,1.000000E+00;
c,9.9,1.000000E+00;
c,10.10,1.000000E+00;
c,11.11,1.000000E+00;
c,12.12,1.000000E+00;
c,13.13,1.000000E+00;
c,14.14,1.000000E+00;
c,15.15,1.000000E+00;
c,16.16,1.000000E+00;

d,P,1.842510E+02,1.162530E+02,6.761150E+01,3.799750E+01,1.785460E+01,1.150510E+01,5.481380E+00,3.191730E+00,1.796510E+00,1.053960E+00,6.886240E-01,3.616740E-01,2.134560E-01,1.295120E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;
c,5.5,1.000000E+00;
c,6.6,1.000000E+00;
c,7.7,1.000000E+00;
c,8.8,1.000000E+00;
c,9.9,1.000000E+00;
c,10.10,1.000000E+00;
c,11.11,1.000000E+00;
c,12.12,1.000000E+00;
c,13.13,1.000000E+00;
c,14.14,1.000000E+00;

f,P,8.503850E+01,4.671720E+01,2.912320E+01,1.395010E+01,6.281070E+00,2.914270E+00,1.710520E+00,7.429880E-01,4.843800E-01,2.883870E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;
c,5.5,1.000000E+00;
c,6.6,1.000000E+00;
c,7.7,1.000000E+00;
c,8.8,1.000000E+00;
c,9.9,1.000000E+00;
c,10.10,1.000000E+00;

g,P,2.810970E+01,1.415300E+01,4.540640E+00,1.896630E+00,8.976870E-01,5.950910E-01,3.336730E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;
c,5.5,1.000000E+00;
c,6.6,1.000000E+00;
c,7.7,1.000000E+00;

h,P,2.728210E+00,9.319940E-01,4.114320E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;

i,P,1.594580E+00;
c,1.1,1.000000E+00;


! cc-pVQZ-F12_MP2FIT

s,H,9.713130E+01,1.893660E+01,5.433230E+00,1.938230E+00,1.273240E+00,6.862960E-01,3.344500E-01,1.429990E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;
c,5.5,1.000000E+00;
c,6.6,1.000000E+00;
c,7.7,1.000000E+00;
c,8.8,1.000000E+00;

p,H,1.079190E+01,3.943920E+00,1.737290E+00,7.238210E-01,3.579070E-01,1.760460E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;
c,5.5,1.000000E+00;
c,6.6,1.000000E+00;

d,H,3.831630E+00,1.548550E+00,6.693950E-01,4.279490E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;

f,H,3.880630E+00,1.538820E+00,7.066550E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;

g,H,3.071320E+00,8.983620E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;

h,H,8.074910E-01;
c,1.1,1.000000E+00;

}


symmetry,nosym;
orient,noorient;

geometry={angstrom;
  P;
  P, 1, OO;
  H, 1, OH1, 2, HOO1;
  H, 2, OH2, 1, HOO2, 3, HOOH;
}


gridr1  = [2.00,2.02,1.98,2.05,1.96,2.10,1.92,2.15,1.90,2.20,1.85,2.50,1.80,2.80,1.70,3.00,4.00]
gridr2  = [1.40,1.42,1.38,1.44,1.36,1.50,1.32,1.55,1.25,1.60,1.20,1.70,1.10,1.80,1.10,2.00,1.00,3.00,4.00] 
grida1  = [90.0,92.0,88.0,96.0,86.0,100.0,82.0,105.0,78.0,110.,74.0,115.0,70.0,120.0,65.0]
#grida2  = [90.0,110.,120.0,150.,180.0] 
grida2  = [0.0,30.0,60.0,90.0,110.,120.0,150.,180.0] 

i1 = 1
field=[-0.005, 0.005]

N1max= 1 !#gridr1
N2max= 1 !#gridr2
N3max= 1 !#grida1
N4max= 1 !#grida2
!
!N1max = 1 ; N2max = 1; N3max = 1; N4max = 1;
!
!
imin =     ${irun}
imax =     ${irun} 
!
i  = 0
i1 = 0
n1 = 0

yyy(1)='0';xOO(1)=0;xOH1(1)=0;xOH2(1)=0;xHOO1(1)=0;xHOO2(1)=0;xHOOH(1)=0;eeee(1)=0;dip6dx(1)=0;dip6dy(1)=0;dip6dz(1)=0

NOGPRINT,VARIABLE

do k3=1,N4max
  do q1=1,N1max
   do q2=1,N2max
    do q3=1,N2max
     if ( gridr2(q2).le.gridr2(q3) ) then
     !
     show,q1,q2,q3
     !
     do k1=1,N3max
      do k2=1,N3max
       if ( grida1(k1).le.grida1(k2) ) then 
       show,k1,k2
       !
       i=i+1
       !
       if (q2.eq.q3.and.q1+k1+k2.le.5) then
        irunjob = 1;
       else if (k1.eq.k2.and.q1+q2+q3.le.6) then
        irunjob = 2;
       else if (q2.eq.q3.and.k1.eq.k2.and.q1+q2+q3+k1+k2.le.18) then
        irunjob = 3;
       else if (q1+q2+q3.le.10.and.k1+k2.le.7) then
        irunjob = 4;
       else if (q1+(q2+q3)*0.9+k1+k2.le.18) then
        irunjob = 6;
       else if (mod(i,200).eq.0) then	
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
            GPRINT,VARIABLE
            !
            n1 = n1 + 1 
            !
            SHOW,n1;
            !
            OO  = gridr1(q1)
            OH1  = gridr2(q2)
            OH2  = gridr2(q3)
            HOO1 = grida1(k1)
            HOO2 = grida1(k2)
            HOOH = grida2(k3)
            !
            xOO(n1) = gridr1(q1)
            xOH1(n1) = gridr2(q2)
            xOH2(n1) = gridr2(q3)
            xHOO1(n1) = grida1(k1)
            xHOO2(n1) = grida1(k2)
            xHOOH(n1) = grida2(k3)
            !
            !gexpec,rel,darwin,massv;
            !
            calc-mp2
            !
            !MP2
            !
            emp2cm = (energy - emp20)*tocm
            !
            if (emp2cm.lt.5000000000.0) then
                !
                CCT-f12
                eeee(n1) = energy
                !
                SHOW,energy;
                !
                xxx(n1) = '222'
                yyy(n1) = '111'
                !
                text ### HOOH dipole
                table,xxx,xOO,xOH1,xOH2,xHOO1,xHOO2,xHOOH,eeee
                DIGITS, 3,  5,   5,   5,   5,     5,    5,  12
                !
                NOGPRINT,VARIABLE
                !
                goto,END:
                !
             endif
          end if
        end  if
       end if 
      end do 
     end do 
    end if
   end do   
  end do   
 end do    
end do 

END: TEXT ### HOOH


text ### HOOH dipole
table,yyy,xOO,xOH1,xOH2,xHOO1,xHOO2,xHOOH,eeee
DIGITS, 3,  5,   5,   5,   5,     5,    5,  12
          
save,$fname.dat,new


---

endb

export OMP_NUM_THREADS=1

#limit
#echo $TMP

export TMP=$TMPDIR
export wdir=$TMPDIR


# Run Molpro
echo "System TMPDIR = $TMPDIR"
echo "wdir = $wdir"

cd $wdir

module load gcc-libs/4.9.2
module unload compilers/intel/2015/update2
module load compilers/gnu/4.9.2
module unload mpi/intel/2015/update3/intel
module load mpi/openmpi/1.8.4/gnu-4.9.2
module load molpro/2012.1.25/gnu-4.9.2


echo "Running molpro -n 1 -d $TMPDIR -W $wdir < $pwd/$fname.com > $pwd/$fname.out"

time molpro -n 1 -I $TMPDIR -d $TMPDIR -W $wdir < $pwd/${fname}.com > $pwd/${fname}.out

#$molproexe/molpro -n 1 -I $TMPDIR -d $TMPDIR -W $wdir < $pwd/${fname}.com > $wdir/${fname}.out



if [ -e $fname.dat ]; then
   /bin/cp  $fname.dat  $pwd/${outdir}-dat
fi

if [ -e $fname.xyz.dat ]; then
   /bin/cp  $fname.xyz.dat  $pwd/${outdir}-dat
fi

if [ -e $pwd/$fname.out ]; then
    gzip $pwd/$fname.out
   /bin/mv  $pwd/$fname.out.gz  $pwd/${outdir}-out
fi


if (-e $pwd/$fname.com) then
   /bin/rm  $pwd/$fname.com
endif


