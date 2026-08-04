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

set irun = $1


set fname = ch4.${irun}.tdz.ccsd.f12

echo $fname
#
cat<<endb> $fname.com

emp20 =  -151.320212816163

gthresh,twoint=1.d-10,prefac=1.d-18,zero=1.d-10,energy=1.d-10;
memory,512,m

PROC calc-mp2
 {hf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,18;wf,orbital}
 {mp2}
ENDPROC


PROC CCT-f12
 {rhf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,17;wf,orbital}
 {ccsd(t)-f12b,thrden=1.0d-10,ri_basis=ri_basis_set,df_basis=AWCV5Z/MP2FIT,df_basis_exch=V5Z/JKFIT,gem_beta=1.2;maxit,100}
 !{mp2}
 !{STATUS,nocheck}
 !{OPTG,grad=5,energy=8,saveact='xxx.opt';print,history;inactive,alpha}
 !{readvar,'xxx.opt'}
ENDPROC
!
ANSATZ=3C,fix=1,hybrid=1,canonical=1
!
BASIS={
!
O=aug-cc-pV(T+d)Z,H=aug-cc-pV(T+d)Z
!
set,ri_basis_set
!
!Auxiliary RI (OptRI) matched to the aug-cc-pV(n+d)Z Quadruple-zeta. Basis Set for O and H
!This should be used within the CABS framework with the above orbital set
!
! aug-cc-pVTZ
!
s,H,1.747055E+00,4.888246E-01,2.076744E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
!
p,H,4.837096E+00,2.114191E+00,7.026262E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
!
d,H,1.784914E+00,6.175128E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
!
f,H,1.918725E+00;
c,1.1,1.000000E+00;
!
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
}



symmetry,nosym;
orient,noorient;

geometry={angstrom;
  O;
  O, 1, OO;
  H, 1, OH1, 2, HOO1;
  H, 2, OH2, 1, HOO2, 3, HOOH;
}


gridr1  = [1.46,1.47,1.45,1.49,1.43,1.52,1.40,1.55,1.35,1.60,1.30,1.70,1.20,1.80,1.10,2.00,3.00]
gridr2  = [0.97,0.98,0.96,1.00,0.94,1.10,0.90,1.20,0.85,1.25,0.80,1.30,0.70,1.40,0.60,1.50,1.60,1.80,2.00] 
grida1  = [90.3,92.0,88.0,96.0,86.0,100.0,82.0,105.0,78.0,110.,74.0,115.0,70.0,120.0,65.0]
grida2  = [0.0,30.0,60.0,90.0,110.,120.0,150.,180.0] 

i1 = 1
field=[-0.005, 0.005]

N1max= 1 !#gridr1
N2max= 1 !#gridr2
N3max= 1 !#grida1
N4max= #grida2
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
       if (q2.eq.q3.and.q1+k1+k2.le.6) then
        irunjob = 1;
       else if (k1.eq.k2.and.q1+q2+q3.le.8) then
        irunjob = 2;
       else if (q2.eq.q3.and.k1.eq.k2.and.q1+q2+q3+k1+k2.le.24) then
        irunjob = 3;
       else if (q1+q2+q3.le.10.and.k1+k2.le.8) then
        irunjob = 4;
       else if (q1+(q2+q3)*0.9+k1+k2.le.20) then
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
            if (emp2cm.lt.50000.0) then
                !
                dip,0.0d0
                CCT-f12
                eeee(n1) = energy
                !
                dip,field(1)
                CCT-f12
                ex6d1=energy
                dip,field(2)
                CCT-f12
                ex6d2=energy
                dipau=(ex6d2-ex6d1)/(field(2)-field(1))
                dip6dx(n1)=dipau*TODEBYE
                !             
                dip,,field(1)
                CCT-f12
                ex6d1=energy
                dip,,field(2)
                CCT-f12
                ex6d2=energy
                dipau=(ex6d2-ex6d1)/(field(2)-field(1))
                dip6dy(n1)=dipau*TODEBYE
                !
                dip,,,field(1)
                CCT-f12
                ex6d1=energy
                dip,,,field(2)
                CCT-f12
                ex6d2=energy
                dipau=(ex6d2-ex6d1)/(field(2)-field(1))
                dip6dz(n1)=dipau*TODEBYE
                !
                put, xyz, $fname.xyz.dat, old
                !
                SHOW,energy;
                !
                xxx(n1) = '222'
                yyy(n1) = '111'
                !
                text ### HOOH dipole
                table,xxx,xOO,xOH1,xOH2,xHOO1,xHOO2,xHOOH,eeee,dip6dx,dip6dy,dip6dz
                DIGITS, 3,  5,   5,   5,   5,     5,    5,  12,    12,    12,    12 
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
table,yyy,xOO,xOH1,xOH2,xHOO1,xHOO2,xHOOH,eeee,dip6dx,dip6dy,dip6dz
DIGITS, 3,  5,   5,   5,   5,     5,    5,  12,    12,    12,    12 
          
save,$fname.dat,new


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

#setenv molproexe /shared/ucl/apps/Molpro/2012.1-bindist/bin/
setenv molproexe /shared/ucl/apps/Molpro/2009.1/bin/
#setenv PATH /shared/ucl/apps/Molpro/2009.1/bin/:"${PATH}"

echo "Running molpro -n 1 -d $TMPDIR -W $wdir < $pwd/$fname.com > $pwd/$fname.out"

$molproexe/molpro -n 1 -I $TMPDIR -d $TMPDIR -W $wdir < $pwd/$fname.com > $pwd/$fname.out


if (-e $fname.dat) then
   /bin/cp  $fname.dat  $pwd/dat-ccsd-f12-tz
endif

if (-e $fname.xyz.dat) then
   /bin/cp  $fname.xyz.dat  $pwd/dat-ccsd-f12-tz-xyz
endif

if (-e $fname.out) then
   /bin/cp  $fname.out  $pwd/out-ccsd-f12-tz
endif

#/bin/rm -rf $HOME/${fname}.codine
#/bin/rm -rf /usr/scratch/$USER/${fname}


