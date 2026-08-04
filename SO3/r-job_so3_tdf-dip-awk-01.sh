#!/bin/csh 
#
# Generation of the input file for MOLPRO 
#

cd $2
set pwd = `pwd`
echo $pwd

set irun = $1


set fname = SO33.${irun}.ccsdt


set outdir = $fname


if ( ! -e $pwd/${outdir}-dat) then
   mkdir $pwd/${outdir}-dat
endif

if ( ! -e $pwd/${outdir}-out) then
   mkdir $pwd/${outdir}-out
endif



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
END



symmetry,nosym;
orient,noorient;

geometry={angstrom;
 S;
 O, 1, R1SO;
 O, 1, R2SO, 2, AA3;
 O, 1, R3SO, 2, AA2, 3, AA1, 1;
}

hartree = 219474.63
emp20 = -623.083175267015

!stretching
jt0 = 1
gRR=[1.43, 1.44, 1.42, 1.46, 1.40, 1.35, 1.60, 1.30, 1.70, 1.10, 2.00]
!bending
gAA=[120., 119., 121.0, 117.0, 123., 115., 110.0, 100.0, 90.0, 70.0 ]


field=[-0.005, 0.005]

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
                dip,0.0d0
                CCT
                eee(i1) = energy
                !
                dip,field(1)
                CCT
                ex6d1=energy
                dip,field(2)
                CCT
                ex6d2=energy
                dipau=(ex6d2-ex6d1)/(field(2)-field(1))
                dip6dx(i1)=dipau*TODEBYE
                !             
                dip,,field(1)
                CCT
                ex6d1=energy
                dip,,field(2)
                CCT
                ex6d2=energy
                dipau=(ex6d2-ex6d1)/(field(2)-field(1))
                dip6dy(i1)=dipau*TODEBYE
                !
                dip,,,field(1)
                CCT
                ex6d1=energy
                dip,,,field(2)
                CCT
                ex6d2=energy
                dipau=(ex6d2-ex6d1)/(field(2)-field(1))
                dip6dz(i1)=dipau*TODEBYE
                !
                put, xyz, dxyz.tz.${irun}.dat, old
                !
                SHOW,energy;
                !
                xxx(n1) = '222'
                yyy(n1) = '111'
                !
                text ### SO3 energies 
                table,xxx,R1SOt,R2SOt,R3SOt,AAt1,AAt2,AAt3,eeee,dip6dx,dip6dy,dip6dz
                DIGITS, 3,    5,    5,    5,   5,   5,   5,  12,    12,    12,    12 
                !
             endif
             !
          end if
      end  if
   end do 
   end do 
   end do 
 end do   
 end do   
end do    
          
          
text ###  SO3 nergies 
table,yyy,R1SOt,R2SOt,R3SOt,AAt1,AAt2,AAt3,eeee,dip6dx,dip6dy,dip6dz
DIGITS, 3,    5,    5,    5,   5,   5,   5,  12,    12,    12,    12 
          
save,so3.dip3.${irun}.dat,new
          
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

module load gcc-libs/4.9.2 module unload compilers/intel/2015/update2
module load compilers/gnu/4.9.2 module unload
mpi/intel/2015/update3/intel module load mpi/openmpi/1.8.4/gnu-4.9.2
module load molpro/2015.1.3


echo "Running molpro -n 1 -d $TMPDIR -W $wdir < $pwd/$fname.com > $pwd/$fname.out"

time molpro -n 1 -I $TMPDIR -d $TMPDIR -W $wdir < $pwd/${fname}.com > $pwd/${fname}.out




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

