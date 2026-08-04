#!/bin/csh 
#
# Generation of the input file for MOLPRO 
#

cd $2
#cd /ds20th/san/scratch/yurchenk/ch4/molpro/
set pwd = `pwd`
echo $pwd

set irun = $1
set ffile = "./inp/ch4_grid1.dat"

#./read_geom.sh 1

set fname = ch4.${irun}.Q.f12c.kroll.1
set datfile = ch4.${irun}.Q.f12c.k.dat
#
cat<<endb> $fname.com

gthresh,twoint=1.d-10,prefac=1.d-18,zero=1.d-10,energy=1.d-10;
memory,512,m
!
dkroll=1
!

PROC CCT-f12
 {hf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,18;wf,orbital}
 {ccsd(t)-f12c; thresh,energy=1.d-10,coeff=1.d-9,thrint=1.d-10,zero=1.d-10;maxit,100}
 !{OPTG,grad=5,energy=8,saveact='xxx.opt';print,history;inactive,alpha}
 !{readvar,'xxx.opt'}
ENDPROC


PROC calc-mp2
 {hf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,18;wf,orbital}
 {mp2}
ENDPROC


BASIS=cc-pVQZ-F12
ANSATZ=3c,fix=1,canonical=1

!BASIS=ACVTZ


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
 eeee(i1)=energy
 !
 !gexpec,rel,darwin,massv
 !relen(i1) =erel
 !
 NZ(i1) = "XX"
 IZ(i1) = ${irun}
 !
 table,NZ,IZ,R1CH,R2CH,R3CH,R4CH,A12,A13,A24,A23,A14,eeee
 DIGITS,0, 0,   8,   8,   8,   8,  8,  8,  8,  8,  8,  12
 save,$fname.dat,new
 !
 title Energies of CH4 

---

endb

setenv OMP_NUM_THREADS 1

if (! $?COD_O_WORKDIR) then
   setenv COD_O_WORKDIR $SGE_O_WORKDIR
endif
limit


mkdir /usr/scratch/$USER/$fname
mkdir $HOME/${fname}.codine
# be sure /tmp will not be used
setenv TMP /usr/scratch/$USER/${fname}
setenv TMPDIR /usr/scratch/$USER/${fname}
#
cd $HOME/${fname}.codine
#/usr/local/bin/molpro -I /usr/scratch/$USER/$1 -d /usr/scratch/$USER/$1 < $COD_O_WORKDIR/$1.com >& $1.out
/usr/local/bin/molpro -I /usr/scratch/$USER/${fname} -d /usr/scratch/$USER/${fname} < $COD_O_WORKDIR/${fname}.com >& $COD_O_WORKDIR/${fname}.out
#/bin/cp $1.out $COD_O_WORKDIR

#/bin/cp $1.out $pwd

if (-e ch4*.dat) then
   /bin/cp  *.dat  $pwd
endif

cd ..
/bin/rm -rf $HOME/${fname}.codine
/bin/rm -rf /usr/scratch/$USER/${fname}

