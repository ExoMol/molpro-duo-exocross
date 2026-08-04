#!/bin/bash -l 
#
# Generation of the input file for MOLPRO 
#

cd $2
export pwd=`pwd`
echo $pwd

export irun=$1
export ijob=$3

export outdir=tz_14.${ijob}

echo $outdir

#5z-7220_noSO_THRNRM_G2_${ijob}


if [ ! -e $pwd/${outdir}-dat ]; then
   mkdir $pwd/${outdir}-dat
fi

if [ ! -e $pwd/${outdir}-out ]; then
   mkdir $pwd/${outdir}-out
fi


export fname=hcn.${irun}.${outdir}

export dname=${outdir}.${irun}


echo $fname
#
cat<<endb> $fname.com

emp20 =  -65.654108026081

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
 !{ccsd(t)-f12b,thrden=1.0d-12,ri_basis=ri_basis_set,df_basis=AWCV5Z/MP2FIT,df_basis_exch=V5Z/JKFIT,gem_beta=1.2;
 !thresh,energy=1.d-10,coeff=1.d-10,thrint=1.d-10,zero=1.d-11;maxit,100}
 {ccsd(t)-f12c,thrden=1.0d-10,ri_basis=optri,df_basis=awcv5z/mp2fit,df_basis_exch=v5z/jkfit,gem_beta=1.0;
 thresh,energy=1.d-10,coeff=1.d-10,thrint=1.d-10,zero=1.d-11;maxit,100}
 !{STATUS,nocheck}
 !{OPTG,grad=5,energy=8,saveact='xxx.opt';print,history;inactive,HCN}
 !{OPTG,grad=5,energy=8,saveact='xxx.opt';print,history;inactive,OO}
 !{readvar,'xxx.opt'}
ENDPROC


ANSATZ=3C,fix=1,hybrid=1,canonical=1
basis=vtz-f12
!basis=vtz
ansatz=3c,fix=1,canonical=1



symmetry,nosym;
orient,noorient;

geometry={angstrom;
  C;
  N, 1, CN;
  H, 1, CH, 2, HCN;
}


grida   = [180,178,176,175,172,170,160,150,140,130,120,110,100,90,80,70,60,50,40,30,20,15,10,8,6,4,2,0]
gridr1  = [1.15484278,1.15486448,1.15492955,1.15497832,1.15518931,1.1553836,1.15698534,1.15958837,1.1630929,1.16734426,1.17210501,1.17701178,1.18153515,1.18505036,1.18725981,1.18893572,1.19203105,1.19614817,1.19552312,1.18955997,1.18093325,1.17682607,1.17350459,1.17248493,1.17167147,1.17107925,1.17071918,1.17059833] Ang
gridr2  = [1.06657279,1.06658792,1.06663335,1.06666746,1.06681569,1.06695313,1.06812166,1.07016601,1.07325651,1.07767977,1.08388784,1.09260394,1.10503604,1.12327081,1.15097892,1.19491678,1.27234685,1.42688154,1.65646605,1.87156776,2.03440777,2.09233499,2.13379951,2.14573688,2.15501616,2.16164041,2.16561286,2.16693658] Ang


!grida   = [180,178,176,175,172,170,165,160,150,140,130,120,110,100,90,80,70,65,60,55,50,40,35,30,25,20,15,10,8,6,4,2,0]
!gridr1  = [1.06657279,1.06658792,1.06663335,1.06666746,1.0667092,1.06681569,1.06695313,1.06743496,1.06812166,1.06902602,1.07016601,1.07156573,1.07325651,1.075278,1.07767977,1.08052373,1.08388784,1.08787166,1.09260394,1.09825216,1.10503604,1.11324718,1.12327081,1.13562596,1.15097892,1.17026152,1.19491678,1.22744944,1.27234685,1.33669199,1.42688154,1.53847032,1.65646605,1.7693084,1.87156776,1.96049118,2.03440777,2.09233499,2.13379951,2.14573688,2.15501616,2.16164041,2.16561286,2.16693658] Ang
!gridr2  = [1.15484278,1.15486448,1.15492955,1.15497832,1.15503788,1.15518931,1.1553836,1.15605473,1.15698534,1.15816687,1.15958837,1.16123608,1.1630929,1.16513769,1.16734426,1.16968014,1.17210501,1.17456913,1.17701178,1.17936079,1.18153515,1.18345375,1.18505036,1.18630499,1.18725981,1.18805679,1.18893572,1.19018926,1.19203105,1.1942973,1.19614817,1.19661894,1.19552312,1.1930789,1.18955997,1.18534974,1.18093325,1.17682607,1.17350459,1.17248493,1.17167147,1.17107925,1.17071918,1.17059833] Ang

gridd =  [0.0,0.005,-0.005,0.01,-0.01,0.015,-0.015,0.02,-0.02,0.025,-0.025,0.03,-0.03,0.04,-0.04,-0.05,0.05,0.06,-0.06,0.08,-0.08,0.10,-0.10,0.15,-0.15,0.20,-0.20,0.3,-0.3,0.4,-0.4,0.5,-0.5,0.6,-0.6,0.7,-0.7.0.8,-0.8,0.9,-0.9,1.0,-1.0] Ang


!gridr1  = [1.066,1.42,1.38,1.44,1.36,1.50,1.32,1.55,1.25,1.60,1.20,1.70,1.10,1.80,1.10,2.00,1.00,3.00,4.00] 
!gridr2  = [1.153,2.02,1.98,2.05,1.96,2.10,1.92,2.15,1.90,2.20,1.85,2.50,1.80,2.80,1.70,3.00,4.00]
!grida   = [0.0,2.0,4.0,6.0,8.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0,95.0,100.0,105.0,110.,115.0,120.0,125.0,130.0,135.0,140.0,145.0,150.,155.0,160.0,165.0,170.0,172.0,174.0,175.0,176.0,178.0,180.0] 

i1 = 1
field=[-0.005, 0.005]

N1max= #gridd
N2max= #gridd
N3max= #grida
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


do qN=0,N1max-1
 do k3=1,N3max
   do q2=1,qN+1
       !
       q1 = qN-(q2-1)+1
       !
       show,k3,q1,q2
       !
       i=i+1
       !
       !if (q1.eq.1.or.q2.le.1) then
       !   irunjob = 0;
       !else if (k1.eq.k2.and.q1+q2+q3.le.8) then
       ! irunjob = 2;
       !else if (q2.eq.q3.and.k1.eq.k2.and.q1+q2+q3+k1+k2.le.22) then
       ! irunjob = 3;
       !else if (q1+q2+q3.le.10.and.k1+k2.le.8) then
       ! irunjob = 4;
       if (q1+q2.le.42) then
        irunjob = 6;
       !if (mod(i,200).eq.0) then	
       ! irunjob = 10;
       else
        irunjob = 0;
       end if;
       !
       SHOW,irunjob;
       !
       if (irunjob.ne.0) then
         !
         i1 = i1+1
         !
         if (imin.le.i1.and.i1.le.imax) then
            !
            n1 = n1 + 1 
            !
            SHOW,n1;
            !
            CN  = gridr1(k3)+gridd(q1)
            CH  = gridr2(k3)+gridd(q2)
            HCN  = grida(k3)

            CCT-f12
            eeee(n1) = energy
            !
            SHOW,energy;
            !
            xCN(n1)  = CN
            xCH(n1)  = CH
            xHCN(n1) = HCN
            !
            xxx(n1) = 'opqr'
            yyy(n1) = 'abcd'
            !
            text ### HCN TZ
            table,xxx,xCN,xCH,xHCN,eeee
            DIGITS, 3,  8,  8,   5,  10
            !save,$fname.dat,new
            !
            goto,END:
            !
         end if
       end if
  end do   
 end do    
end do 


END: TEXT ### HCN

text ### HCN 
table,yyy,xCN,xCH,xHCN,eeee
DIGITS, 3,  8,  8,   5,  10

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

#module load molpro/2015.1.5/intel-2015-update2
module load molpro/2020.1/openmp

#df

echo "Running molpro -n 1 -d $TMPDIR -W $wdir < $pwd/$fname.com > $pwd/$fname.out"

time molpro -n 1 -I $TMPDIR -d $TMPDIR -W $wdir < $pwd/$fname.com > $pwd/$fname.out

#$molproexe/molpro -n 1 -I $TMPDIR -d $TMPDIR -W $wdir < $pwd/${fname}.com > $wdir/${fname}.out



if [ -e $pwd/$fname.dat ]; then
   /bin/cp  $fname.dat  $pwd/${outdir}-dat
fi

if [ -e $pwd/$fname.xyz.dat ]; then
   /bin/cp  $fname.xyz.dat  $pwd/${outdir}-dat
fi

if [ -e $pwd/$fname.out ]; then
    gzip $pwd/$fname.out
   /bin/mv  $pwd/$fname.out.gz  $pwd/${outdir}-out
fi

if [ -e $pwd/$fname.com ]; then
     /bin/rm $pwd/$fname.com
fi


if [ -e $pwd/opt.act ]; then
   /bin/mv  $pwd/opt.act  $pwd/${outdir}-dat
fi

if [ -e $pwd/xxx.opt ]; then
   /bin/mv  $pwd/xxx.opt  $pwd/${outdir}-dat/$fname.opt
fi
