#!/bin/bash
#
# Generation of the input file for MOLPRO
#

echo "Start Molpro"

export pwd=`pwd`


export name=`echo $1 | sed -e 's/\.com//'`


export fname=$name

export outdir=${fname}

if [ -e $fname.out ]; then
   /bin/cp $fname.out $fname.tmp
fi


if [ ! -e $pwd/${outdir} ]; then
   mkdir $pwd/${outdir}
fi


export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1



# be sure /tmp will not be used

export TMP=$pwd/${fname}
export TMPDIR=$pwd/${fname}
cd $TMPDIR
#


module load molpro/2024

molpro -I $TMPDIR -d $TMPDIR -W $TMPDIR < $pwd/${fname}.com > $pwd/${fname}.out

cd $pwd

grep " xxxx " $pwd/${fname}.out | cut -c 6-35  > pec.txt
grep " yyyy " $pwd/${fname}.out | cut -c 6-35  > dmc.txt

#if (-e $fname.dat) then
#   /bin/mv  $fname.dat     $pwd/${outdir}-dat
#endif

#if (-e $fname.molden) then
#   /bin/mv  $fname.molden  $pwd/${outdir}-out
#endif

#if (-e $pwd/$fname.out) then
#    gzip $pwd/$fname.out
#   /bin/mv  $pwd/$fname.out.gz  $pwd/${outdir}-out
#endif


#if (-e $fname.xyz.dat) then
#   /bin/mv  $fname.xyz.dat  $pwd/${outdir}-dat
#endif


#if (-e $pwd/$fname.com) then
#   /bin/rm  $pwd/$fname.com
#endif


cd $pwd


/bin/rm -rf $TMPDIR


