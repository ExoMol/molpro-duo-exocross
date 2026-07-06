#!/bin/bash  -l
#


export name=`echo $1 | sed -e 's/\.states//'`
export name=`echo $name | sed -e 's/\.trans//'`


echo $name

export fname=i-${name}_g1_20000_1cm-1_T$2


echo $fname

cat<<endb> $fname.inp

mem 190.0 gb

Temperature $2
Range 0 20000

npoints 20000

(abundance)

absorption 
gaussian

hwhm 1

output $fname

NRAM 1000000

NPROCS 4

verbose 4

States $name.states

transitions $name.trans

endb

./j-xcross.x  <$fname.inp > $fname.out
