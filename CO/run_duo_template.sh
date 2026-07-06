#!/bin/bash 

export pwd=`pwd`


export fname=`echo $1 | sed -e 's/\.template//'`


export template=$fname

cat $template.template >  $fname.inp
cat pec.txt >>  $fname.inp
cat<<endb>> $fname.inp
end
endb

cd $pwd

./j-duo.x <$fname.inp > $fname.out 
