#!/bin/bash 

export pwd=`pwd`


export fname=`echo $1 | sed -e 's/\.template//'`

export dname=`echo $2 | sed -e 's/\.template//'`

export template=$fname

cat $template.template >  $fname.inp
cat pec.txt >>  $fname.inp
cat<<end1>> $fname.inp
end
end1

cat $dname.template >>  $fname.inp
cat dmc.txt >>  $fname.inp
cat<<end2>> $fname.inp
end
end2

export linelist=${fname}

cat<<end3>> $fname.inp
intensity 
absorption
THRESH_INTES 1e-60
THRESH_LINE 1e-40
TEMPERATURE 2000
thresh_dipole 1e-9
linelist ${linelist}
J 0 - 40
freq-window 0,20000
energy low 0, 20000 upper 0, 40000
end
end3

cd $pwd

./j-duo.x <$fname.inp > $fname.out 

./h_g1_T.sh ${linelist}.states  1000
