#!/bin/csh

cd $3
set pwd = `pwd`
echo $pwd
set script = $4
echo $script

set step = 0
set omp = 4

set i = $1
set k = $2
set k = `expr $k + 1`
set j = `expr $i + $step`

echo $i,$k
while ($i<$k);
 echo $i $j "("{$i}")";
 ./$script   $i $pwd $5 $6 $7
 set i = `expr $j + 1`
 set j = `expr $i + $step`
end

