#!/bin/sh


function do_run ()
{
 echo $*
 if $* ; then
   true
 else
   exit 1
 fi
}


function maketp () # create transsys program
{

  num_topology=$1
  num_parameter=$2
  nn=$3
  ne=$4
  nw=$5

  for (( i=1; i<=${num_topology}; i++ ))
  do
    for (( j=1; j<=${num_parameter}; j++ ))
    do
      for k in 0 1 2 3 4 5 6 7 9 11 13 15 18 22 27 32 38 46 55 66
      # for k in 0 1
      do
        echo $k
        tname=`printf 'jobER100%02d%02d%02d' $i $j $k`
        sed -e '/-t / s/'${nn}'/'$i'/' -e '/-p / s/'${ne}'/'$j'/' -e '/-w / s/'${nw}'/'$k'/' jobmaster.sh > $tname'.sh'
        do_run qsub $tname'.sh'
     done
   done 
  done
}

function deletejob () {

  for (( i=0; i<40; i++))
    do
    num=`expr 233683 + $i`
    echo $num
    qdel $num
  done
}

num_topology=10
num_parameter=5
nn=$(grep -w "t" jobmaster.sh | awk '{print$3}')
ne=$(grep -w "p" jobmaster.sh | awk '{print$5}')
nw=$(grep -w "w" jobmaster.sh | awk '{print$7}')

maketp $num_topology $num_parameter $nn $ne $nw
#deletejob
