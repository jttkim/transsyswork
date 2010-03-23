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


function maketp () # create jobs
{

  num_topology=$1
  num_parameter=$2
  nn=$3
  ne=$4
  nw=$5
  nettype=$6
  
  for (( i=1; i<=${num_topology}; i++ ))
  do
    for (( j=1; j<=${num_parameter}; j++ ))
    do
      for k in 0 1 2 3 4 5 6 7 9 11 13 15 18 22 27 32 38 46 55 66
      do
        echo $k
        let "rndseed = $rndseed + 1"
        tname=`printf 'job%s%02d%02d%02d' $nettype $i $j $k`
        sed -e '/-t / s/'${nn}'/'$i'/' -e '/-p / s/'${ne}'/'$j'/' -e '/-w / s/'${nw}'/'$k'/' -e '/-s / s/'${rd}'/'$rndseed'/' $jobmaster > $tname'.sh'
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


function deletefiles () {

  for (( i=0; i<11; i++))
    do
    tname=`printf 'yn30e60%02d*.tra' $i`
    tnameexp=`printf 'yn30e60%02d*_expr.txt' $i`
    tnamephe=`printf 'yn30e60%02d*_pheno.txt' $i`
    tnamefea=`printf 'yn30e60%02d*_feature.txt' $i`
    rm -fr $tname
    rm -fr $tnameexp
    rm -fr $tnamephe
    rm -fr $tnamefea
    rm -fr *logo*
  done
}


function fixerror()
{

  num_topology=$1
  num_parameter=$2
  nn=$3
  ne=$4
  nw=$5
  nettype=$6

  for (( i=1; i<=${num_topology}; i++ ))
  do
    for (( j=1; j<=${num_parameter}; j++ ))
    do
      for k in 0 1 2 3 4 5 6 7 9 11 13 15 18 22 27 32 38 46 55 66
      do
        tname=`printf 'yn30e60%02d%02d%02dlogo.txt' $i $j $k`
        tname1=`printf 'yn30e60%02d%02d%02dlogo_op.txt' $i $j $k`
        key=`printf '%02d%02d' $i $j`
        t=`expr '0'$i'0'$j$k`
        key1=`printf '%02d\t%02d\t' $i $j`
        echo $t
        sed -n -e '/${t}/p' $tname
     done
   done 
  done
}


nettype=ER

while getopts m:t:d opt
do
  case "$opt" in
    m) mode="$OPTARG";;
    t) nettype="$OPTARG";;
    d) isdef=1;;
    
    \?) help_ani;;
  esac
done


num_topology=1
num_parameter=1
jobmaster=$HOME/makemodel-1.0/jobmaster.sh

nn=$(grep -w "t" $jobmaster | awk '{print$3}')
ne=$(grep -w "p" $jobmaster | awk '{print$5}')
nw=$(grep -w "w" $jobmaster | awk '{print$7}')
rd=$(grep -w "s" $jobmaster | awk '{print$9}')
rndseed=12

if [ $mode -eq 1 ]
  then
  maketp $num_topology $num_parameter $nn $ne $nw $nettype $rndseed
elif [ $mode -eq 2 ]
  then
  deletefiles
elif [ $mode -eq 3 ]
  then
  deletejob
else
  fixerror $num_topology $num_parameter $nn $ne $nw
fi

#$HOME/makemodel-1.0/modelone.sh -m 1
