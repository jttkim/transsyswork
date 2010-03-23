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


function checkpython()
{
  if [ -z "$PYTHONPATH" ]
    then
    export PYTHONPATH=$HOME/lib64/python
  fi
}

function cleanFolder()
{

  echo "Cleaning folder"
  rm -f yn*
}


function maketp () # create transsys program
{

  num_topology=$1
  transsys_name=$2
  mxnode=$3
  nn=$4
  nnp=$5
  ne=$6
  nep=$7

  for (( i=${nn}; i<=${mxnode}; i++ ))
  do
    nnp=`expr $i`
    let "nep=$i*$indegree"
    sed -e '/n: / s/'${nn}'/'$nnp'/' -e '/num_edges: / s/'${ne}'/'$nep'/' $transsys_file > file.dat
    tname=`printf '%sn%de%d' $transsys_name $i $nep`
    do_run $HOME/makemodel-1.0/transsysrandomprogram -t $tname -m $num_topology -p file.dat
  done
}


function changeparameter () # Change transsys parameters
{

  num_topology=$1
  transsys_name=$2
  nn=$3
  rng_parameter=$4
  rndseed=$5

  for (( i=${nn}; i<=${mxnode}; i++ ))
  do
    for (( j=1; j<=${num_topology}; j++ ))
    do
      let "rndseed=$rndseed + 1"
      echo $rndseed
      tname=`printf '%sn%de%d%02d' $transsys_name $i $nep $j`
      do_run $HOME/makemodel-1.0/transsysreparam -r ${rndseed} -n ${tname} -p ${rng_parameter} -T transformerfile.dat ${tname}'.tra'
    done
  done

}


function makesimdata () # create simulated data
{
  
  num_topology=$1
  equilibration_length=$2
  transsys_name=$3
  mxnode=$4
  nn=$5
  nnp=$6
  rng_parameter=$7
  noise_rate=$8

  for (( i=$nn; i<=${mxnode}; i++ ))
  do
    let "nep=$i*$indegree"
    for (( j=1; j<=${num_topology}; j++))
    do
      for (( k=1; k<=${rng_parameter}; k++ ))
      do
        tname=`printf '%sn%de%d%02d%02d' $transsys_name $i $nep $j $k`
        #do_run $HOME/makemodel-1.0/transsyscreatedata -g $i expr.txt pheno.txt feature.txt
        val=$(echo $noise_rate | sed 's/0.//')
        if [[ $val -eq 0 ]]
        then
          do_run $HOME/makemodel-1.0/transsyswritesimset -e ${equilibration_length} -x expr.txt -p pheno.txt -f feature.txt ${tname}'.tra' ${tname}
        else 
          do_run $HOME/makemodel-1.0/transsyswritesimset -r 1 -s ${noise_rate} -e ${equilibration_length} -x expr.txt -p pheno.txt -f feature.txt ${tname}'.tra' ${tname}
        fi
      done    
    done
  done 
}


function rewire_net ()
{

  num_topology=$1
  transsys_name=$2
  mxnode=$3
  nn=$4
  nnp=$5
  rng_parameter=$6
  rnd_repetition=$7
  rw_operation_vector=$8

  for (( i=$nn; i<=${mxnode}; i++ ))
  do
    let "nep=$i*$indegree"
    for (( j=1; j<=${num_topology}; j++))
    do
      for (( k=1; k<=${rng_parameter}; k++ ))
      do
        tname=`printf '%sn%de%d%02d%02d' $transsys_name $i $nep $j $k`
        do_run $HOME/makemodel-1.0/transsysrewire -v $rw_operation_vector -t ${tname} -r ${rnd_repetition} $tname'.tra'
      done
    done
  done

}

function runfunc ()
{

  num_topology=$1
  equilibration_length=$2
  random_start=$3
  transsys_name=$4
  mxnode=$5
  nn=$6
  nnp=$7
  rng_parameter=$8
  rnd_repetition=$9

 
  for (( i=$nn; i<=${mxnode}; i++ ))
  do
    let "nep=$i*$indegree"
    for (( j=${num_topology}; j<=${num_topology}; j++))
    do
      for (( k=${rng_parameter}; k<=${rng_parameter}; k++ ))
      do
        echo $rndseedop
        tname=`printf '%sn%de%d%02d%02d' $transsys_name $i $nep $j $k`
        expr=`printf "$tname"_expr.txt""`
        pheno=`printf "$tname"_pheno.txt""`
        feature=`printf "$tname"_feature.txt""`
        if [ $logratio_mode -eq 1 ] 
          then
          do_run $HOME/makemodel-1.0/netoptrew -N ${j} -P ${k} -r ${rndseedop} -l TRUE -o $offset -w $rw_operation -R ${random_start} -e ${equilibration_length} -n ${rnd_repetition} -t ${tname} -x $expr -p $pheno -f $feature -u correlation -L $log_file -g $optimiser_file -T $transformer_file
        else 
          do_run $HOME/makemodel-1.0/netoptrew -N ${j} -P ${k} -r ${rndseedop} -w $num_operation -R ${random_start} -e ${equilibration_length} -n ${rnd_repetition} -t ${tname} -x $expr -p $pheno -f $feature -u correlation -L $log_file -g $optimiser_file -T $transformer_file
        fi
      done
    done
  done

}

num_topology=10
rnd_parameter=5
rnd_repetition=10
noise_rate=0
logratio_mode=1
indegree=2
equilibration_length=1000
random_start=5
verbose=0
filename='x'
transsys_name='y'
rw_operation_vector=24
transsys_file="None"
transformer_file="None"
optimiser_file="None"
rndseedop=1
mode=1
offset=0.01
rw_operation=0
log_file=None

while getopts s:t:p:w:v:n:r:l:o:i:L:P:T:O:m:W:d opt
do
  case "$opt" in
    t) num_topology="$OPTARG";;
    p) rnd_parameter="$OPTARG";;
    w) rw_operation="$OPTARG";;
    r) rnd_repetition="$OPTARG";;
    n) noise_rate="$OPTARG";;
    l) logratio_mode="$OPTARG";;
    i) indegree="$OPTARG";;
    m) mode="$OPTARG";;
    s) rndseedop="$OPTARG";;
    P) transsys_file="$OPTARG";;
    T) transformer_file="$OPTARG";;
    O) optimiser_file="$OPTARG";;
    L) log_file="$OPTARG";;
    o) offset="$OPTARG";;
    v) rw_operation_vector="$OPTARG";;
    d) isdef=1;;
    \?) help_ani;;
  esac
done


nn=$(grep -w "n:" $transsys_file | awk '{print$2}')
nnp=`expr ${nn}`
ne=$(grep -w "num_edges:" $transsys_file | awk '{print$2}')
nep=`expr ${ne}`
mxnode=${nn}
rndseed=$(grep -w "rndseed:" $transsys_file | awk '{print$2}')


if [ $mode -eq 1 ] 
  then
  checkpython
  maketp $num_topology $transsys_name $mxnode $nn $nnp $ne $nep
  changeparameter $num_topology $transsys_name $nn $rng_parameter $rndseed
  makesimdata $num_topology $equilibration_length $transsys_name $mxnode $nn $nnp $rng_parameter $noise_rate
  rewire_net $num_topology $transsys_name $mxnode $nn $nnp $rng_parameter $rnd_repetition $rw_operation_vector
 else
  checkpython
  runfunc $num_topology $equilibration_length $random_start $transsys_name $mxnode $nn $nnp $rnd_parameter $rnd_repetition
fi

