#!/bin/sh


function do_run()
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

function dummy()
{
  dummyL=1
  for i in $dummy ; do
    dummyL=`expr $dummyL + 1 `
  done
}

function dummy2()
{
  dummyL=1
  if [[ $dummy -eq 1 ]]
  then
    dummyL=`expr $dummyL + 1 `
  fi
}

function dummy3()
{
  dummyL=1
  while test $dummy -le $dummyL ; do
    dummyL=`expr $dummyL + 1 `
  done
}

function createEmpiricalData()
{
  tp_name=`printf '%s' ${true_model} `
  do_run ./transsyswritesimsetOF -o  ${simgenex_model} -s ${rndseed} -N ${signal_to_noise} ${true_model}'_m00.tra' ${tp_name} ${tp_name}'_trace.txt'
}

function optimiseModelSynthetic()
{
  tp_name=`printf '%s' ${true_model} `
  for (( model=0; model<=${num_model}; model++ )) ; do
    model_name=`printf '%s_m%02d' ${tp_name} ${model}`
    do_run ./netopt -x ${tp_name}'_expr.txt' -o ${simgenex_model} -R ${num_optimisations} -g ${gradientfile} -T ${transformerfile} -L ${logfile} -s ${rndseed} -c ${model_name}
  done
}


function optimiseModelEmpiric()
{

  for model in ${tp_candidate[@]} 
  do
    do_run ./netopt -x  ${empirical_data} -o ${simgenex_model} -R ${num_optimisations} -g ${gradientfile} -T ${transformerfile} -L ${logfile} -s ${rndseed} -c ${model} 'opt_'${model}'.txt' 'tp_'${model}'.txt'
  done
}


function mergeFile()
{
  tp_name=`printf '%s' ${true_model} `
  for (( model=0; model<=${num_model}; model++ )) ; do
    candidate_topology=`printf '%s_m%02d' ${tp_name} ${model}`
    candidate_topology_logo=`printf '%s_logo.txt' ${candidate_topology}`
    echo $candidate_topology
    cp $candidate_topology_logo $candidate_topology_logo'.bk'
    sed 's/rst/'$model'\t/g' $candidate_topology_logo'.bk' > clean.txt
    mv clean.txt $candidate_topology_logo
    rm -f $candidate_topology_logo'.bk'
    name=`printf '%s %s ' $name $candidate_topology_logo `
  done
  fitnessname=`printf 'fitnesstable' `
  rm -fr $fitnessname'.txt'
  rm -f 'tt.txt'
  label=`printf 'model\trestart\tfitness'`
  cat $name > 'tt.txt'
  echo $label >> $fitnessname'.txt'
  sed '/restart/d' 'tt.txt' >> $fitnessname'.txt'
  cat $fitnessname'.txt'
  rm -f 'tt.txt'
}


num_model=0
simgenex_model=jasmonate_model.txt
signal_to_noise=0
rndseed=2
true_model=jasmonate
num_optimisations=3
transformerfile=transformerfile.dat
logfile=logo
gradientfile=optspec.dat
empirical_data=emp_data_jas.txt
tp_candidate=none

while getopts m:e:s:o:d opt
do
  case "$opt" in
    m) tp_candidate="$OPTARG";;
    e) empirical_data="$OPTARG";;
    s) simgenex_model="$OPTARG";;
    o) num_optimisations="$OPTARG";;
    d) isdef=1;;
    \?) help_ani;;
  esac
done


checkpython
#createEmpiricalData
#optimiseModelSynthetic
optimiseModelEmpiric
#mergeFile
