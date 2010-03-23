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


function getData()
{

  do_run ./testonenet.sh -t $num_topology -p $rnd_parameter -r $rnd_repetition -v $rw_operation_vector -l $logratio_mode -i $indegree -n $noise_rate -P $transsys_file -T $transformer_file -O $optimiser_file -m $mode

}


function runOptimiser()
{
  echo $transsys_file
  do_run ./testonenet.sh -t $num_topology -p $rnd_parameter -w $rw_operation -s $rndseedop -r $rnd_repetition -l $logratio_mode -i $indegree -L $log_file -P $transsys_file -T $transformer_file -O $optimiser_file -o $offset -m $mode 

}


num_topology=1
rnd_parameter=1
rw_operation=0
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
offset=0.01
mode=1
transsys_file="parapl.dat"
transformer_file="transformerfile.dat"
optimiser_file="optspec.dat"
log_file="logo"


while getopts t:p:w:r:v:n:l:i:s:e:P:T:O:L:m:o:d opt
do
  case "$opt" in
    t) num_topology="$OPTARG";; # Topology
    p) rnd_parameter="$OPTARG";; # Subnetworks with parameters changed at random
    r) rnd_repetition="$OPTARG";; # Random repetitions
    w) rw_operation="$OPTARG";; # Rewiring operation - Test on taht rewiring operation
    n) noise_rate="$OPTARG";; # Noise rate
    l) logratio_mode="$OPTARG";; #Logratio mode 1. Logration 0. Non logratio
    i) indegree="$OPTARG";; # Indegree
    s) rndseedop="$OPTARG";; # Random seed
    e) equilibration_length="$OPTARG";; # Equilibration length
    P) transsys_file="$OPTARG";; # Transsys file
    T) transformer_file="$OPTARG";; # Transformer file
    O) optimiser_file="$OPTARG";; # Optimiser file
    L) log_file="$OPTARG";; # Log file
    m) mode="$OPTARG";; #Testing mode
    o) offset="$OPTARG";; # Offset for data shiftting
    v) rw_operation_vector="$OPTARG";; # rewiring operation vector
    d) isdef=1;;
    \?) help_ani;;
  esac
done

echo $transformer_file
if [ $mode -eq 1 ] 
  then
  getData
 else
  runOptimiser
fi


getData -t 1 -p 1 -n 0 -m 1
runOptimiser -t 1 -p 1 -w 0 -s 13 -r 10 -l 1 -i 2 -m 2
