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

  $HOME/makemodel-1.0/testonenet.sh -t $num_topology -p $rnd_parameter -r $rnd_repetition -v $rw_operation_vector -l $logratio_mode -i $indegree -n $noise_rate -P $transsys_file -T $transformer_file -O $optimiser_file -m $mode

}


function runOptimiser()
{

  $HOME/makemodel-1.0/testonenet.sh -t $num_topology -p $rnd_parameter -w $rw_operation -s $rndseedop -r $rnd_repetition -l $logratio_mode -i $indegree -L $log_file -P $transsys_file -T $transformer_file -O $optimiser_file -o $offset -m $mode 
 # $HOME/makemodel-1.0/testonenet.sh -l 1 -o $offset -w $rewiring_op -r $random_start -e $equilibration_length -n $rnd_repetition -L $log_file -t yn30e600101 -x yn30e600101_expr.txt -p yn30e600101_pheno.txt -f yn30e600101_feature.txt -u correlation -g optspec.dat -T transformerfile.dat

}


num_topology=1
rnd_parameter=5
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
transsys_file="$HOME/makemodel-1.0/parapl.dat"
transformer_file="$HOME/makemodel-1.0/transformerfile.dat"
optimiser_file="$HOME/makemodel-1.0/optspec.dat"
log_file="logo"


while getopts t:p:w:r:v:n:l:i:s:e:P:T:O:L:m:o:d opt
do
  case "$opt" in
    t) num_topology="$OPTARG";; # Topology
    p) rnd_parameter="$OPTARG";; # Subnetworks with parameters changed at random
    r) rnd_repetition="$OPTARG";; 
    w) rw_operation="$OPTARG";; # Number of rewiring operations
    n) noise_rate="$OPTARG";;
    l) logratio_mode="$OPTARG";;
    i) indegree="$OPTARG";;
    s) rndseedop="$OPTARG";;
    e) equilibration_length="$OPTARG";;
    P) transsys_file="$OPTARG";;
    T) transformer_file="$OPTARG";;
    O) optimiser_file="$OPTARG";;
    L) log_file="$OPTARG";;
    m) mode="$OPTARG";;
    o) offset="$OPTARG";;
    v) rw_operation_vector="$OPTARG";;
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


# $HOME/makemodel-1.0/runprog.sh -t 1 -p 1 -n 0 -m 1
#$HOME/makemodel-1.0/runprog.sh -t 1 -p 1 -w 0 -s 13 -r 10 -l 1 -i 2 -m 2
