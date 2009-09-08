#!/bin/sh

function do_run ()
{
  echo $*
  # return
  if $* ; then
    true
  else
    exit 1
  fi
}


function jtprndname ()
{
  printf '%srnd%03d' ${1} ${2}
}


function jtpmutname ()
{
  printf '%smut%02d_%03d' ${1} ${2} ${3}
}


function makeprgset_regroupings ()
{
  basename=$1
  nrnd=$2
  genes=$3
  control_opts=''
  if test -n "$genes" ; then
    control_opts="$control_opts -n $genes"
  fi
  do_run ./makejasmonate $control_opts prg ${basename}prg.tra
  do_run ./makejasmonate -d $control_opts prg ${basename}prg_dummy.tra
  do_run ./makejasmonate $control_opts blank ${basename}blank.tra
  rndseed=1
  while test $rndseed -le $nrnd ; do
    tfilename=`jtprndname $basename $rndseed`.tra
    do_run ./makejasmonate $control_opts -s $rndseed rnd $tfilename
    rndseed=`expr $rndseed + 1`
  done
}


function makeprgset_perturbations ()
{
  basename=$1
  nrnd=$2
  num_perturbations_list="$3"
  genes=$4
  control_opts=''
  if test -n "$genes" ; then
    control_opts="$control_opts -n $genes"
  fi
  rndseed=1
  while test $rndseed -le $nrnd ; do
    echo $rndseed
    for num_perturbations in $num_perturbations_list ; do
      tfilename=`jtpmutname $basename $num_perturbations $rndseed`.tra
      do_run ./makejasmonate $control_opts -m $num_perturbations -s $rndseed mut $tfilename
    done
    rndseed=`expr $rndseed + 1`
  done
}


function write_cshscript_header ()
{
  scriptname=$1
  echo '#!/bin/csh' > $scriptname
  echo '#$ -q long.q' >> $scriptname
  echo '#$ -j y' >> $scriptname
  echo '#$ -cwd' >> $scriptname
  echo >> $scriptname
  echo 'setenv PYTHONPATH ${HOME}/lib/python' >> $scriptname
  echo "cd $PWD" >> $scriptname
}


function add_command ()
{
  scriptname=$1
  cmd="$2"
  # echo adding command to $scriptname
  echo "echo $cmd" >> $scriptname
  echo "if ( { $cmd } ) then" >> $scriptname
  echo "  echo completed \`date\`" >> $scriptname
  echo "else" >> $scriptname
  echo "  echo ERROR" >> $scriptname
  echo "  exit 1" >> $scriptname
  echo "endif" >> $scriptname
}


function add_qsub ()
{
  last_scriptname=$1
  scriptname=$2
  echo "if ( { qsub $scriptname } ) then" >> $last_scriptname
  echo "  echo submitted $qsubscriptname \`date\`" >> $last_scriptname
  echo "else" >> $last_scriptname
  echo "  echo ERROR submitting $qsubscriptname" >> $last_scriptname
  echo "  exit 1" >> $last_scriptname
  echo "endif" >> $last_scriptname
}


function make_objective_chain ()
{
  jtpname=$1
  control_opts="$2"
  outfile_prefix=$3
  objectivefunction_list="squaresum correlation"
  last_scriptname=''
  for objectivefunction in $objectivefunction_list ; do
    outfile_basename=${outfile_prefix}_${objectivefunction}
    scriptname=${outfile_basename}.csh
    write_cshscript_header $scriptname
    echo $scriptname
    if test -n "$last_scriptname" ; then
      add_qsub $last_scriptname $scriptname
    fi
    last_scriptname=$scriptname
    output_opts="-L ${outfile_basename}_log.dat -p ${outfile_basename}_profiles.dat"
    add_command $scriptname "./jasmonopt $control_opts -f $objectivefunction $output_opts ${jtpname}.tra ${outfile_basename}_opt.tra"
    echo $scriptname
    add_command $scriptname "gzip ${outfile_basename}_log.dat"
    echo $scriptname

    # logratio portion of loop body
    outfile_basename=${outfile_prefix}_${objectivefunction}_logratio
    scriptname=${outfile_basename}.csh
    write_cshscript_header $scriptname
    if test -n "$last_scriptname" ; then
      add_qsub $last_scriptname $scriptname
    fi
    last_scriptname=$scriptname
    output_opts="-L ${outfile_basename}_log.dat -p ${outfile_basename}_profiles.dat"
    add_command $scriptname "./jasmonopt $control_opts -f $objectivefunction -l -o 10.0 $output_opts ${jtpname}.tra ${outfile_basename}_opt.tra"
    add_command $scriptname "gzip ${outfile_basename}_log.dat"
  done
  if test -n "$last_scriptname" ; then
    echo "echo chain complete" >> $last_scriptname
  fi
}


function make_chains_regroupings ()
{
  jbase=$1
  optspec=$2
  transformspec=$3
  equilibration_length=$4
  minutes_per_timestep=$5
  num_restarts=$6
  factorlist=$7
  genelist=$8
  nrnd=$9
  control_opts="-x jdata_aggr_final.txt -m $minutes_per_timestep -e $equilibration_length -P ${optspec}.txt -T ${transformspec}.txt -r ${num_restarts}"
  if test -n "$factorlist" ; then
    control_opts="$control_opts -F $factorlist"
  fi
  if test -n "$genelist" ; then
    control_opts="$control_opts -G $genelist"
  fi
  for jtpext in prg blank ; do
    jtpname=${jbase}${jtpext}
    outfile_prefix="test${jtpname}_${optspec}_${transformspec}"
    make_objective_chain $jtpname "$control_opts" $outfile_prefix
  done
  rndseed=1
  while test $rndseed -le $nrnd ; do
    jtpname=`jtprndname ${jbase} ${rndseed}`
    outfile_prefix="test${jtpname}_${optspec}_${transformspec}"
    make_objective_chain $jtpname "$control_opts" $outfile_prefix
    rndseed=`expr $rndseed + 1`
  done
}


function make_chains_perturbations ()
{
  jbase=$1
  optspec=$2
  transformspec=$3
  equilibration_length=$4
  minutes_per_timestep=$5
  num_restarts=$6
  factorlist=$7
  genelist=$8
  nrnd=$9
  nmut_list="${10}"
  control_opts="-x jdata_aggr_final.txt -m $minutes_per_timestep -e $equilibration_length -P ${optspec}.txt -T ${transformspec}.txt -r ${num_restarts}"
  if test -n "$factorlist" ; then
    control_opts="$control_opts -F $factorlist"
  fi
  if test -n "$genelist" ; then
    control_opts="$control_opts -G $genelist"
  fi
  for nmut in $nmut_list ; do
    rndseed=1
    while test $rndseed -le $nrnd ; do
      jtpname=`jtpmutname ${jbase} ${nmut} ${rndseed}`
      outfile_prefix="test${jtpname}_${optspec}_${transformspec}"
      rndseed=`expr $rndseed + 1`
      make_objective_chain $jtpname "$control_opts" $outfile_prefix
    done
  done
}


num_restarts_regroupings=50
num_restarts_perturbations=5
num_regrouped_programs=20
num_perturbed_programs=10
optimiser=grad

num_perturbations_list='0 2 4 6 8 10 12 14 17 20 24 28 33'
fullfactorlist=At1g18710,At1g23080,At1g74020,At2g34810,At4g11290,At4g22880,At4g23600,At5g13930,At5g19890,At5g57090,At3g23050,At1g17740,At2g21130,At2g28900,At2g42530,At4g24360,At5g63790,At1g20440,At2g26690,At2g40900,At3g61890,At4g34000,At4g37760,At5g06760,At2g23320,At2g29500,At2g38470,At4g11280,At4g17880,At5g62470,At4g35770,At1g49570,At4g11320
fullgenelist=`echo $fullfactorlist | sed -e 's/,/_gene,/g' -e 's/$/_gene/'`


makeprgset_regroupings jfull $num_regrouped_programs $fullfactorlist
makeprgset_perturbations jfull $num_perturbed_programs "$num_perturbations_list" $fullfactorlist

make_chains_regroupings jfull $optimiser tfhypplus 100 5.0 $num_restarts_regroupings $fullfactorlist $fullgenelist $num_regrouped_programs
make_chains_perturbations jfull $optimiser tfhypplus 100 5.0 $num_restarts_perturbations $fullfactorlist $fullgenelist $num_perturbed_programs "$num_perturbations_list"

