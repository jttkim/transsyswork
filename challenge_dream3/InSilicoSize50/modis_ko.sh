#!/bin/sh


modeldiscrimination=$HOME/transsyswork/modeldiscrimination/trunk
transsysrewire=$modeldiscrimination/transsysrewire
netoptrewGold=$modeldiscrimination/netoptrewGold
echo $transsysrewire

function do_run ()
{
  echo $*
  if $* ; then
  #if echo $* ; then
    true
  else
    echo '*** ERROR ***' >& 2
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


# naming scheme:
# target topology: target_<gentype>##.tra
# target program:  target_<gentype>##_p##.tra


function generate_candidate_programs ()
{
  for gentype in ${gentype_list} ; do
    target_topology_number=1
    while test ${target_topology_number} -le ${num_target_topologies} ; do
      target_topology_name=`printf 'target_%s%02d' ${gentype} ${target_topology_number}`
      target_topology_file=${target_topology_name}.tra
      candidate_topology_basename=`printf 'candidate_%s%02d' ${gentype} ${target_topology_number}`
      tp_name=`printf '%s' ${candidate_topology_basename}`
      for num_rewirings in ${num_rewirings_list} ; do
        do_run $transsysrewire -w ${num_rewirings} -n ${tp_name} -r ${num_rewired_topologies} -s ${rndseed} ${target_topology_file}
        rndseed=`expr ${rndseed} + 1`
      done
      target_topology_number=`expr ${target_topology_number} + 1`
    done
  done
}


function optimise_numrewired ()
{
  rewired_topology_number=1

  while test ${rewired_topology_number} -le ${num_rewired_topologies} ; do
    candidate_topology=`printf '%s_w%02d_r%02d' ${tp_name} ${num_rewirings} ${rewired_topology_number}`
    echo $candidate_topology
    target_parameterisation_number=1
    while test ${target_parameterisation_number} -le ${num_target_parameterisations} ; do
      tp_c_name=`printf '%s%02d' ${tp_basename} ${target_parameterisation_number}`
      do_run ./netoptrew -s ${rndseed} -l TRUE -o ${offset} -R ${num_optimisations} -e ${equilibration_timesteps} -n ${tp_c_name} -c ${candidate_topology} -u ${distance_measurement} -L $logfile -T $transformerfile -g ${gradientfile}
      rndseed=`expr ${rndseed} + 1`
      target_parameterisation_number=`expr ${target_parameterisation_number} + 1`
    done
    rewired_topology_number=`expr ${rewired_topology_number} + 1`
  done
}


function write_cshscript_header ()
{
  scriptname=$1
  echo '#!/bin/csh' > $scriptname
  echo '#$ -q long.q' >> $scriptname
  echo '#$ -j y' >> $scriptname
  echo '#$ -m e' >> $scriptname
  echo '#$ -cwd' >> $scriptname
  echo >> $scriptname
  echo 'setenv PYTHONPATH ${HOME}/lib64/python' >> $scriptname
  echo "cd $PWD" >> $scriptname
  echo 'echo start time: `date`' >> $scriptname
}


function write_cshscript_footer ()
{
  echo 'echo success: `date`' >> $scriptname
}


function add_command ()
{
  scriptname=$1
  cmd="$2"
  echo "if ( { $cmd } ) then" >> $scriptname
  echo "  echo completed \`date\`" >> $scriptname
  echo "else" >> $scriptname
  echo "  echo ERROR" >> $scriptname
  echo "  exit 1" >> $scriptname
  echo "endif" >> $scriptname
  #if [[ $rewired_topology_number -eq $num_rewired_topologies ]]
  #  then
  #  echo "qsub ${c_t}.csh" >> $scriptname
  #fi
}


function optimise_numrewired_cluster ()
{
  for num_rewirings in ${num_rewirings_list} ; do
    scriptname=`printf '%s_w%02d.csh' ${candidate_topology_basename} ${num_rewirings}`
    write_cshscript_header $scriptname
    rewired_topology_number=1
    while test ${rewired_topology_number} -le ${num_rewired_topologies} ; do
      candidate_topology=`printf '%s_w%02d_r%02d' ${candidate_topology_basename} ${num_rewirings} ${rewired_topology_number}`
      c_t=`printf '%s_w%02d' ${c_t_b} ${num_rewirings}`
      add_command $scriptname "./netoptrewGold -s ${rndseed} -o ${simgenex} -x ${empiricaldata} -R ${num_optimisations} -c ${candidate_topology} -L $logfile -T $transformerfile -g ${gradientfile}"
      rndseed=`expr ${rndseed} + 1`
      rewired_topology_number=`expr ${rewired_topology_number} + 1`
    done
    if [[ $target_topology_number -eq 1 ]]
      then
      do_run qsub ${scriptname}
    fi
  done
}


function optimiseGold ()
{
  for gentype in ${gentype_list} ; do
    target_topology_number=1
    while test ${target_topology_number} -le ${num_target_topologies} ; do
      topology_name=`printf 'target_%s%02d' ${gentype} ${target_topology_number}`
      candidate_topology_basename=`printf 'candidate_%s%02d' ${gentype} ${target_topology_number}`
      c_t_b=`printf 'candidate_%s%02d' ${gentype} ${target_topology_number}`
      optimise_numrewired_cluster
      target_topology_number=`expr ${target_topology_number} + 1`
    done
  done
}

function optimise ()
{
  for gentype in ${gentype_list} ; do
    target_topology_number=1
    while test ${target_topology_number} -le ${num_target_topologies} ; do
      target_parameterisation_number=1
      while test ${target_parameterisation_number} -le ${num_target_parameterisations} ; do
        topology_name=`printf 'target_%s%02d_p%02d' ${gentype} ${target_topology_number} ${target_parameterisation_number} `
        candidate_topology_basename=`printf 'candidate_%s%02d_p%02d' ${gentype} ${target_topology_number} ${target_parameterisation_number}`
        t=`expr $target_parameterisation_number + 1`
        c_t_b=`printf 'candidate_%s%02d_p%02d' ${gentype} ${target_topology_number} ${t}`
        optimise_numrewired_cluster
        target_parameterisation_number=`expr ${target_parameterisation_number} + 1`
      done
      target_topology_number=`expr ${target_topology_number} + 1`
    done
  done
}



function maketable () # create transsys program
{
  for gentype in ${gentype_list} ; do
    target_topology_number=1
    while test ${target_topology_number} -le ${num_target_topologies} ; do
      name=""
      target_parameterisation_number=1
        candidate_topology_basename=`printf 'candidate_%s%02d' ${gentype} ${target_topology_number}`
        for num_rewirings in ${num_rewirings_list} ; do
          rewired_topology_number=1
            while test ${rewired_topology_number} -le ${num_rewired_topologies} ; do
              candidate_topology=`printf '%s_w%02d_r%02d' ${candidate_topology_basename} ${num_rewirings} ${rewired_topology_number}`
              candidate_topology_logo=`printf '%s_logo.txt' ${candidate_topology}`
	      cp $candidate_topology_logo $candidate_topology_logo'.bk'
              sed 's/rst/'$target_topology_number'\t'$num_rewirings'\t'$rewired_topology_number'\t/g' $candidate_topology_logo'.bk' > clean.txt
              mv clean.txt $candidate_topology_logo
              rm -f $candidate_topology_logo'.bk'
              name=`printf '%s %s ' $name $candidate_topology_logo ` 
              rewired_topology_number=`expr ${rewired_topology_number} + 1`
            done
      done
    fitnessname=`printf 'fitnesstable%02d' $target_topology_number` 
    target_topology_number=`expr ${target_topology_number} + 1`
    cat $name > 'tt.txt'
    sed '/restart/d' 'tt.txt' > $fitnessname'.txt'
    rm -f 'tt.txt'
    done
    #cat $name > 'tt.txt'
    #sed '/restart/d' 'tt.txt' > $fitnessname'.txt'
    #rm -f 'tt.txt'
  done
}


function count_rw () # create transsys program
{
  for gentype in ${gentype_list} ; do
    target_topology_number=1
    while test ${target_topology_number} -le ${num_target_topologies} ; do
      name=""
      target_parameterisation_number=1
      #while test ${target_parameterisation_number} -le ${num_target_parameterisations} ; do
        candidate_topology_basename=`printf 'candidate_%s%02d' ${gentype} ${target_topology_number}`
        for num_rewirings in ${num_rewirings_list} ; do
          rewired_topology_number=1
            while test ${rewired_topology_number} -le ${num_rewired_topologies} ; do
              target_topology=`printf '%s_w00_r%02d.tra' ${candidate_topology_basename} ${rewired_topology_number}`
              candidate_topology=`printf '%s_w%02d_r%02d.tra' ${candidate_topology_basename} ${num_rewirings} ${rewired_topology_number}`
              diff $target_topology $candidate_topology > tt.txt
              t=$(grep -c '>' tt.txt)
              echo "$num_rewirings\t$t" >> rw.txt 
              rewired_topology_number=`expr ${rewired_topology_number} + 1`
            done
	done
        #target_parameterisation_number=`expr ${target_parameterisation_number} + 1`
      #done
    target_topology_number=`expr ${target_topology_number} + 1`
    done
  done
}

#control parameters
num_target_topologies=1
num_target_parameterisations=1
num_rewirings_list='0 1 2 3 4 5 6 7 9 11 13 15 18 22 27 32 38 46 55 66 100 1000'
gentype_list='er'
num_rewired_topologies=10
num_optimisations=5
signal_to_noise=0
transformerfile=transformerfile.dat
logfile=logo
gradientfile=optspec.dat
simgenex=testDream.sgx
empiricaldata=procdata.tgs

# initial rndseed, incremented each time a rndseed parameter is required
rndseed=2

# run the show
checkpython
#count_rw
generate_candidate_programs
#optimise_numrewired
optimiseGold
#maketable

