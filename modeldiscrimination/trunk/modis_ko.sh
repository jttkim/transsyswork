#!/bin/sh

function do_run ()
{
  # echo $*
  if echo $* ; then
    true
  else
    echo '*** ERROR ***' >& 2
    exit 1
  fi
}



# naming scheme:
# target topology: target_<gentype>##.tra
# target program:  target_<gentype>##_p##.tra


function generate_target_programs ()
{
  # generate target topologies
  # note: rndseed for topology generation is in transsysgen files
  # jtk: is parapl.dat the right file? which one is the ER file?
  for gentype in er pl ; do
    do_run transsysrandomprogram -t target_${gentype} -m ${num_target_topologies} -p para${gentype}.dat
  done
  # parameterise target topologies
  for gentype in er pl ; do
    target_topology_number=1
    while test ${target_topology_number} -le ${num_target_topologies} ; do
      topology_name=`printf 'target_%s%02d' ${gentype} ${target_topology_number}`
      topology_file=${topology_name}.tra
      tp_basename=${topology_name}_p
      do_run transsysreparam -T ${transformerfile} -r ${rndseed} -p ${num_target_parameterisations} -n ${tp_basename} ${topology_file}
      rndseed=`expr ${rndseed} + 1`
      target_topology_number=`expr ${target_topology_number} + 1`
    done
  done
}


function generate_target_expressionsets ()
{
  for gentype in er pl ; do
    target_topology_number=1
    while test ${target_topology_number} -le ${num_target_topologies} ; do
      topology_name=`printf 'target_%s%02d' ${gentype} ${target_topology_number}`
      topology_file=${topology_name}.tra
      tp_basename=${topology_name}_p
      target_parameterisation_number=1
      while test ${target_parameterisation_number} -le ${num_target_parameterisations} ; do
	tp_name=`printf '%s%02d' ${tp_basename} ${target_parameterisation_number}`
	do_run echo generate expression data from ${tp_name} with equilibration ${equilibration_timesteps}
        target_parameterisation_number=`expr ${target_parameterisation_number} + 1`
      done
      target_topology_number=`expr ${target_topology_number} + 1`
    done
  done
}


function generate_candidate_programs ()
{
  for gentype in er pl ; do
    target_topology_number=1
    while test ${target_topology_number} -le ${num_target_topologies} ; do
      target_topology_name=`printf 'target_%s%02d' ${gentype} ${target_topology_number}`
      target_topology_file=${target_topology_name}.tra
      candidate_topology_basename=`printf 'candidate_%s%02d' ${gentype} ${target_topology_number}`
      do_run transsysrewire -t ${candidate_topology_basename} -v ${distinct_num_rewirings} -r ${num_rewired_topologies} ${target_topology_file} --should-use-rndseed ${rndseed}
      rndseed=`expr ${rndseed} + 1`
      target_topology_number=`expr ${target_topology_number} + 1`
    done
  done
}


function optimise ()
{
  for gentype in er pl ; do
    target_topology_number=1
    while test ${target_topology_number} -le ${num_target_topologies} ; do
      topology_name=`printf 'target_%s%02d' ${gentype} ${target_topology_number}`
      tp_basename=${topology_name}_p
      candidate_topology_basename=`printf 'candidate_%s%02d' ${gentype} ${target_topology_number}`
      for num_rewirings in ${num_rewirings_list} ; do
	rewired_topology_number=1
	while test ${rewired_topology_number} -le ${num_rewired_topologies} ; do
	  candidate_topology=`printf '%s_%02d%02d' ${candidate_topology_basename} ${num_rewirings} ${rewired_topology_number}`
	  target_parameterisation_number=1
	  while test ${target_parameterisation_number} -le ${num_target_parameterisations} ; do
	    tp_name=`printf '%s%02d' ${tp_basename} ${target_parameterisation_number}`
	    do_run echo optimise ${candidate_topology} against expdata from ${tp_name} with ${num_optimisations} restarts using transformers from ${transformerfile} with equilibration ${equilibration_timesteps}, rndseed = ${rndseed}
	    rndseed=`expr ${rndseed} + 1`
	    target_parameterisation_number=`expr ${target_parameterisation_number} + 1`
	  done
	  rewired_topology_number=`expr ${rewired_topology_number} + 1`
	done
        rewired_topology_number=`expr ${rewired_topology_number} + 1`
      done
      target_topology_number=`expr ${target_topology_number} + 1`
    done
  done
}


#control parameters
num_target_topologies=10
num_target_parameterisations=5
distinct_num_rewirings=24
# jtk: redundant to distinct_num_rewirings -- remove list generation from transsysrewire?
num_rewirings_list='0 1 2 3 4 5 6 7 9 11 13 15 18 22 27 32 38 46 55 66'
num_rewired_topologies=5
num_optimisations=5
equilibration_timesteps=100
transformerfile=transformerfile.dat

# initial rndseed, incremented each time a rndseed parameter is required
rndseed=2

# run the show
generate_target_programs
generate_target_expressionsets
generate_candidate_programs
optimise

