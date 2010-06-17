#!/bin/sh

function do_run ()
{
  echo $*
  if $* ; then
    true
  else
    echo '*** ERROR ***' >&2
    exit 1
  fi
}


#do_run ./transsysrandomprogram -n test -m 1 -p rtptest.dat
#do_run mv test01.tra test_targettopology.tra
#do_run ./transsysreparam -T transformerfile.dat -s 1 -p 1 -n test test_targettopology.tra
#do_run mv test01.tra test_target.tra
do_run ./transsyswritesimsetOF -o jinjut_obj.txt -s 2 -N 0.1 jinjut_toy_01.tra test_target dummy.txt
#do_run ./transsysrewire -w 1 -n test_candidate -r 1 -s 3 test_targettopology.tra
#do_run ./transsysrewire -w 3 -n test_candidate -r 1 -s 3 test_targettopology.tra
#do_run ./transsysrewire -w 4 -n test_candidate -r 1 -s 3 test_targettopology.tra

#do_run ./netopt -x test_target_expr.txt -o jinjut_obfunctoy.txt -R 5 -g optspec.dat -T transformerfile.dat -s 1 -c jinjut_toy.tra t1.tzxt t2.txt

#do_run ./netoptrew -l -o 0.01 -s 4 -R 5 -e 100 -n test_target -c test_targettopology -g optspec.dat -T transformerfile.dat -u correlation -L test_rew0 testlog_rew0.txt testfinalparam_rew0.txt
#do_run ./netoptrew -l -o 1e-10 -s 4 -R 1 -e 100 -n test_target -c test_candidate_w01_r01 -g optspec.dat -T transformerfile.dat -u correlation -L test_rew1 testlog_rew1.txt testfinalparam_rew1.txt
#do_run ./netoptrew -l -o 0.01 -s 4 -R 1 -e 100 -n test_target -c test_candidate_w03_r01 -g optspec.dat -T transformerfile.dat -u correlation -L test_rew1 testlog_rew1.txt testfinalparam_rew1.txt
#do_run ./netoptrew -l -o 0.01 -s 4 -R 1 -e 100 -n test_target -c test_candidate_w10_r01 -g optspec.dat -T transformerfile.dat -u correlation -L test_rew1 testlog_rew1.txt testfinalparam_rew1.txt
