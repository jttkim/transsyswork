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


#do_run ./transsyswritesimsetOF -o jasmonate_model_v1.txt -s 2 -N 0 jasmonate.tra testkojas testkojas_trace.txt
#do_run ./transsyswritesimsetOF -o jasmonate_model_snapshot.txt -s 2 -N 0 jasmonate_m00.tra jasmonatem00 testkojas_trace.txt
#do_run ./transsyswritesimsetOF -o jasmonate_model_snapshot.txt -s 2 -N 0 jasmonate_m01.tra jasmonatem01 testkojas_trace.txt
#do_run ./transsyswritesimsetOF -o jasmonate_model_snapshot.txt -s 2 -N 0 jasmonate_m02.tra jasmonatem02 testkojas_trace.txt
#do_run ./transsyswritesimsetOF -o jasmonate_model_snapshot.txt -s 2 -N 0 jasmonate_m03.tra jasmonatem03 testkojas_trace.txt
#do_run ./netopt -x testkojas_expr.txt -o jasmonate_model_snapshot.txt -R 5 -g optspec.dat -T transformerfile.dat -L logo -s 1 -c jasmonate.tra t1.txt t2.txt
do_run ./netopt -x test_target_expr_temp01.txt -o simgenex_temp01.txt -R 5 -g optspec.dat -T transformerfile.dat -L logo -s 1 -c test_temp01 t1.txt t2.txt
#do_run ./netopt -x testkojas_expr.txt -o jasmonate_model_snapshot.txt -R 5 -g optspec.dat -T transformerfile.dat -L logo -s 1 -c jasmonate_m01.tra t1.txt t2.txt
#do_run ./netopt -x testkojas_expr.txt -o jasmonate_model_snapshot.txt -R 5 -g optspec.dat -T transformerfile.dat -L logo -s 1 -c jasmonate_m02.tra t1.txt t2.txt
#do_run ./netopt -x testkojas_expr.txt -o jasmonate_model_snapshot.txt -R 5 -g optspec.dat -T transformerfile.dat -L logo -s 1 -c jasmonate_m03.tra t1.txt t2.txt
