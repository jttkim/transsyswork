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


do_run ./transsyswritesimsetOF -o simgenex_temp01.txt test_temp01.tra temp01 dummy.txt
do_run ./netopt -x temp01_expr.txt -o simgenex_temp01.txt -R 1 -g optspec.dat -T transformerfile.dat -L logo -s 1 -c test_temp01 t1.txt t2.txt

