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


do_run ./transsyswritesimsetOF -o sgx.txt sgx.tra out.txt
do_run ./netopt -x out.txt -o sgx.txt -R 1 -g optspec.dat -T transformerfile.dat -L logo -s 1 -c sgx t1.txt t2.txt

