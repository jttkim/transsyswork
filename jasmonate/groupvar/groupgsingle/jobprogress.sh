#!/bin/sh

for chainhead in `ls *squaresum.csh` ; do
  chainname=`echo $chainhead | sed -e 's/_squaresum.*//'`
  chaintail=${chainname}_correlation_logratio.csh
  chaintail_log=`ls ${chaintail}.o*`
  if test -f "$chaintail_log" ; then
    if grep -q 'chain complete' $chaintail_log ; then
      echo $chainname complete
    else
      echo $chainname running / aborted
    fi
  else
    echo $chainname not started
  fi
done

