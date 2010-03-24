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


function maketp () # create transsys program
{
   for (( i=1; i<=1; i++)) 
   do
     for (( j=1; j<=5; j++)) 
     do
       name=""
       for k in 0 1 2 3 4 5 6 7 9 11 13 15 18 22 27 32 38 46 55 66
       do
         tname=`printf '%s%02d%02d%02dlogo.txt' $net_name $i $j $k`
         t=$(ls -ls $tname | awk '{print$1}')
         if [ $t -eq 0 ]
           then
           t_p=`printf '%02d%02d%02d;  ' $i $j $k` 
           namenull=`printf '%s %s \n' $namenull $tname ` 
           t_f=`printf '%s %s ' $t_f $t_p ` 
         fi
         fname=`printf '%s%02d%02dopt.txt' $net_name $i $j`
         name=`printf '%s %s ' $name $tname ` 
       done
       cat $name > $fname
       sed '/net_top/d' $fname > clean.txt
       mv clean.txt $fname
     done
   done
   echo "None existing data tp/rp/rw/:" $t_f
}

function MakeTable()
{



}

function makeonetp () # create transsys program
{

   name=""
   fname=`printf '%sN%s%s.txt' $net_name $noise_ratio $net_type`
   for (( i=1; i<=1; i++)) 
   do
     for (( j=1; j<=5; j++)) 
     do
       tname=`printf '%s%02d%02dopt.txt' $net_name $i $j`
       name=`printf '%s %s ' $name $tname ` 
     done
   done
   cat $name > $fname
}


while getopts n:o:s:d opt
do
  case "$opt" in
    n) net_name="$OPTARG";;
    o) net_type="$OPTARG";;
    s) noise_ratio="$OPTARG";;
    d) isdef=1;;
    \?) help_ani;;
  esac
done

net_name=yn30e60
net_type=ER
noise_ratio=0

maketp
makeonetp
