#!/bin/csh
#$ -q long.q
#$ -m e
#$ -cwd

pwd

date

echo -- start time --

setenv PYTHONPATH $HOME/lib64/python
$HOME/makemodel-1.0/testonenet.sh -t z -p x -w b -s q -r 10 -l 1 -i 2 -L "logo" -P "$HOME/makemodel-1.0/parapl.dat" -T "$HOME/makemodel-1.0/transformerfile.dat" -O "$HOME/makemodel-1.0/optspec.dat" -o 0.01 -m 2


echo -- end time --

date
