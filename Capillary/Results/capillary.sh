#!/bin/bash

for a in `seq 1 16`
do
  echo $a
  mkdir $a
  cp capillary.pbs $a/capillary.pbsold
  cp main.out $a/
  cd $a
  mkdir tmp2
  
  sed 's/TOCHANGE1/'"$a"'/g' capillary.pbsold > capillary.pbs
  rm -f capillary.pbsold

  qsub capillary.pbs
  #rm main.out
  #rm capillary.pbs 
  cd ..
done
exit 0
