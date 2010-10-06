#!/bin/bash

for a in `seq -10 2 10`
do
  echo $a
  mkdir grad$a
  cp wall_gradient.pbs grad$a/wall_gradient.pbsold
  cp main.out grad$a/
  cd grad$a
  mkdir tmp2
  
  sed 's/TOCHANGE1/'"$a"'/g' wall_gradient.pbsold > wall_gradient.pbs
  rm -f wall_gradient.pbsold
  nohup ./main.out $a &
  #qsub wall_gradient.pbs
  #rm main.out
  #rm capillary.pbs 
  cd ..
done
exit 0
