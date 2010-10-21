#!/bin/bash

for a in `seq 1 4`
do
  echo $a
  b=49*$a
  echo $b
  #mkdir Omega${a*49+2}
  #cp main.out Omega$a/
  #cp athena_palettes/LeeLoo.txt Omega$a/
  
  #cd Omega$a
 
  #for b in `seq 20 2 30`
  #do
  #  mkdir R$b
  #  mkdir R$b/tmp
  #  mkdir R$b/athena_palettes

  #  cp LeeLoo.txt R$b/athena_palettes/
  #  cp main.out R$b/
  #  cd R$b
  #  ./main.out $a $b
  #  rm main.out
  #  rm -rf athena_palettes
  #  cd ..
  #done
    
  #rm main.out
  #rm LeeLoo.txt
  #cd ..
done
exit 0