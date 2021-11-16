#!/bin/bash

# python my_synthetic.py -g 0 -e 200 -ld 16 -bs 8 -nl 2 -hd 32 -lr 0.0005 -nb 5 -ng 250
max_gamma=0
for ld in 8 16
do
  for bs in 2 4 
    do
      for nl in 2 4
      do
        for hd in 64 128
        do
	  for lr in 1e-04 5e-04
	  do
             gamma=$(python ../../../NASA/GAN/my_synthetic.py -g 0 -e 100 -ld $ld -bs $bs -nl $nl -hd $hd -lr $lr -nb 5 -ng 0 | grep "Gamma(Dx, Dz):" | awk -F: '{print $2}' | xargs)
             echo "gamma:" $gamma "ld:" $ld "bs:" $bs "nl:" $nl "hd:" $hd "lr:" $lr >> ./params.txt
             if [ $gamma -gt $max_gamma ]
             then
                echo "gamma:" $gamma "ld:" $ld "bs:" $bs "nl:" $nl "hd:" $hd "lr:" $lr >> ./best_params.txt
                max_gamma=$gamma
             fi
           done
        done
      done
    done
done
