#!/bin/bash

# python my_synthetic.py -g 0 -e 200 -ld 16 -bs 8 -nl 2 -hd 32 -lr 0.0005 -nb 5 -ng 250
max_gamma=0
epochs=30
re='^[0-9]+([.][0-9]+)?$'
for ld in 2 4
do
  for bs in 4 8
    do
      for nl in 4 8
      do
        for hd in 128 256
        do
	  for lr in 5e-04 1e-03
	  do
             #python ../../../NASA/GAN/my_synthetic.py -g 0 -e 100 -ld $ld -bs $bs -nl $nl -hd $hd -lr $lr -nb 5 -ng 0 
	     echo "new job: ld=$ld bs=$bs nl=$nl hd=$hd lr=$lr"
             gamma=$(python ../../../NASA/GAN/my_synthetic.py -g 0 -e $epochs -ld $ld -bs $bs -nl $nl -hd $hd -lr $lr -nb 5 -ng 0 2>/dev/null | grep "Gamma(Dx, Dz):" | awk -F: '{print $2}' | xargs)
	     if ! [[ $gamma =~ $re ]]
	     then
	        rm -rf checkpoints
	        mkdir checkpoints
		echo "got a nan ... continuing to next"
                continue 
 	     fi
	     echo "got a number $gamma ... seeing if it's max"
	     gamma_ten=$(echo $gamma | cut -d. -f2)
             echo "gamma:" $gamma "ld:" $ld "bs:" $bs "nl:" $nl "hd:" $hd "lr:" $lr >> ./params.txt
             if [ $gamma_ten -gt $max_gamma ]
             then
		echo "found max $gamma"
                echo "gamma:" $gamma "ld:" $ld "bs:" $bs "nl:" $nl "hd:" $hd "lr:" $lr >> ./best_params.txt
		cp checkpoints/models/gen_liver.h5 MODELS/gamma_${gamma}_ld_${ld}_bs_${bs}_nl_${nl}_hd_${hd}_lr_${lr}.h5
                max_gamma=$gamma_ten
             fi
           done
        done
      done
    done
done
