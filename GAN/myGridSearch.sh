#!/bin/bash

# -e 200 -ld 8 -bs 1 -nl 2 -hd 64 -lr 5e-03 -nb 5 -ng 0  -ie ./Proj2_Normalized_Counts.csv -im ./all_metadata_Proj2.csv -cd checkpoints/ -gf top-liver-genes.txt
max_gamma=0
epochs=100
re='^[0-9]+([.][0-9]+)?$'
for ld in 7 8 9 
do
  for bs in 1 16
    do
      for nl in 1 2 3 
      do
        for hd in 128 192
        do
	  for lr in 1e-03 5e-02
	  do
             #python ../../../NASA/GAN/my_synthetic.py -g 0 -e 100 -ld $ld -bs $bs -nl $nl -hd $hd -lr $lr -nb 5 -ng 0 
	     echo "new job: ld=$ld bs=$bs nl=$nl hd=$hd lr=$lr"
             gamma=$(python ../../../NASA/GAN/my_synthetic.py -g 0 -e $epochs -ld $ld -bs $bs -nl $nl -hd $hd -lr $lr -nb 5 -ng 0 -s 23 2>/dev/null | grep "Gamma(Dx, Dz):" | awk -F: '{print $2}' | xargs)
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
