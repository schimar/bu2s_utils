#! /bin/bash
# cd /home/mschilling/Desktop/gbs15/mstr/
for i in Run*/; do 
	bunzip2 $i/parameters.m.bz2; 
	~/flaxmans/bu2s/bu2s_utils/xtrc_params.py $i; 
	bzip2 $i/parameters.m;
done



