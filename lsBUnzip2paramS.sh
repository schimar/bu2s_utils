#! /bin/bash

for i in Run215*/; do 
	#echo $i
	id=$(echo $i | cut -f1 -d/)
	if ls ${i}parameters.m.bz2; then
		bunzip2 ${i}parameters.m.bz2; 
	fi
	~/flaxmans/bu2s/bu2s_utils/xtrct_params.py ${i}parameters.m $i params.txt; 
	#version=$(grep CodeVersion ${i}parameters.m | cut -f2 -d= | cut -f2 -d\') 
	#runVersion=$(echo $id $version)    #cut -f2 -d=; echo $version
	#echo $id	$version
	if ls ${i}parameters.m; then
		bzip2 ${i}parameters.m
	fi
done



