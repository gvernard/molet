#!/bin/bash

infile=$1
in_path=$2
out_path=$3


myinput=`grep -o '^[^//]*' $infile`

band_name=`echo $myinput | jq '.instrument.bands[0].name'`



#### Find and set number of intrinsic light curves (Nin)
lc_in_type=`echo $myinput | jq '.point_source.variability.intrinsic.type' | sed -e 's/^"//' -e 's/"$//'`
if [ $lc_in_type = "custom" ]
then
    lc_in_file=${in_path}"input_files/intrinsic_light_curves.json"
else
    lc_in_file=${out_path}"output/intrinsic_light_curves.json"
fi
Nin=`cat $lc_in_file | jq ".$band_name | length"`




#### Find and set number of extrinsic light curves (Nex)
lc_ex_type=`echo $myinput | jq '.point_source.variability.extrinsic.type' | sed -e 's/^"//' -e 's/"$//'`
if [ $lc_ex_type = "custom" ]
then
    lc_ex_file=${in_path}"input_files/extrinsic_light_curves.json"
else
    lc_ex_file=${out_path}"output/extrinsic_light_curves.json"
fi
Nq=`cat $lc_ex_file | jq ". | length"`

for (( q=0; q<$Nq; q++ ))
do  
    len=`cat $lc_ex_file | jq ".[$q].$band_name | length"`
    if [ $len -gt 0 ]
    then
	Nex=$len
    fi
done



#### Loop and create mock directories, remove beforehand if they exist
for (( i=0; i<$Nin; i++ ))
do
    for (( j=0; j<$Nex; j++ ))
    do
	printf -v mock "mock_%04d_%04d" $i $j
	mydir=${out_path}$mock
	if [ -d "$mydir" ]
	then
	    rm -r $mydir
	fi
	mkdir $mydir
    done
done
