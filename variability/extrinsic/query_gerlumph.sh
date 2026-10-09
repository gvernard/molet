#!/bin/bash


out_path=$1

json_maps=`grep -o '^[^//]*' ${out_path}output/multiple_images.json`

gerlumph_maps='[]'
Nq=`echo $json_maps | /usr/bin/jq '. | length'`
for (( q=0; q<$Nq; q++ ))
do  
    k=`echo $json_maps | /usr/bin/jq ".[$q].k "`
    g=`echo $json_maps | /usr/bin/jq ".[$q].k "`
    s=`echo $json_maps | /usr/bin/jq ".[$q].k "`

    query_args="{\"k\":$k,\"g\":$g,\"s\":$s}"
    response=`curl -s -u "Giorgos:123micro" -X POST http://127.0.0.7:7777/api/query-single-map/ -H "Content-Type: application/json" -d $query_args`
    #echo $response
    
    check=`echo $response | /usr/bin/jq '.message' | /usr/bin/sed -e 's/^"//' -e 's/"$//' `
    if [[ $check == *"No map"* ]];
    then
	map_prop=`jq -n --arg id "none" --argjson k "0" --argjson g "0" --argjson s "0" --argjson ds "0" --argjson dkg "0" '$ARGS.named'`
    else
	id=`echo $response | /usr/bin/jq '.map.id'`
	kk=`echo $response | /usr/bin/jq '.map.k'`
	gg=`echo $response | /usr/bin/jq '.map.g'`
	ss=`echo $response | /usr/bin/jq '.map.s'`
	ds=`echo $response | /usr/bin/jq '.map.ds'`
	dkg=`echo $response | /usr/bin/jq '.map.dkg'`
	# resolution
	# width
	# Navg
	# mu_avg
	map_prop=`jq -n --arg id "$id" --argjson k "$kk" --argjson g "$gg" --argjson s "$ss" --argjson ds "$ds" --argjson dkg "$dkg" '$ARGS.named'`
    fi

    gerlumph_maps=`echo "$json_data" | jq --argjson map "$map_prop" '. += [$map]'`
done


echo $json_data > ${out_path}output/gerlumph_maps_2.json

