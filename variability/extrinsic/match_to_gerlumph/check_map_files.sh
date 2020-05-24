#!/bin/bash


map_path=$1
out_path=$2

json_maps=`grep -o '^[^//]*' ${out_path}output/gerlumph_maps.json`

missing=()
Nq=`echo $json_maps | jq '. | length'`
for (( q=0; q<$Nq; q++ ))
do  
    id=`echo $json_maps | jq ".[$q].id " | sed -e 's/^"//' -e 's/"$//'`
    if [ $id != "none" ]
    then
	mydir=${map_path}$id
	if [ ! -d "$mydir" ]
	then
	    missing+=("$id")
	fi	
    fi
done




urlencode() {
    # urlencode <string>
    old_lc_collate=$LC_COLLATE
    LC_COLLATE=C
    
    local length="${#1}"
    for (( i = 0; i < length; i++ )); do
	local c="${1:i:1}"
	case $c in
            [a-zA-Z0-9.~_-]) printf "$c" ;;
            *) printf '%%%02X' "'$c" ;;
	esac
    done
    
    LC_COLLATE=$old_lc_collate
}


if [ ${#missing[@]} -ne 0 ]
then
    get_args=""
    for id in "${missing[@]}"
    do
	get_args=${get_args}"ids[]="$id"&"
    done
    get_args=${get_args::-1}
    
    printf "ATTENTION: missing GERLUMPH maps from '%s' !\n" $map_path 1>&2
    printf "Click the following link to download them:" 1>&2
    printf "\n\n" 1>&2
    printf "      %s" "https://gerlumph.swin.edu.au/inc/generic/put_to_cart.php?"$get_args 1>&2
    printf "\n\n" 1>&2
    printf "and place them at: %s " $map_path 1>&2
fi


