#!/bin/bash


myprocess () {
    msg=$1
    cmd_list=$2
    log_file=$3

    
    echo "COMMAND: $cmd_list" >> $log_file
    dum="=========================================="
    echo $dum$dum$dum$dum$dum >> $log_file

    printf "%-100s" "$msg"

    error_file=$(mktemp)
    out=$($cmd_list 2>$error_file)
    exit_code=$?
    err=$(< $error_file)
    rm $error_file
    
    if [ ! -z "$err" ]
    then
        printf "\n===> Failed with error: \n%s" "$err"
	printf "\n===> Use the following command on its own to debug:\n"
	printf "$cmd_list\n"
	exit
    fi

    if [ "$exit_code" -ne "0" ]
    then
	printf "\n===> Failed with unspecified error\n===> Use the following command on its own to debug:\n"
	printf "$cmd_list\n"
    fi
           
    if [ ! -z "$out" ]
    then
       echo "$out" >> $log_file
    fi
       
    printf "%-10s\n" "...done"
}




infile=$1
infile=`realpath $infile`
injson=`grep -o '^[^//]*' $infile`
in_path=`dirname $infile`"/"
molet_home=`pwd`"/"


# Check if input_files directory exists (it needs to contain at least the intrinsic variability file)
if [ ! -d ${in_path}"input_files" ]
then
    printf "Input files must be in a directory named 'input_files', at the same path as the 'molet_input.json' file!\n"
    exit
fi


# Check instrument compatibility
available=($(ls -d ${molet_home}instrument_modules/*))
Ninstruments=`echo $injson | jq '.instruments | length'`
instruments=()
for (( b=0; b<$Ninstruments; b++ ))
do
    name=`echo $injson | jq ".instruments[$b].name" | sed -e 's/^"//' -e 's/"$//'`
    instruments+=( $name )
done
# Check if instrument name exists in the modules
for (( b=0; b<$Ninstruments; b++ ))
do
    check=false
    for (( i=0; i<${#available[@]}; i++ ))
    do  
	avail=`basename ${available[$i]}`
	if [ $avail == ${instruments[$b]} ]
	then
	    check=true
	    break
	fi
    done
    if [ "$check" = false ]
    then
	printf "Instrument named \"${instruments[$b]}\" does not exist!\n"
	exit
    fi
done
# Check if instruments match the INTRINSIC light curve files (one file per instrument)
for (( b=0; b<$Ninstruments; b++ ))
do  
    if [ ! -f ${in_path}"input_files/"${instruments[$b]}"_LC_intrinsic.json" ]
    then
	printf "Input INTRINSIC light curves don't exist for instrument \"${instruments[$b]}\"!\n"
	exit
    fi
done
# Check if instruments match the EXTRINSIC light curve files (one file per instrument)
extrinsic=`echo $injson | jq ".point_source.variability.extrinsic.type" | sed -e 's/^"//' -e 's/"$//'`
if [ $extrinsic = custom ]
then
    for (( b=0; b<$Ninstruments; b++ ))
    do  
	if [ ! -f ${in_path}"input_files/"${instruments[$b]}"_LC_extrinsic.json" ]
	then
	    printf "Input EXTRINSIC light curves don't exist for instrument \"${instruments[$b]}\"!\n"
	    exit
	fi
    done
fi


# Check if optional output path argument is present
if [ $# -eq 2 ]
then
    out_path=$2
    i=$((${#out_path}-1))
    if [ "${out_path:$i:1}" != "/" ]
    then
	out_path=${out_path}"/"
    fi
else
    out_path=$in_path
fi
if [ ! -d ${out_path}"output" ]
then
    mkdir ${out_path}"output"
fi
log_file=${out_path}"output/log.txt"


# Get map path
map_path=`${molet_home}"variability/extrinsic/get_map_path/bin/get_map_path"`




    
# Step 1:
# Get angular diameter distances
####################################################################################
msg="Getting angular diameter distances..."
cmd=$molet_home"cosmology/angular_diameter_distances/bin/angular_diameter_distances "$infile" "$out_path
myprocess "$msg" "$cmd" "$log_file"



# Step 2:
# Get extended lensed images of the source
####################################################################################
msg="Getting extended source lensed features..."
cmd=$molet_home"lensed_extended_source/vkl_fproject/bin/fproject "$infile" "$in_path" "$out_path
myprocess "$msg" "$cmd" "$log_file"



# Intermediate step:
# Get point source images, their locations are needed for the following
####################################################################################
check=`echo $injson | jq '. | select(.point_source)'`
if [ ! -z "${check}" ]
then
    msg="Getting point-like source lensed images..."
    cmd=$molet_home"lensed_point_source/vkl_point_source/bin/point_source "$infile" "$in_path" "$out_path
    myprocess "$msg" "$cmd" "$log_file"
fi



# Step 3:
# Get light profile of the lens (and compact matter if required)
####################################################################################
msg="Getting light profile of the lens..."
cmd=$molet_home"lens_light_mass/vkl_llm/bin/llm "$infile" "$in_path" "$out_path
myprocess "$msg" "$cmd" "$log_file"


    
# Intermediate step:
# Get variability
####################################################################################
check=`echo $injson | jq '. | select(.point_source)'`
if [ ! -z "${check}" ]
then
    # Extrinsic
    ex_type=`echo $injson | jq '.point_source.variability.extrinsic.type' | sed -e 's/^"//' -e 's/"$//'`
    if [ $ex_type != "custom" ]
    then
	msg="Matching macro-images to GERLUMPH maps..."
	cmd=$molet_home"variability/extrinsic/match_to_gerlumph/bin/match_to_gerlumph "$molet_home"data/gerlumph.db "$out_path
	myprocess "$msg" "$cmd" "$log_file"

	msg="Checking if GERLUMPH maps exist locally..."
	cmd=$molet_home"variability/extrinsic/match_to_gerlumph/check_map_files.sh "$map_path" "$out_path
	myprocess "$msg" "$cmd" "$log_file"

	if [ $ex_type = moving_disc ]
	then
	    msg="Getting 'moving_disc' microlensing variability for each image..."
	    cmd=$molet_home"variability/extrinsic/moving_disc/bin/moving_disc "$infile" "$out_path
	    myprocess "$msg" "$cmd" "$log_file"
	elif [ $ex_type = expanding_supernova ]
	then
	    msg="Getting 'expanding_supernova' microlensing variability for each image..."
	    cmd=$molet_home"variability/extrinsic/expanding_supernova/bin/expanding_supernova "$infile" "$out_path
	    myprocess "$msg" "$cmd" "$log_file"
	fi  
    fi    
fi


# Step 4:
# Combine different light components
####################################################################################
# Create output directories if necessary
check=`echo $injson | jq '. | select(.point_source)'`
if [ ! -z "${check}" ]
then
   msg="Mock output directories created..."
   cmd=$molet_home"combined_light/setup_dirs.sh "$infile" "$in_path" "$out_path
   myprocess "$msg" "$cmd" "$log_file"
fi


# Combine light
msg="Combining light components and including instrumental effects..."
cmd=$molet_home"combined_light/bin/combine_light "$infile" "$in_path" "$out_path
myprocess "$msg" "$cmd" "$log_file"
    

printf "\nCompleted successfully!\n\n"
printf "Output in: %s\n" $out_path
