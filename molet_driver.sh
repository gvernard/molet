#!/bin/bash


myprocess () {
    msg=$1
    cmd_list=$2
    log_file=$3

    
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
       
    echo "COMMAND: $cmd_list" >> $log_file
    dum="=========================================="
    echo $dum$dum$dum$dum$dum >> $log_file
    
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
    printf "Input files must be in a directory named 'input_files', at the same path as the 'input_molet.json' file!\n"
    exit
fi

# Check if instruments match the intrinsic and extrinsic light curve files (one file per instrument)


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
cmd=$molet_home"lensed_extended_source/vkl_fproject/bin/fproject "$infile" "$out_path
myprocess "$msg" "$cmd" "$log_file"



# Intermediate step:
# Get point source images, their locations are needed for the following
####################################################################################
check=`echo $injson | jq '. | select(.point_source)'`
if [ ! -z "${check}" ]
then
    msg="Getting point-like source lensed images..."
    cmd=$molet_home"lensed_point_source/vkl_point_source/bin/point_source "$infile" "$out_path
    myprocess "$msg" "$cmd" "$log_file"
fi



# Step 3:
# Get light profile of the lens (and compact matter if required)
####################################################################################
msg="Getting light profile of the lens..."
cmd=$molet_home"lens_light_mass/vkl_llm/bin/llm "$infile" "$out_path
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
	
	msg="Getting microlensing variability for each image..."
	cmd=$molet_home"variability/extrinsic/gerlumph_moving_source/bin/moving_source "$infile" "$out_path
	myprocess "$msg" "$cmd" "$log_file"	
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
