import re
import os
import sys
import json
from subprocess import Popen,PIPE


def myprocess(msg,cmd_list,log_file):
    print("%-100s" % msg,end="",flush=True)
    p = Popen(cmd_list,universal_newlines=True,stdout=PIPE,stderr=PIPE)
    out,err = p.communicate()
    if len(err) != 0:
        sys.exit("\n===> Failed with error: \n%s" % (str(err))+"\n===> Use the following command on its own to debug:\n"+" ".join(cmd_list))
    if p.returncode != 0:
        sys.exit("\n===> Failed with unspecified error\n===> Use the following command on its own to debug:\n"+" ".join(cmd_list))
    log_file.write("COMMAND: "+" ".join(cmd_list)+"\n")
    log_file.write("===="*50+"\n")
    if len(out) != 0:
        log_file.write(out)
    print("%-10s" % "...done",flush=True)




infile = sys.argv[1]
infile = os.path.abspath(infile) # make sure this is the absolute path

dum  = infile.split("/")
in_path = "/".join(dum[:-1]) + "/"

# Check if input_files directory exists (it needs to contain at least the psf files)
if not os.path.isdir(in_path+"input_files"):
    print("Input files must be in a directory named 'input_files', at the same path as the 'input_molet.json' file!")
    sys.exit()

if len(sys.argv) > 2:
    out_path = sys.argv[2]
    if out_path[-1] != "/":
        out_path += "/"
else:
    out_path = in_path    
    
molet_home = os.path.dirname(os.path.abspath(__file__)) + "/"


# Get map path
p = Popen([molet_home+"variability/extrinsic/get_map_path/bin/get_map_path"],universal_newlines=True,stdout=PIPE,stderr=PIPE)
out,err = p.communicate()
p.kill()
map_path = out.strip()


f          = open(infile,'r')
input_str  = f.read()
input_str  = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str  = re.sub(re.compile("//.*?\n" ),"",input_str)
json_in    = json.loads(input_str)

if os.path.isdir(out_path+"output") == False:
    os.mkdir(out_path+"output")

log_file = open(out_path+"log.txt",'w')


    

# Step 1:
# Get angular diameter distances
####################################################################################
cmd_list = [
    molet_home+"cosmology/angular_diameter_distances/bin/angular_diameter_distances",
    infile,
    out_path
]
#print(" ".join(cmd_list))
myprocess("Getting angular diameter distances...",cmd_list,log_file)



# Step 2:
# Get extended lensed images of the source
####################################################################################
cmd_list = [
    molet_home+"lensed_extended_source/vkl_fproject/bin/fproject",
    infile,
    out_path
]
#print(" ".join(cmd_list))
myprocess("Getting extended source lensed features...",cmd_list,log_file)



# Intermediate step:
# Get point source images, locations are needed for the following
####################################################################################
if "point_source" in json_in:
    cmd_list = [
        molet_home+"lensed_point_source/vkl_point_source/bin/point_source",
        infile,
        out_path
    ]
    #print(" ".join(cmd_list))
    myprocess("Getting point-like source lensed images...",cmd_list,log_file)


    
# Step 3:
# Get light profile of the lens (and compact matter if required)
####################################################################################
cmd_list = [
    molet_home+"lens_light_mass/vkl_llm/bin/llm",
    infile,
    out_path
]
#print(" ".join(cmd_list))
myprocess("Getting light profile of the lens...",cmd_list,log_file)



# Intermediate step:
# Get variability
####################################################################################
if "point_source" in json_in:

    # Intrinsic
    if json_in["point_source"]["variability"]["intrinsic"]["type"] != "custom":
        dummm = "yo"
    
    # Extrinsic
    if json_in["point_source"]["variability"]["extrinsic"]["type"] != "custom":
        cmd_list = [
            "python3",
            molet_home+"variability/extrinsic/match_to_gerlumph/match_to_gerlumph.py",
            molet_home+"data/gerlumph.db",
            map_path,
            out_path
        ]
        #print(" ".join(cmd_list))
        myprocess("Matching macro-images to GERLUMPH maps...",cmd_list,log_file)
    
        cmd_list = [
            molet_home+"variability/extrinsic/gerlumph_moving_source/bin/moving_source",
            infile,
            out_path
        ]
        #print(" ".join(cmd_list))
        myprocess("Getting microlensing variability for each image...",cmd_list,log_file)



# Step 4:
# Combine different light components
####################################################################################
# Create output directories if necessary
if "point_source" in json_in:
    cmd_list = [
        "python3",
        molet_home+"combined_light/setup_dirs.py",
        infile,
        in_path,
        out_path
    ]
    #print(" ".join(cmd_list))
    myprocess("Mock output directories created...",cmd_list,log_file)

# Combine light
cmd_list = [
    molet_home+"combined_light/bin/combine_light",
    infile,
    in_path,
    out_path
]
#print(" ".join(cmd_list))
myprocess("Combining light components and including instrumental effects...",cmd_list,log_file)
    

log_file.close()
print("\nCompleted successfully!\n")
print("Output in: ",out_path)
