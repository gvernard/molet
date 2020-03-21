import re
import os
import sys
import json
from subprocess import Popen,PIPE


def myprocess(msg,cmd_list):
    print("%-100s" % msg,end="",flush=True)
    p = Popen(cmd_list,universal_newlines=True,stdout=PIPE,stderr=PIPE)
    out,err = p.communicate()
    if len(err) != 0:
        sys.exit("\n===> Failed with error: \n%s" % (str(err)))
    print("%-10s" % "...done",flush=True)




infile = sys.argv[1]
infile = os.path.abspath(infile) # make sure this is the absolute path

dum  = infile.split("/")
path = "/".join(dum[:-1]) + "/"
name = dum[-1]

molet_home = "/home/george/myCodes/molet/"
f          = open(infile,'r')
input_str  = f.read()
input_str  = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str  = re.sub(re.compile("//.*?\n" ),"",input_str)
json_in    = json.loads(input_str)






# Step 1:
# Get angular diameter distances
####################################################################################
cmd_list = [
    "python",
    molet_home+"cosmology/angular_diameter_distances.py",
    infile,
    path
]
#print(" ".join(cmd_list))
myprocess("Getting angular diameter distances...",cmd_list)




# Step 2:
# Get extended lensed images of the source
####################################################################################
cmd_list = [
    molet_home+"lensed_extended_source/vkl_fproject/bin/fproject",
    infile,
    path
]
#print(" ".join(cmd_list))
myprocess("Getting extended source lensed features...",cmd_list)


# Intermediate step:
# Get point source images, locations are needed for the following
####################################################################################
if "point_source" in json_in:
    cmd_list = [
        molet_home+"lensed_point_source/vkl_point_source/bin/point_source",
        infile,
        path
    ]
    #print(" ".join(cmd_list))
    myprocess("Getting point-like source lensed images...",cmd_list)


    
# Step 3:
# Get light profile of the lens (and compact matter if required)
####################################################################################
cmd_list = [
    molet_home+"lens_light_mass/vkl_llm/bin/llm",
    infile,
    path
]
#print(" ".join(cmd_list))
myprocess("Getting light profile of the lens...",cmd_list)



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
            "python",
            molet_home+"variability/extrinsic/match_to_gerlumph/match_to_gerlumph.py",
            molet_home+"data/gerlumph_database/gerlumph.db",
            path
        ]
        #print(" ".join(cmd_list))
        myprocess("Matching macro-images to GERLUMPH maps...",cmd_list)  
    
        cmd_list = [
            molet_home+"variability/extrinsic/gerlumph_moving_source/bin/moving_source",
            infile,
            path
        ]
        #print(" ".join(cmd_list))
        myprocess("Getting microlensing variability for each image...",cmd_list)  


cmd_list = []
cmd_list = [
    molet_home+"instrument_combined_light/bin/combine_light",
    infile,
    path
]
print(" ".join(cmd_list))
myprocess("Combining light components and including instrumental effects...",cmd_list) 
    



    



print("\nCompleted successfully!\n")
print("Output in: ",path+"output/")
