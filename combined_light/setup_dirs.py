import os
import re
import sys
import json
import shutil

infile     = sys.argv[1]
in_path    = sys.argv[2]
out_path   = sys.argv[3]

f          = open(infile,'r')
input_str  = f.read()
input_str  = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str  = re.sub(re.compile("//.*?\n" ),"",input_str)
myinput    = json.loads(input_str)

# Choose the first band (number of light curves, in or ex, should be the same in all bands
band_name = myinput["instrument"]["bands"][0]["name"]

lc_path = ""
if myinput["point_source"]["variability"]["intrinsic"]["type"] == "custom":
    lc_path = in_path + "input_files/intrinsic_light_curves.json" 
else:
    lc_path = out_path + "output/intrinsic_light_curves.json"
f = open(lc_path,'r')
input_str  = f.read()
input_str  = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str  = re.sub(re.compile("//.*?\n" ),"",input_str)
in_lc      = json.loads(input_str)
Nin = len(in_lc[band_name])



if myinput["point_source"]["variability"]["extrinsic"]["type"] == "custom":
    lc_path = in_path + "input_files/extrinsic_light_curves.json" 
else:
    lc_path = out_path + "output/extrinsic_light_curves.json" 
f = open(lc_path,'r')
input_str  = f.read()
input_str  = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str  = re.sub(re.compile("//.*?\n" ),"",input_str)
ex_lc      = json.loads(input_str)

for q in range(0,len(ex_lc)):
    if len(ex_lc[q][band_name]) > 0:
        Nex = len(ex_lc[q][band_name])
        break


for i in range(0,Nin):
    for j in range(0,Nex):
        mock = "mock_%04d_%04d" % (i,j)
        if os.path.isdir(out_path+mock):
            shutil.rmtree(out_path+mock)
        os.mkdir(out_path+mock)

        

