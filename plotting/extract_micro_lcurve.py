import json
import sys
import re


# Read file with the extrinsic (only microlensing) light curves
lc_path = sys.argv[1]
f = open(lc_path,'r')
input_str  = f.read()
input_str  = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str  = re.sub(re.compile("//.*?\n" ),"",input_str)
ex_lc      = json.loads(input_str)


# Get light curve index
index_lc_ex = int(sys.argv[2])

images = []
for i in range(0,len(ex_lc)):
    image = []
    lc = {}
    lc["time"]   = ex_lc[i][index_lc_ex]["time"]
    lc["signal"] = ex_lc[i][index_lc_ex]["signal"]
    #print(lc["signal"][0])
    image.append(lc)
    images.append(image)

with open('extrinsic_lc.json','w') as outfile:
    json.dump(images,outfile)
