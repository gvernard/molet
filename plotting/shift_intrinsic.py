import re
import sys
import os
import json
import numpy as np


in_json = sys.argv[1]
f          = open(in_json,'r')
input_str  = f.read()
input_str  = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str  = re.sub(re.compile("//.*?\n" ),"",input_str)
myinput    = json.loads(input_str)

print(len(myinput[0]["time"]))
print(len(myinput[0]["signal"]))

myinput[0]["signal"] = list(np.add(myinput[0]["signal"],9))

print(np.amin(myinput[0]["signal"]))
print(np.mean(myinput[0]["signal"]))
print(np.amax(myinput[0]["signal"]))


with open('json_data.json', 'w') as outfile:
    json.dump(myinput,outfile)
