import os
import re
import sys
import json
import math
import sqlite3
import urllib.parse

'''
Reuiqres:
- multiple_images.json
'''

dbfile     = sys.argv[1]
map_path   = sys.argv[2]
out_path   = sys.argv[3]

f          = open(out_path+"output/multiple_images.json",'r')
input_str  = f.read()
input_str  = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str  = re.sub(re.compile("//.*?\n" ),"",input_str)
images     = json.loads(input_str)

outfile = out_path+"output/gerlumph_maps.json"



conn = sqlite3.connect(dbfile)
cur = conn.cursor()

gerlumph_maps = []
kg_sep = 0.05
for image in images:
    
    k0 = image["k"]
    g0 = image["g"]
    s0 = image["s"]
    cur.execute("SELECT id,k,g,s,(k-?)*(k-?)+(g-?)*(g-?) AS d,ABS(s-?) as ds FROM gerlumph GROUP BY d,ds HAVING d < ?*? ORDER BY d ASC,ds ASC LIMIT 10;",[k0,k0,g0,g0,s0,kg_sep,kg_sep])
    cases = cur.fetchall()

    mymap = {}
    if len(cases) == 0:
        mymap["id"]  = "none"
        mymap["k"]   = 0
        mymap["g"]   = 0
        mymap["s"]   = 0
        mymap["dkg"] = 0
        mymap["ds"]  = 0
    else: 
        mymap["id"]  = cases[0][0]
        mymap["k"]   = float(cases[0][1])
        mymap["g"]   = float(cases[0][2])
        mymap["s"]   = float(cases[0][3])
        mymap["dkg"] = math.sqrt(cases[0][4])
        mymap["ds"]  = cases[0][5]
    gerlumph_maps.append(mymap)
    
conn.close()


with open(outfile,"w") as file:
    json.dump(gerlumph_maps,file,indent=4)

missing = []
for mymap in gerlumph_maps:
    if mymap["id"] != "none":
        if not os.path.isdir(map_path+str(mymap["id"])):
            missing.append(str(mymap["id"]))
            
if len(missing) > 0:
    missing = list(dict.fromkeys(missing))
    missing = list( map(lambda x: "ids[]="+x,missing) )
    print("ATTENTION: missing GERLUMPH maps from '"+map_path+"' !",file=sys.stderr)
    print("Click the following link to download them:",file=sys.stderr)
    print("\n",file=sys.stderr)
    myurl = "gerlumph.swin.edu.au/inc/generic/put_to_cart.php?"+"&".join(missing)
    print("      http://"+urllib.parse.quote(myurl),file=sys.stderr)
    print("\n",file=sys.stderr)

