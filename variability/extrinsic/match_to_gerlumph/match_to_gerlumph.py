import re
import sys
import json
import math
import sqlite3

'''
Reuiqres:
- multiple_images.json
'''

dbfile     = sys.argv[1]
path       = sys.argv[2]

f          = open(path+"output/multiple_images.json",'r')
input_str  = f.read()
input_str  = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str  = re.sub(re.compile("//.*?\n" ),"",input_str)
images     = json.loads(input_str)

outfile = path+"output/gerlumph_maps.json"



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

    map = {}
    if len(cases) == 0:
        map["id"]  = "none"
        map["k"]   = 0
        map["g"]   = 0
        map["s"]   = 0
        map["dkg"] = 0
        map["ds"]  = 0
    else: 
        map["id"]  = cases[0][0]
        map["k"]   = float(cases[0][1])
        map["g"]   = float(cases[0][2])
        map["s"]   = float(cases[0][3])
        map["dkg"] = math.sqrt(cases[0][4])
        map["ds"]  = cases[0][5]
    gerlumph_maps.append(map)
    
conn.close()



with open(outfile,"w") as file:
    json.dump(gerlumph_maps,file,indent=4)
