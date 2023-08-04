import re
import sys
import json
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from matplotlib.gridspec import GridSpec
from astropy.io import fits
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import ConnectionPatch




path = sys.argv[1]
band_name = sys.argv[2]
index = int(sys.argv[3])


# Read the intrinsic light curves
f = open(path+'input_files/'+band_name+'_LC_intrinsic.json','r')
input_str  = f.read()
input_str  = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str  = re.sub(re.compile("//.*?\n" ),"",input_str)
lc_in    = json.loads(input_str)

# Read the unmicrolensed light curves
f = open(path+'output/'+band_name+'_LC_unmicro.json','r')
input_str  = f.read()
input_str  = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str  = re.sub(re.compile("//.*?\n" ),"",input_str)
lc_un    = json.loads(input_str)







prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
    
def make_rgb_transparent(rgb,bg_rgb,alpha):
    return [(alpha * c1 + (1 - alpha) * c2)/255.0 for (c1, c2) in zip(rgb, bg_rgb)]
faded_colors = []
for i in range(0,len(colors)):
    hex = colors[i].lstrip('#')
    rgb = tuple(int(hex[k:k+2], 16) for k in (0, 2, 4))
    new_rgb = make_rgb_transparent(rgb,[255.0,255.0,255.0],0.6)
    faded_colors.append(tuple(new_rgb))




    



    
#fig = plt.figure(figsize=(20,8))
fig,ax = plt.subplots(1,figsize=(18,7))

ax.invert_yaxis()
ax.set_xlabel("t [days]")
ax.set_ylabel("Mag")
#ax.set_xlim(1230,1350)
#ax.set_ylim(15.9,12.5)

ax.plot(lc_in[index]["time"],lc_in[index]["signal"],color=colors[0],zorder=2)
ax.scatter(lc_un[index]["time"],lc_un[index]["signal"],color=colors[1],zorder=1,s=5)

plt.savefig('in_un_light_curves.png')

 
