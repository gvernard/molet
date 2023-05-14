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
mock = sys.argv[2]
band_name = sys.argv[3]

# Read main input that produced the mocks
f          = open(path + "molet_input.json",'r')
input_str  = f.read()
input_str  = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str  = re.sub(re.compile("//.*?\n" ),"",input_str)
myinput    = json.loads(input_str)
#band_name = myinput["instruments"][0]["name"]

# Read the multiple images
f          = open(path + "output/multiple_images.json",'r')
input_str  = f.read()
input_str  = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str  = re.sub(re.compile("//.*?\n" ),"",input_str)
images     = json.loads(input_str)

# Read the continuous light curves
f = open(path+mock+'/'+band_name+'_LC_continuous.json','r')
input_str  = f.read()
input_str  = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str  = re.sub(re.compile("//.*?\n" ),"",input_str)
lc_cont    = json.loads(input_str)

# Read the sampled light curves
f = open(path+mock+'/'+band_name+'_LC_sampled.json','r')
input_str  = f.read()
input_str  = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str  = re.sub(re.compile("//.*?\n" ),"",input_str)
lc_samp    = json.loads(input_str)






# Image indices, omitting the maximum image that is demagnified
td_max = 0
td_max_index = 0
for i in range(0,len(images)):    
    if float(images[i]["dt"]) > td_max:
        td_max = float(images[i]["dt"])
        td_max_index = i
img_indices = []
for i in range(0,len(images)):
    #img_indices.append(i)
    if i != td_max_index:
        img_indices.append(i)


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
#ax.set_xlim(61150,61340)
#ax.set_ylim(0.5,-6)

for i in range(0,len(img_indices)):
    index = img_indices[i]
    ax.plot(lc_cont[index]["time"],lc_cont[index]["signal"],color=faded_colors[i],zorder=1)
    ax.scatter(lc_samp[index]["time"],lc_samp[index]["signal"],color=colors[i],zorder=2)

    max_index = lc_cont[index]["signal"].index(np.amin(lc_cont[index]["signal"]))
    ax.text(lc_cont[index]["time"][max_index],lc_cont[index]["signal"][max_index],index,horizontalalignment='right',verticalalignment='bottom',color=colors[i],fontsize=22,zorder=3)

curves_ymin = ax.get_ylim()[0]
curves_ymax = ax.get_ylim()[-1]


plt.savefig('light_curves.png')

 
