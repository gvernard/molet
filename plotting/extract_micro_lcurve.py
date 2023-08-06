import json
import sys
import re
import numpy as np
import matplotlib.pyplot as plt


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
    if len(ex_lc[i]) != 0:
        lc = {}
        lc["time"]    = ex_lc[i][index_lc_ex]["time"]
        lc["signal"]  = ex_lc[i][index_lc_ex]["signal"]
        lc["dsignal"] = ex_lc[i][index_lc_ex]["dsignal"]
        image.append(lc)
        #print(lc["signal"][0])
    images.append(image)

with open('extrinsic_lc.json','w') as outfile:
    json.dump(images,outfile)





    
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

#ax.invert_yaxis()
ax.set_xlabel("t [days]")
ax.set_ylabel(r"$\mu$")
#ax.set_xlim(61150,61340)
#ax.set_ylim(0.5,-6)

for q in range(0,len(images)):
    if len(images[q]) != 0:
        time    = np.array(images[q][0]["time"])
        signal  = np.array(images[q][0]["signal"])
        dsignal = np.array(images[q][0]["dsignal"])
        
        ax.plot(time,signal,color=colors[q],zorder=1)
        ax.fill_between(time,signal-dsignal,signal+dsignal,color=faded_colors[q],zorder=1,alpha=0.5)
        #ax.scatter(time,lc_samp[index]["signal"],color=colors[i],zorder=2)
        
        max_index = list(signal).index(np.amax(signal))
        ax.text(time[max_index],signal[max_index],q,horizontalalignment='right',verticalalignment='bottom',color=colors[q],fontsize=22,zorder=3)

curves_ymin = ax.get_ylim()[0]
curves_ymax = ax.get_ylim()[-1]


plt.savefig('micro_light_curves.png')
