import os
import re
import sys
import json
import fnmatch
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from astropy.io import fits




def updatefig(img_path):
    image = fits.getdata(img_path,ext=0)    
    image = np.flipud(np.array(image))
    im.set_array(image)
    for i in range(0,len(lab)):
        lab[i].set_text(chr(65+i))
    return tuple([im]) + tuple(lab)
#    return tuple([im])

f          = open(sys.argv[1],'r')
input_str  = f.read()
input_str  = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str  = re.sub(re.compile("//.*?\n" ),"",input_str)
myinput    = json.loads(input_str)

f          = open(sys.argv[2],'r')
input_str  = f.read()
input_str  = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str  = re.sub(re.compile("//.*?\n" ),"",input_str)
images     = json.loads(input_str)

output  = sys.argv[3]

mycmap = matplotlib.cm.get_cmap('jet')




for band in myinput["instrument"]["bands"]:    
    width  = band["field-of-view_x"]
    height = band["field-of-view_y"]
    img_xmin = -width/2.0
    img_xmax =  width/2.0
    img_ymin = -height/2.0
    img_ymax =  height/2.0

    targets = []
    for file in os.listdir(output):
        if fnmatch.fnmatch(file,"OBS_"+band["name"]+"_*.fits"):
            targets.append(file)
    sorted = sorted(targets)
    for i in range(0,len(sorted)):
        sorted[i] = output + sorted[i]
    #print(sorted)
    
    maxs = []
    for img_path in sorted:
        image = fits.getdata(img_path,ext=0)
        maxs.append(np.amax(image))
    limit_up = np.amax(maxs)
    limit_do = 0
    
    fig = plt.figure()
    ax = fig.add_subplot(111,aspect='equal',xlim=(img_xmin,img_xmax),ylim=(img_ymin,img_ymax))
    ax.set_xlabel('arcsec')
    ax.set_ylabel('arcsec')
        
    dum   = fits.getdata(sorted[0],ext=0)
    image = np.full_like(dum,0)
    im    = ax.imshow(image,interpolation='none',cmap=mycmap,extent=[img_xmin,img_xmax,img_ymin,img_ymax],vmin=limit_do,vmax=limit_up)

    # Find the maximum image
    td_max = 0
    td_max_index = 0
    for i in range(0,len(images)):    
        if float(images[i]["dt"]) > td_max:
            td_max = float(images[i]["dt"])
            td_max_index = i

    lab = []
    for i in range(0,len(images)):    
        if i != td_max_index:
            lab.append( ax.text(float(images[i]["x"])+0.25,float(images[i]["y"])+0.25,str(i),color='white',fontweight='bold') )

    ani = animation.FuncAnimation(fig,updatefig,sorted[0:],interval=10,repeat_delay=500,blit=True)
#    ani.save('basic_animation.mp4',fps=30)
    
    plt.show()
