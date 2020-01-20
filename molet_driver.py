import os
import sys
from subprocess import Popen,PIPE


def myprocess(msg,cmd_list):
    print("%-60s" % msg,end="",flush=True)
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







# Step 1:
# Get angular diameter distances
cmd_list = [
    "python",
    molet_home+"cosmology/angular_diameter_distances.py",
    infile,
    path
]
myprocess("Getting angular diameter distances...",cmd_list)


# Step 2:
# Get extended lensed images of the source
cmd_list = [
    molet_home+"lensed_extended_source/fproject/bin/fproject",
    infile,
    path+"angular_diameter_distances.json",
    path
]
myprocess("Getting extended source lensed features...",cmd_list)



