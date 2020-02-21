import re
import sys
import json

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u


f          = open(sys.argv[1],'r')
input_str  = f.read()
input_str  = re.sub(re.compile("/\*.*?\*/",re.DOTALL),"",input_str)
input_str  = re.sub(re.compile("//.*?\n" ),"",input_str)
myinput    = json.loads(input_str)


outfile = sys.argv[2]+"output/angular_diameter_distances.json"
H0      = myinput["cosmology"]["H0"]  # in km s^-1 Mpc^-1
Om0     = myinput["cosmology"]["Om0"] # Omega matter at t=t_0 (now)
cosmo   = FlatLambdaCDM(H0=H0 * u.km / u.s / u.Mpc, Om0=Om0)


Ds_arr  = []
Dl_arr  = []
Dls_arr = []

zl   = float( myinput["lenses"][-1]["redshift"] )
zs   = float( myinput["source"]["redshift"] )
Dl_arr.append( cosmo.angular_diameter_distance(zl).value )           # in Mpc
Ds_arr.append( cosmo.angular_diameter_distance(zs).value )           # in Mpc
Dls_arr.append( cosmo.angular_diameter_distance_z1z2(zl,zs).value )  # in Mpc

# The following assumes that the different lensing planes are in order of increasing redshift
if len(myinput["lenses"]) > 1:
    for i in range(1,len(myinput["lenses"])+1):
        zl = myinput["lenses"][-1-i]["redshift"]
        zs = myinput["lenses"][-1-(i-1)]["redshift"]
        Dl_arr.append( cosmo.angular_diameter_distance(zl).value )           # in Mpc
        Ds_arr.append( cosmo.angular_diameter_distance(zs).value )           # in Mpc
        Dls_arr.append( cosmo.angular_diameter_distance_z1z2(zl,zs).value )  # in Mpc

output = []
for i in range(len(Dl_arr)):
    output.append( {
        "Dl": Dl_arr[i],
        "Ds": Ds_arr[i],
        "Dls": Dls_arr[i]
    } )


with open(outfile,"w") as file:
    json.dump(output,file)
