# Mock Lenses in Time

This is a collection of software creating simulated galaxy-galaxy lenses as a function of time.




## Installing

There are two different ways to setup a working installation of MOLET.
The quickest is by using the docker provided images, however, the docker virtualization makes the computations signigicantly slower.
This approach is recommended for getting to know the code, test runs, even some development.
Setting up a local installation is a bit more time consuming, but the code can run considerably faster in this way.
The latter is the recommended approach for using MOLET in production mode, e.g. when one needs to simulate more than a handful of systems.
Furthermore, enabling GPU support dramatically increase the speed of microlensing related computations.

Here are some timing results when running two examples found in the 'tests' directory on a generic laptop:  

Test 2 (including microlensing):
- docker CPU: 120s
- local CPU: 80s (x1.5 faster)
- local GPU: 22s (x5.5 faster)

Test 5 (no microlensing, producing just a cutout):
- docker CPU: 8.8s
- local CPU: 1.8s (x5 faster)



### Docker container
[Docker](https://docs.docker.com/get-started/overview/) images provide a pre-compilled environment to run MOLET with all required third-party software and libraries installed.
These include [gerlumphpp](https://github.com/gvernard/gerlumphpp) and [vkl_lib](https://github.com/gvernard/vkl_lib).
After installing docker through a package manager, type:

```
docker pull gvernard/molet:production
docker run -it --rm -v </path/to/gerlumph/maps/>:/home/gerlumph_maps:ro gvernard/molet:production
```

which pulls the production docker image of MOLET and runs a container.
If you want to do some development, i.e. being able to re-compile the code, then use the gvernard/molet:development docker image.

Notice that in the command above the host directory (i.e. outside the docker virtual environment) */path/to/gerlumph/maps/* containing the GERLUMPH maps is mounted in the docker environment as an external volume.
Because whatever is written inside a docker container is deleted when the container exits, it is recommended to mount an additional volume to hold all the user input and output files (i.e. making a results-directory in your system visible to docker).
This can be achieved by modifying the above command as follows:

```
docker run -it --rm -v </path/to/gerlumph/maps/>:/home/gerlumph_maps:ro -v </path/to/molet/output/>:/home/molet_output/:rw gvernard/molet:production
```

You can find more information on this standard docker usage [here](https://docs.docker.com/storage/volumes/).


### Installing manually

#### Prerequisites
The following packages/programs are required in order to proceed with the installation of MOLET: **CMake**, **Autotools** (for the 'autoreconf' command), an existing **C/C++** compiler compatible with the c++17 standard (e.g. GCC version 8 and above should be fine), and **CUDA** (i.e. the nvcc compiler) for GPU support.
These programs need to be available in the command line but are standard packages that can be easily installed via package managers.

#### Third-party libraries
MOLET uses third party libraries that must be pre-installed in the system in order to run it.
For convenience, the list of commands that downloads and installs all these libraries is provided in the *third_party/install_all.sh* script.

To use the install_all.sh script, the user needs to provide three command line arguments:
 - -g: either 'yes' or 'no', enables GPU support
 - -m: the absolute path to a directory with the GERLUMPH maps (can be dummy if not used)
 - -s: the absolute path to the directory where all the third party libraries will be installed
 
The first two arguments are relevant only for the microlensing part of the code, but they have to be specified anyway (one can set '-g no' and provide any directory as a path to maps if microlensing is not relevant for their simulations).
Also, these two arguments can be changed at any time, requiring only the re-compillation of the gerlumph++ library (see its [README](https://github.com/gvernard/gerlumphpp) on how to achieve this).
The third argument should be a directory anywhere in the system where the third lbraries will be installed.
MOLET will eventually be linked to the libraries installed by this script, even if some or all of them (unlikely) already exist in the system.  
Pointing to different libraries is a somewhat advanced task, which the users are advised not to undertake unless they are familiar with installing software through autotools (the './configure, make, make install' way).
*Note*: If the installation produces errors, then you can edit the script, e.g. comment out the part that downloads the source files for the libraries, and re-run it to debug.

#### Installing MOLET

Once the script has finished its dull but necessary task, it will have produced a './configure' command with the right options that link the installed libraries to MOLET.
Here is an example of what this will look like:

```
./configure --with-jq=/path/to/libraries/jq --with-fftw3=/path/to/libraries/fftw --with-cfitsio=/path/to/libraries/cfitsio --with-CCfits=/path/to/libraries/CCfits --with-gmp=/path/to/libraries/gmp --with-CGAL=/path/to/libraries/CGAL --with-jsconcpp=/path/to/libraries/jsoncpp --with-png=/path/to/libraries/libpng --with-sqlite3=/path/to/libraries/sqlite3 --with-vkl=/path/to/libraries/vkl_lib --with-gerlumph=/path/to/libraries/gerlumphpp
```

where '/path/to/libraries/' will corredpond to the directory the user passed as argument to the install_all.sh script. 
Switching to the root installation directory of MOLET and running 'autoreconf -i' followed by this command and the usual 'make, make install' will complete the installation.
MOLET can be run now, and the code can even be modified by the user and recompiled using just 'make, make install' at the root directory.
The './configure' command needs to be run only once during installation to link the libraries - only if the location of the libraries changes it will have to be modified accordingly (i.e. with the new paths provided as options) and rerun.




## Running
The main script calling the various programms and managing input and output files is *molet_driver.sh*.
Each simulation is contained in a single directory that needs to include a *.json* input file, an *input_files* directory and the eventual output.


### Input
A *.json* file is required with all the relevant parameters.
The *input_files* directory needs to include **at least** the intrinsic light curves, which are stored as an associative array with two keys, "time" in days and "signal" in aparent magnitudes, each one being a coma-separated list of values.
Similarly for an unmicrolensed intrinsic light curve component, if present.
If custom microlensing light curves are provided, another *.json* file needs to be provided in *input_files* with the same structure as the intrinsic and unmicrolensed light curve files, but having a list of light curves per image (run the tests to see examples, see below).


### Output
The output consists of an *output* directory containing separate images of the static image components and other quantities of interest.
One or more *mock_<index_in>_<index_ex>* directories are also produced, which contain the results for each realization using the provided intrinsic and extrinsic light curves.
These directories are named after the corresponding indices in the *.json* file input lists.
The purely microlensing light curves will be located in the *output* directory.
Inside each *mock_<index_in>_<index_ex>* directory there is a file with the resulting continuous (daily cadence) light curves and another one with those sampled according to the provided time vector.
These output files are in the same format as the input light curve *.json* files.
Finally, any requested image cutouts will also be located these directories.


### Running the tests

The root installation contains a directory named 'tests', where a few standard use cases are given as examples.
Inside the tests directory there is a [README](tests/README.txt) file describing the various tests.
Each test can be run by typing:

```
./molet_driver.sh tests/test_<X>/molet_input.json
```

The user can move and rename these directories anywhere in their system and use them to build their own simulations (using a mounted volume is encouraged if MOLET is run within docker, see above). 
The molet_driver.sh program will then need to be invoked with the full path to the *.json* input file.


### Plot observed light curves
A [visualization program](plotting) that uses python to plot the observed light curves is provided, which can be run (once all the required python packages are installed) as:

```
cd plotting/
python plot_light_curves.py /full/path/to/simulation/dir/ mock_<index_in>_<index_ex>
```

to produce a plot similar to the following one:

![Alt text](plotting/light_curves.png?raw=true "Example observed light curves")



## Adding new instruments

A few basic examples of instruments (i.e. a psf, some resolution, bandwidth, etc) are provided with MOLET, but the user can easily add their own.
From the root installation directory, the user needs to run:
```
./bin/aux_create_instrument /path/to/instrument/specs.json /path/to/instrument/psf.fits
```
passing the path to a .json and .fits files that contain the information on the instrument.
Examples of such files can be found in data/instrument_data/.
The .fits file is just an image of the PSF, but its dimensions (physical and in pixels) need to be provided in the specs.json file.  



## Note on using magnification maps

### Using GERLUMPH maps

MOLET uses [gerlumphpp](https://github.com/gvernard/gerlumphpp) that needs to be installed with an option specifying a fixed path to a directory containing GERLUMPH maps: **MAP_PATH**=*/path/to/gerlumph/maps/*
If the required maps are not present, MOLET automatically provides a download link.
The user is notified via email to download the compressed collection of map files in a tarball named *output.tar.gz*.
One needs to simply extract the contents of the tarball in */path/to/gerlumph/maps/*.
The maps are then organized in sub-directories named after the map ID in the gerlumph database, each containing a **map.bin** and **mapmeta.dat** file and ready to be used by MOLET.
There is a plan to replace this manual download by an automated command line download in the future. 

### Using custom maps

A list of the ID,&kappa;,&gamma;, and s of the all GERLUMPH maps (last updated on 16/06/2020) is stored in the *data/gerlumph.db* file.
This file is parsed to match the macromodel computed &kappa;,&gamma;,s to a GERLUMPH map ID, and then use this ID to find the actual map data in */path/to/gerlumph/maps/*.
One can use custom maps as long as they update the *data/gerlumph.db* and use the same format as GERLUMPH.
