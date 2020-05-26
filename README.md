# Mock Lenses in Time

This is a collection of software creating simulated galaxy-galaxy lenses as a function of time.



## Using GERLUMPH maps

MOLET uses [gerlumphpp](https://github.com/gvernard/gerlumphpp) that needs to be installed with an option specifying a fixed path to a directory containing GERLUMPH maps: **MAP_PATH**=*/path/to/gerlumph/maps/*
If the required maps are not present, MOLET automatically provides a download link.
The user is notified via email to download the compressed collection of map files in a tarball named *output.tar.gz*.
One needs to simply extract the contents of the tarball in */path/to/gerlumph/maps/*.
The maps are then organized in sub-directories named after the map ID in the gerlumph database, each containing a **map.bin** and **mapmeta.dat** file and ready to be used by MOLET.


## Install

### Docker container
Docker images provide an environment to run MOLET with all required third-party software and libraries installed, including [gerlumphpp](https://github.com/gvernard/gerlumphpp) and [vkl_lib](https://github.com/gvernard/vkl_lib).
After installing docker, type:

```
docker pull gvernard/molet:production
docker run -it --rm -v </path/to/gerlumph/maps/>:/home/gerlumph_maps:ro gvernard/molet:production
```

which pulls the production docker image of MOLET and runs a container.

Notice that the host directory */path/to/gerlumph/maps/* containing the GERLUMPH maps is mounted as an external volume.

An additional volume can be mounted to hold all the user input and output files (whatever is written inside a docker container is deleted when the container exits):

```
docker run -it --rm -v </path/to/gerlumph/maps/>:/home/gerlumph_maps:ro -v </path/to/molet/output/>:/home/molet_output/:rw gvernard/molet:production
```


### Installing manually
Check the [dockerfiles](dockerfiles) to see how the docker images are created and get a list of all third-party software.
One will need to install manually all the pre-requisites and set the environment variables.




## Run

The main script calling the various programms and managing input and output files is *molet_driver.sh*.
A single system simulation can in fact contain many different intrinsic and extrinsic variability light curve realizations.
Each simulation is contained in a single directory that needs to include a *.json* input file, an *input_files* directory and the eventual output.


### Input
A *.json* file is required with all the relevant parameters.
The *input_files* directory needs to include **at least** the intrinsic light curves, which are stored as an associative array with two keys, "time" in days and "signal" in aparent magnitudes, each one being a coma-separated list of values.
Similarly for an unmicrolensed intrinsic light curve component.
If custom microlensing light curves are provided, another *.json* file needs to be provided in *input_files* with the same structure as the intrinsic and unmicrolensed light curve files, but having a list of light curves per image.


### Output
The output consists of an *output* directory containing separate images of the static image components and other quantities of interest, and one or more *mock_<index_in>_<index_ex>* directories containing the results for each realization using the provided intrinsic and extrinsic light curves, named after the corresponding indices in the *.json* file input lists.
If MOLET is asked to extract microlensing light curves then these will be located in the *output* directory.
Inside each *mock_<index_in>_<index_ex>* realization directory there is a file with the final (observed) continuous (daily cadence) light curves and another one with those sampled according to the provided time vector.
These output files are in the same format as the input light curve *.json* files.
Finally, if image cutouts are requested they will be located there along with the light curve files.

### Run the tests

Inside the tests directory there is a [README](tests/README.txt) file describing the various tests.
Each test can be run by typing:

```
./molet_driver.sh tests/test_<X>/molet_input.json
```

or by copying the simulation directory in any other path and calling molet_driver.sh with the full path to the *.json* input file.


### Plot observed light curves
A [visualization program](plotting) that uses python to plot the observed light curves is provided, which can be run (once all the required python packages are installed) as:

```
cd plotting/
python plot_light_curves.py /full/path/to/simulation/dir/ mock_<index_in>_<index_ex>
```

to produce a plot similar to the following one:

![Alt text](plottin/example_observed_lc.png?raw=true "Example observesd light curves")