# Mock Lenses in Time

This is a collection of software creating simulated galaxy-galaxy lenses as a function of time.



## Accessing GERLUMPH maps

MOLET uses [gerlumphpp](https://github.com/gvernard/gerlumphpp) that needs to be installed with an option specifying a fixed path to a directory containing GERLUMPH maps: **MAP_PATH**=*/path/to/gerlumph/maps/*
If the required maps are not present, MOLET automatically provides a download link.
The user is notified via email to download the compressed collection of map files in a tarball named *output.tar.gz*.
One needs to simply extract the contents of the tarball in */path/to/gerlumph/maps/*.
The maps are then organized in sub-directories named after the map ID in the gerlumph database, each containing a **map.bin** and **mapmeta.dat** file and ready to be used by MOLET.


## Install

### Running using Docker containers

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
Check [dockerfiles](dockerfiles) to see how the docker images are created and get a list of all third-party software.
One will need to install manually all the pre-requisites and set the environment variables.


## Run

### Input

### Output

### Run the tests

