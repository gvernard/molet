#!/bin/bash
MYLIBS=/home/giorgos/myLibraries/
g++ -std=c++17 \
    -I${MYLIBS}Catch2/ \
    -I${MYLIBS}cfitsio/include \
    -I${MYLIBS}CCfits/include \
    -I${MYLIBS}jsoncpp/include \
    -L${MYLIBS}cfitsio/lib \
    -L${MYLIBS}CCfits/lib \
    -L${MYLIBS}jsoncpp/lib64 \
    -Wl,-rpath,${MYLIBS}cfitsio/lib \
    -Wl,-rpath,${MYLIBS}CCfits/lib \
    tests.cpp -o tests -lcfitsio -lCCfits -ljsoncpp
