#!/bin/bash
MYLIBS=/home/giorgos/myLibraries/libraries/
g++ -std=c++17 \
    -I${MYLIBS}Catch2/include \
    -I${MYLIBS}cfitsio/include \
    -I${MYLIBS}CCfits/include \
    -I${MYLIBS}jsoncpp/include \
    -L${MYLIBS}Catch2/lib \
    -L${MYLIBS}cfitsio/lib \
    -L${MYLIBS}CCfits/lib \
    -L${MYLIBS}jsoncpp/lib \
    -Wl,-rpath,${MYLIBS}Catch2/lib \
    -Wl,-rpath,${MYLIBS}cfitsio/lib \
    -Wl,-rpath,${MYLIBS}CCfits/lib \
    -Wl,-rpath,${MYLIBS}jsoncpp/lib \
    tests.cpp -o tests -lcfitsio -lCCfits -ljsoncpp -lCatch2Main -lCatch2
