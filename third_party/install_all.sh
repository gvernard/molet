#!/bin/bash
while getopts ":g:m:s:" opt; do
  case $opt in
    g) WITH_GPU="$OPTARG"
    ;;
    m) MAP_PATH="$OPTARG"
    ;;
    s) SCRIPT_DIR="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac
done

if [ -z "$WITH_GPU" ] || [ -z "$MAP_PATH" ] || [ -z "$SCRIPT_DIR" ]
then
    echo "You need to specify all of the following command line options:"
    echo " -g: either 'yes' or 'no', enables GPU support"
    echo " -m: the absolute path to a directory with the GERLUMPH maps (can be dummy if not used)"
    echo " -s: the absolute path to the directory where all the third party libraries will be installed"
    exit 1;
fi



LIBDIR=$SCRIPT_DIR/libraries
SRCDIR=$SCRIPT_DIR/src


mkdir -p $LIBDIR
mkdir -p $SRCDIR
cd $SRCDIR



# Download, extract, and rename all the external libraries
#########################################################################################################
wget ftp://ftp.fftw.org/pub/fftw/fftw-3.3.8.tar.gz && tar -xvf fftw-3.3.8.tar.gz && mv fftw-3.3.8 fftw3
wget http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-3.47.tar.gz && tar -xvf cfitsio-3.47.tar.gz && mv cfitsio-3.47 cfitsio 
wget https://heasarc.gsfc.nasa.gov/fitsio/CCfits-2.5/CCfits-2.5.tar.gz && tar -xvf CCfits-2.5.tar.gz
wget --no-check-certificate https://sourceforge.net/projects/libpng/files/libpng16/1.6.37/libpng-1.6.37.tar.gz && tar -xvf libpng-1.6.37.tar.gz && mv libpng-1.6.37 libpng
wget https://gmplib.org/download/gmp/gmp-6.1.2.tar.xz && tar -xvf gmp-6.1.2.tar.xz && mv gmp-6.1.2 gmp
wget http://www.mpfr.org/mpfr-4.0.2/mpfr-4.0.2.tar.gz && tar -xvf mpfr-4.0.2.tar.gz && mv mpfr-4.0.2 mpfr
wget https://boostorg.jfrog.io/artifactory/main/release/1.66.0/source/boost_1_66_0.tar.gz && tar -xvf boost_1_66_0.tar.gz && mv boost_1_66_0 boost
wget https://github.com/CGAL/cgal/archive/releases/CGAL-4.11.2.tar.gz && tar -xvf CGAL-4.11.2.tar.gz && mv cgal-releases-CGAL-4.11.2 cgal
wget https://github.com/stedolan/jq/releases/download/jq-1.5/jq-1.5.tar.gz && tar -xvf jq-1.5.tar.gz && mv jq-1.5 jq
wget https://www.sqlite.org/2022/sqlite-autoconf-3380500.tar.gz && tar -xvf sqlite-autoconf-3380500.tar.gz && mv sqlite-autoconf-3380500 sqlite3
git clone https://github.com/open-source-parsers/jsoncpp.git
git clone https://github.com/gvernard/vkl_lib.git
git clone https://github.com/gvernard/gerlumphpp.git



# Install all libraries
#########################################################################################################
# Install fftw
cd fftw3 \
    &&	./configure --prefix=$LIBDIR/fftw --enable-shared \
    &&	make CFLAGS=-fPIC \
    &&	make install \
    &&  cd $SRCDIR

# Install cfitsio
cd cfitsio \
    &&	./configure --prefix=$LIBDIR/cfitsio \
    &&	make \
    &&	make install \
    &&  cd $SRCDIR
	
# Install CCfits
cd CCfits \
    &&	./configure --prefix=$LIBDIR/CCfits --with-cfitsio=$LIBDIR/cfitsio \
    &&	make \
    &&	make install \
    &&  cd $SRCDIR

# Install libpng
cd libpng \
    &&	./configure --prefix=$LIBDIR/libpng \
    &&	make check \
    &&	make install \
    &&  cd $SRCDIR

# Install gmp
cd gmp \
    &&	./configure --prefix=$LIBDIR/gmp \
    &&	make \
    &&	make install \
    &&  cd $SRCDIR

# Install mpfr
cd mpfr \
    &&	./configure --with-gmp=$LIBDIR/gmp --prefix=$LIBDIR/mpfr \
    &&	make \
    &&	make install \
    &&  cd $SRCDIR

# Install boost
cd boost \
    &&	./bootstrap.sh --with-libraries=thread \
    &&	./b2 --prefix=$LIBDIR/boost install \
    &&  cd $SRCDIR

# Install CGAL
cd cgal \
    &&	mkdir -p build_CGAL \
    &&	cd build_CGAL \
    &&  cmake \
	    -DGMP_INCLUDE_DIR=$LIBDIR/gmp/include \
	    -DGMP_LIBRARIES=$LIBDIR/gmp/lib \
	    -DMPFR_INCLUDE_DIR=$LIBDIR/mpfr/include \
	    -DMPFR_LIBRARIES=$LIBDIR/mpfr/lib \
	    -DBOOST_INCLUDEDIR=$LIBDIR/boost/include \
	    -DCMAKE_INSTALL_PREFIX=$LIBDIR/CGAL \
	    .. \
    &&	make \
    &&	make install \
    &&  cd $SRCDIR

# Install jsoncpp
cd jsoncpp \
    &&	mkdir -p build \
    &&	cd build \
    &&  cmake \
	    -DCMAKE_BUILD_TYPE=release \
	    -DBUILD_SHARED_LIBS=ON \
	    -DARCHIVE_INSTALL_DIR=. \
	    -DCMAKE_INSTALL_PREFIX=$LIBDIR/jsoncpp \
	    -G "Unix Makefiles" \
	    .. \
    &&	make \
    &&	make install \
    &&  cd $SRCDIR

# Install jq
cd jq \
    &&  autoreconf -i \
    && ./configure --disable-maintainer-mode --prefix=$LIBDIR/jq \
    &&	make \
    &&	make install \
    &&  cd $SRCDIR

# Install sqlite3
cd sqlite3 \
    && ./configure --prefix=$LIBDIR/sqlite3 \
    &&	make \
    &&  make install \
    &&  cd $SRCDIR

# Install vkl_lib
cd vkl_lib \
    &&  autoreconf -i \
    &&	./configure --prefix=$LIBDIR/vkl_lib --with-cfitsio=$LIBDIR/cfitsio --with-CCfits=$LIBDIR/CCfits --with-gmp=$LIBDIR/gmp --with-CGAL=$LIBDIR/CGAL --with-jsoncpp=$LIBDIR/jsoncpp \
    &&  make \
    &&  make install \
    &&  cd $SRCDIR

# Install gerlumphpp
cd gerlumphpp \
    &&  autoreconf -i \
    &&  ./configure --prefix=$LIBDIR/gerlumphpp --with-map-path=$MAP_PATH --with-cfitsio=$LIBDIR/cfitsio --with-CCfits=$LIBDIR/CCfits --with-png=$LIBDIR/libpng --with-fftw3=$LIBDIR/fftw --enable-gpu=$WITH_GPU \
    &&  make \
    &&  make install \
    &&  cd $SRCDIR


# Need to change directory
MY_PATH=$(dirname "$0")            # relative
MY_PATH=$(cd "$MY_PATH" && pwd)

configure="./configure --with-jq=${LIBDIR}/jq --with-fftw3=${LIBDIR}/fftw --with-cfitsio=${LIBDIR}/cfitsio --with-CCfits=${LIBDIR}/CCfits --with-gmp=$LIBDIR/gmp --with-CGAL=$LIBDIR/CGAL --with-jsoncpp=${LIBDIR}/jsoncpp --with-png=${LIBDIR}/libpng --with-sqlite3=${LIBDIR}/sqlite3 --with-vkl=${LIBDIR}/vkl_lib --with-gerlumph=${LIBDIR}/gerlumphpp"
echo $configure | cat > ${MY_PATH}/configure_command.txt


echo ""
echo $(tput setaf 2)Success! Third-party libraries have been installed$(tput sgr0) - but check above anyway for any missed errors
echo "Run 'cd ..; autoreconf -i' and then the following command to configure MOLET (also saved in $(tput smul)${MY_PATH}/configure_command.txt$(tput rmul)):"
echo ""
echo $(tput setaf 1)$configure$(tput sgr0)
echo ""
echo "Once you are sure that MOLET can compile, you can also delete the $SRCDIR directory"
