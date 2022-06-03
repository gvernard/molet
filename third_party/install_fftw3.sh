wget ftp://ftp.fftw.org/pub/fftw/fftw-3.3.8.tar.gz
tar -xvf fftw-3.3.8.tar.gz
cd fftw-3.3.8
./configure --prefix=$LIBDIR/fftw --enable-shared
make CFLAGS=-fPIC
make install
