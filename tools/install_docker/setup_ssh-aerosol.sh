set -e
# Install Blitz
wget https://github.com/blitzpp/blitz/archive/refs/heads/master.zip
unzip master.zip
cd blitz-main
mkdir build
cd build
cmake ..
make lib
make install
cd /root
# Install SWIG
wget http://prdownloads.sourceforge.net/swig/swig-4.0.2.tar.gz
tar xvf swig-4.0.2.tar.gz
cd swig-4.0.2
./configure --without-pcre --prefix=/usr/local
make
make install
cd /root
# Install NetCDF C
wget https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.7.3.tar.gz
tar xvf v4.7.3.tar.gz
cd netcdf-c-4.7.3
./configure LDFLAGS="-L/usr/local/lib" --disable-hdf5 --prefix=/usr/local
make
make install
cd /root
# Set environment
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
# Path to blitz and netcdf files
export CPLUS_INCLUDE_PATH=/usr/local/include:$CPLUS_INCLUDE_PATH
# path to netcdf.mod
export F90FLAGS=-I/usr/local/include:$F90FLAGS
# path to swig
export SWIG=/usr/local/bin/swig
# Clean up
rm -rf blitz-main swig-4.0.2 netcdf-fortran-4.5.2 netcdf-c-4.7.3 master.zip v4.5.2.tar.gz v4.7.3.tar.gz swig-4.0.2.tar.gz
