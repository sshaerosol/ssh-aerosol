#!/bin/bash
set -e

# Save the current directory
WORKDIR=$(pwd)

# Update and install dependencies
sudo apt-get -y update && sudo apt-get install -y \
    python3 \
    python3-pip \
    cmake \
    vim \
    wget \
    unzip \
    curl \
    build-essential \
    zlib1g-dev \
    libcurl4-openssl-dev \
    m4 \
    gfortran \
    scons \
    libnetcdff-dev

# Blitz installation
echo "Installing Blitz..."
cd "$WORKDIR"
wget https://github.com/blitzpp/blitz/archive/refs/heads/master.zip
unzip master.zip
cd blitz-main
mkdir build
cd build
cmake ..
make lib
sudo make install
cd "$WORKDIR"

# SWIG installation
echo "Installing SWIG..."
wget http://prdownloads.sourceforge.net/swig/swig-4.0.2.tar.gz
tar xvf swig-4.0.2.tar.gz
cd swig-4.0.2
./configure --without-pcre --prefix=/usr/local
make
sudo make install
cd "$WORKDIR"

# NetCDF C installation
echo "Installing NetCDF C..."
wget https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.7.3.tar.gz
tar xvf v4.7.3.tar.gz
cd netcdf-c-4.7.3
./configure LDFLAGS="-L/usr/local/lib" --disable-hdf5 --prefix=/usr/local
make
sudo make install
cd "$WORKDIR"

# Set environment variables
echo "Setting environment variables..."
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
export CPLUS_INCLUDE_PATH=/usr/local/include:$CPLUS_INCLUDE_PATH
export F90FLAGS=-I/usr/local/include:$F90FLAGS
export SWIG=/usr/local/bin/swig

# Add environment variables
echo "Adding environment variables to ~/.bashrc..."
echo "export LD_LIBRARY_PATH=/usr/local/lib:\$LD_LIBRARY_PATH" >> ~/.bashrc
echo "export CPLUS_INCLUDE_PATH=/usr/local/include:\$CPLUS_INCLUDE_PATH" >> ~/.bashrc
echo "export F90FLAGS=-I/usr/local/include:\$F90FLAGS" >> ~/.bashrc
echo "export SWIG=/usr/local/bin/swig" >> ~/.bashrc

# Clean up installation files
echo "Cleaning up..."
rm -rf blitz-main swig-4.0.2 netcdf-c-4.7.3 master.zip swig-4.0.2.tar.gz v4.7.3.tar.gz

echo "Installation completed."
