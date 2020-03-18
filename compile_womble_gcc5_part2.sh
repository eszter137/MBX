#!/bin/bash

#source /home/es732/miniconda3/etc/profile.d/conda.sh
#conda activate mbx-gcc5

export MBX_HOME=$PWD # $HOME/codes/MBX

# export MODULEPATH=/projects/builder-group/jpg/modulefiles/applications:$MODULEPATH
# module load cmake/3.13.4

# module load intel/2018.1.163 gsl openmpi_ib fftw
# module load gnu

g++-5 --version
#rm -rf build install
#cmake -DCMAKE_BUILD_TYPE=Debug -DUSE_OPENMP=True -DCMAKE_CXX_FLAGS=" -fPIC -O2 -Wall" -DCMAKE_CXX_COMPILER=g++-5 -DCMAKE_C_COMPILER=gcc-5 -H. -Bbuild
#cd build
#make -j 8 CXX=g++-5 CC=gcc-5
#make install
#cd ../

# Compile the driver
cd plugins/i-pi/src/main
cp Makefile_womble Makefile
make
cd -


