#!/bin/bash

#source /home/es732/miniconda3/etc/profile.d/conda.sh
#conda activate mbx-gcc5

export MBX_HOME=$PWD # $HOME/codes/MBX

# export MODULEPATH=/projects/builder-group/jpg/modulefiles/applications:$MODULEPATH
# module load cmake/3.13.4

# module load intel/2018.1.163 gsl openmpi_ib fftw
# module load gnu


# Compile the driver
cd plugins/i-pi/src/main
cp Makefile_csd3 Makefile
make
cd -


