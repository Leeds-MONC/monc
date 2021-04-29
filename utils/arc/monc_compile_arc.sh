#! /bin/bash

#Script to compile Monc on ARC4:

module purge
module load user
module switch intel gnu
module switch openmpi mvapich2
module load fftw netcdf hdf5 fcm

fcm make -j4 -f fcm-make/monc-arc4-gnu.cfg
