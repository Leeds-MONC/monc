export PATH=$PATH:/work/y07/shared/umshared/bin
export PATH=$PATH:/work/y07/shared/umshared/software/bin
. mosrs-setup-gpg-agent

module restore PrgEnv-cray
module load cpe-gnu
module load gcc/9.3.0
module load cray-netcdf-hdf5parallel
module load cray-hdf5-parallel
module load cray-fftw/3.3.8.7
module load petsc/3.13.3

fcm make -j4 -f fcm-make/monc-cray-gnu-safe.cfg
