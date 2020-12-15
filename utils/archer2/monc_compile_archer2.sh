wget https://github.com/metomi/fcm/archive/2019.09.0.zip
unzip 2019.09.0.zip 
export PATH=$HOME/fcm-2019.09.0/bin:$PATH

git clone https://github.com/Leeds-MONC/monc
cd monc/

module load cray-netcdf
module load cray-fftw
module load cray-hdf5
# PETSC is not yet available, but will be made available according to
# instructors on "Introduction to ARCHER2 for software developers" course
# line below will make sure build can continue without PETSC
export PETSC_DIR=""

fcm make -f fcm-make/monc-cray-cray.cfg 
