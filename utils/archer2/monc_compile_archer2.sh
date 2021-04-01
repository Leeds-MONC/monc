export PATH=$PATH:/work/y07/shared/umshared/bin
export PATH=$PATH:/work/y07/shared/umshared/software/bin
. mosrs-setup-gpg-agent

module use --append /work/y07/shared/archer2-modules/modulefiles-cse-pyvenvs
module use --append /work/y07/shared/archer2-modules/modulefiles-cse-pymods
module use --append /work/y07/shared/archer2-modules/modulefiles-cse-utils
module use --append /work/y07/shared/archer2-modules/modulefiles-cse-libs
module use --append /work/y07/shared/archer2-modules/modulefiles-cse-apps
module use --append /opt/cray/pe/perftools/20.10.0/modulefiles
module use --append /opt/cray/pe/perftools/20.09.0/modulefiles
module use --append /opt/cray/pe/craype/2.7.0/modulefiles
module use --append /usr/local/Modules/modulefiles
module use --append /opt/cray/pe/cpe-prgenv/7.0.0
module use --append /opt/cray/pe/modulefiles
module use --append /opt/cray/pe/craype-targets/default/modulefiles
module use --append /opt/modulefiles
module use --append /opt/cray/modulefiles
module load cpe-gnu
module load gcc/9.3.0
module load craype
module load craype-x86-rome
module load --notuasked libfabric
module load craype-network-ofi
module load cray-dsmml
module load perftools-base
module load xpmem
module load cray-mpich
module load cray-libsci
module load --notuasked bolt
module load --notuasked /work/y07/shared/archer2-modules/modulefiles-cse/epcc-setup-env
module load /usr/local/share/epcc-module/epcc-module-loader
module load cray-netcdf-hdf5parallel
module load cray-hdf5-parallel
module load cray-fftw/3.3.8.7
module load petsc/3.13.3
module load atp
export ATP_ENABLED=1

fcm make -j4 -f fcm-make/monc-cray-gnu.cfg
