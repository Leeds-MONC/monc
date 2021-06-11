#!/usr/bin/env bash

keywrdfile=$( sed "s|= |= $( pwd -P )\/fcm-make\/|g" fcm-make/keyword.cfg )
if [ ! -f ~/.metomi/fcm/keyword.cfg ]; then
  mkdir -p ~/.metomi/fcm
  echo ${keywrdfile} > ~/.metomi/fcm/keyword.cfg
else
  if cat ~/.metomi/fcm/keyword.cfg | grep -q -v "${keywrdfile}"; then
    echo ${keywrdfile} >> ~/.metomi/fcm/keyword.cfg
  fi
fi

compiler=gnu
#compiler=cray

export PATH=$PATH:/work/y07/shared/umshared/bin
export PATH=$PATH:/work/y07/shared/umshared/software/bin
. mosrs-setup-gpg-agent

module purge
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
if [ $compiler == "gnu" ]; then
  module load cpe-gnu
  module load gcc/9.3.0
elif [ $compiler == "cray" ]; then
  module load cpe-cray
  module load cce
fi
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

if [ $compiler == "gnu" ]; then
  fcmconfig="fcm-make/monc-cray-gnu.cfg"
elif [ $compiler == "cray" ]; then
  fcmconfig="fcm-make/monc-cray-cray.cfg"
fi

echo "Compile options: "
echo "(1)   MONC Standalone,"
echo "(2)   MONC with CASIM,"
echo "(3)   MONC with SOCRATES,"
echo "(4)   MONC with CASIM and SOCRATES"
echo ""
echo "Select which option [1-4]: "
read compileoption

case $compileoption in
1)
  fcm make -j4 -f $fcmconfig
  ;;
2)
  fcm make -j4 -f $fcmconfig -f fcm-make/casim.cfg
  ;;
3)
  fcm make -j4 -f $fcmconfig -f fcm-make/socrates.cfg
  ;;
4)
  fcm make -j4 -f $fcmconfig -f fcm-make/casim_socrates.cfg
  ;;
*)
  echo "Unexpected compilation option. Should be an integer in the range 1-4"
  ;;
esac
