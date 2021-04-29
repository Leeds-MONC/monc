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

#Script to compile Monc on ARC4:

module purge
module load user
module switch intel gnu
module switch openmpi mvapich2
module load fftw netcdf hdf5 fcm
. /nobackup/cemac/cemac.sh

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
  fcm make -j4 -f fcm-make/monc-arc4-gnu.cfg
  ;;
2)
  fcm make -j4 -f fcm-make/monc-arc4-gnu.cfg -f fcm-make/casim.cfg
  ;;
3)
  fcm make -j4 -f fcm-make/monc-arc4-gnu.cfg -f fcm-make/socrates.cfg
  ;;
4)
  fcm make -j4 -f fcm-make/monc-arc4-gnu.cfg -f fcm-make/casim_socrates.cfg
  ;;
*)
  echo "Unexpected compilation option. Should be an integer in the range 1-4"
  ;;
esac
