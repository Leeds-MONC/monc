#!/bin/bash

# Set the directory in which this script is stored
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
SET_ENV=${1:-$DIR/set_env_meto_desktop.sh}
CONFIG=${2:-https://puma.nerc.ac.uk/svn/MONC_svn/MONC/branches/dev/bshipw/r726_use_fcm_make/fcm-make/monc_x86_64_gfortran_opt.cfg}

unset MODULESHOME

. $SET_ENV

fcm make -j4 -f ${CONFIG}
