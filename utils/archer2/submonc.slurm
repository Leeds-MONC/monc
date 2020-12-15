#!/bin/bash
#SBATCH --partition=standard
#SBATCH --qos=standard

#SBATCH --job-name=monc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5

#SBATCH --time=00:10:00

# set up default ARCHER2 job environment
module load epcc-job-env

module list

if [ ! -d checkpoint_files ]; then mkdir checkpoint_files; fi
if [ ! -d diagnostic_files ]; then mkdir diagnostic_files; fi

# set env variables for submission command
config_path='tests/straka_short.mcf'
checkpoint_fn="checkpoint_files/straka_dump.nc"

srun ./build/bin/monc_driver.exe --config=$config_path  --checkpoint_file=$checkpoint_fn
