#!/bin/bash
#SBATCH --partition=standard
#SBATCH --qos=standard

#SBATCH --job-name=monc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4

#SBATCH --time=00:10:00

echo Starting job

# set up default ARCHER2 job environment
module restore PrgEnv-cray
module load cpe-gnu
module load gcc/9.3.0
module load cray-netcdf-hdf5parallel
module load cray-hdf5-parallel
module load cray-fftw/3.3.8.7
module load petsc/3.13.3
module load atp
export ATP_ENABLED=1
ulimit -c unlimited

module list

if [ ! -d checkpoint_files ]; then mkdir checkpoint_files; fi
if [ ! -d monc_stdout ]; then mkdir monc_stdout; fi
if [ ! -d diagnostic_files ]; then mkdir diagnostic_files; fi

# set env variables for submission command
config_path='tests/straka_short.mcf'
checkpoint_fn="checkpoint_files/straka_dump.nc"

srun ./build/bin/monc_driver.exe --config=$config_path  --checkpoint_file=$checkpoint_fn

srun --unbuffered --cpu-bind=cores --distribution=block:block --hint=nomultithread ./build/bin/monc_driver.exe --config=$config_path  --checkpoint_file=$checkpoint_fn
