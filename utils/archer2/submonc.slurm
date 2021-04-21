#!/bin/bash
#SBATCH --job-name=MONC_straka_4core
#SBATCH --output=%x.o%j
  # %x gives job-name (SLURM_JOB_NAME)
  # %j gives jobid (individual SLURM_JOB_ID)
  # %A gives jobid (master     SLURM_ARRAY_JOB_ID)
  # %a gives array task id number
  #  https://slurm.schedmd.com/sbatch.html
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --account=n02-NEP013406
#SBATCH --partition=standard
#SBATCH --reservation=shortqos
#SBATCH --qos=short

echo Starting job

compiler=gnu
#compiler=cray

# set up default ARCHER2 job environment
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
export OMP_NUM_THREADS=1
ulimit -c unlimited

module list

if [ ! -d checkpoint_files ]; then mkdir checkpoint_files; fi
if [ ! -d monc_stdout ]; then mkdir monc_stdout; fi
if [ ! -d diagnostic_files ]; then mkdir diagnostic_files; fi

# set env variables for submission command
config_path='tests/straka_short.mcf'
checkpoint_fn="checkpoint_files/straka_dump.nc"

srun --unbuffered --cpu-bind=cores --distribution=block:block --hint=nomultithread ./build/bin/monc_driver.exe --config=$config_path  --checkpoint_file=$checkpoint_fn
