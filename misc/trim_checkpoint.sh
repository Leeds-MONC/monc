#!/bin/bash

# ------------------------------------------------------------------------------------- #
#  Routine to remove stored diagnostic data (serialised fields) and configuration 
#    information from checkpoints.
#  Good for saving space (if the eliminated data is not needed).
#  Good for creating tidy reconfiguration starter checkpoints.
#
#  Script tested on Monsoon.  
#
#  Prerequisites:
#    nco module available
#    script expected to be run from directory containing checkpoint file to be trimmed
#  Input:
#    checkpoint file name from working directory
#  Output:
#    trimmed file written to working directory
#
#  Example usage:
#    host:.../checkpoint_files> ../misc/trim_checkpoint.sh RCE_dump_7950.nc
#    RESULT: RCE_dump_trimmed_7950.nc
# ------------------------------------------------------------------------------------- #


# Load nco
#module load nco # (available on JASMIN)
module load nco/4.6.7-netcdf_4.4.1.1.3 # better, available on Monsoon

# Receive input
data=$1     # checkpoint file name from an old run

# Create trimmed checkpoint name
stamp=`echo $data | rev | cut -f 2 -d '.' | cut -f 1 -d '_' | rev` 
ndata=${data/$stamp/trimmed_$stamp} 

# Report
echo "Trimming $data   to   $ndata"

# Set up variable exclusion list (those not to be overwritten)
  # The configuration fields
excluded=io_configuration,options_database
  # The serial terms (not always the same/present, so handle on fly)
serial=`ncks -m $data | grep 'type' | cut -f 1 -d ':' | grep 'serialised_'`
for snc in $serial 
do 
  excluded=$snc,$excluded
done

# Report
echo "Excluding ${excluded}"

# (-A)ppend all terms e(-x)cept the (-v)ariables in the $excluded list
# 
ncks -A -x -v $excluded $data $ndata #> quiet ; rm quiet

# Make sure it has a 'created' global attribute
cdate=`date`
ncatted -h -a created,global,c,c,"${cdate}" $ndata

# All done
echo "Good to go!"

exit 0
