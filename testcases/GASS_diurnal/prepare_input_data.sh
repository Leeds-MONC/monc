#!/bin/bash

#  Script works on Monsoon.  
#  Input only input data file and base day value
#  Output written to working directory.


#module load nco
module load nco/4.6.7-netcdf_4.4.1.1.3 # better, available on Monsoon


data=$1     # input GASS data file
base_day=$2 # day of year value to become zero time

stamp=`echo $data | rev | cut -f 2- -d '.' | rev `
ndata=${stamp}_reformed_${base_day}.nc 

echo "Upgrading $data   to   $ndata"

ncap2 -s 'lev=lev*100;theta_tendency=(T_adv_h+s_adv_v)/3600.;q_tendency=(q_adv_h+q_adv_v)/3600./1000.;surface_temperature=T_skin+273.15;' $data temp_$ndata
[ -f $ndata ] && rm $ndata
[ -f B$ndata ] && rm B$ndata
ncks -O -v time,lev,theta_tendency,q_tendency,surface_temperature,SH,LH,RH_srf temp_$ndata $ndata
ncrename -v RH_srf,surface_humidity -v SH,surface_sensible_heat_flux -v LH,surface_latent_heat_flux $ndata

ncatted -h -a long_name,time,o,c,"time since - day ${base_day} -" $ndata
ncatted -h -a units,time,o,c,"seconds since - day ${base_day} -" $ndata

ncatted -h -a units,lev,o,c,"Pa" $ndata

ncatted -h -a long_name,q_tendency,o,c,"total advective water vapour tendency" $ndata
ncatted -h -a units,q_tendency,o,c,"kg/kg/s" $ndata

ncatted -h -a units,surface_temperature,o,c,"K" $ndata

ncatted -h -a units,theta_tendency,o,c,"K/s" $ndata



ncap2 -O -s 'time=(time-'${base_day}')*86400.;' $ndata B$ndata



rm temp_$ndata
mv B$ndata $ndata

#All done
echo "Good to go!"

exit 0
