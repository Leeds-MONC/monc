!> Saturation physics functionality which is used throughout the code
module saturation_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  implicit none

   real(kind=DEFAULT_PRECISION), parameter :: &
        tk0c = 273.15_DEFAULT_PRECISION,      &! Temperature of freezing in Kelvin
        qsa1 = 3.8_DEFAULT_PRECISION,         &! Top in equation to calculate qsat
        qsa2 = -17.2693882_DEFAULT_PRECISION, &! Constant in qsat equation
        qsa3 = 35.86_DEFAULT_PRECISION,       &! Constant in qsat equation
        qsa4 = 6.109_DEFAULT_PRECISION       ! Constant in qsat equation

   public qsaturation, dqwsatdt
contains

  !> Function to return the saturation mixing ratio over water based on tetans formular
  !! QS=3.8/(P*EXP(-17.2693882*(T-273.15)/(T-35.86))-6.109)
  !! @param temperature The temperature
  !! @param pressure The pressure
  !! @returns The mixing ratio
  real(kind=DEFAULT_PRECISION) function qsaturation(temperature, pressure)
    real(kind=DEFAULT_PRECISION), intent(in) :: temperature, pressure

    qsaturation = qsa1/(pressure*exp(qsa2*(temperature - tk0c)/(temperature - qsa3)) - qsa4)
  end function qsaturation  

  !> Calculated the rate of change with temperature of saturation mixing ratio over liquid water. Based on tetans formular
  !! QS=3.8/(P*EXP(-17.2693882*(T-273.15)/(T-35.86))-6.109)
  !! @param saturation_mixing_ratio The saturation mixing ratio
  !! @param temperature The temperature  
  !! @returns The mixing ratio
  real(kind=DEFAULT_PRECISION) function dqwsatdt(saturation_mixing_ratio, temperature)
    real(kind=DEFAULT_PRECISION), intent(in) :: saturation_mixing_ratio, temperature
    
    dqwsatdt = -qsa2*(tk0c-qsa3)*(1.0_DEFAULT_PRECISION+qsa4*saturation_mixing_ratio/qsa1)*&
         saturation_mixing_ratio*(temperature-qsa3)**(-2.0_DEFAULT_PRECISION)
  end function dqwsatdt  
end module saturation_mod
