module naming_conventions_mod

  implicit none

#ifndef TEST_MODE
  private
#endif

  character(len=5), parameter :: k_per_day='K/day'
  character(len=3), parameter :: k_per_second='K/s'
  character(len=7), parameter :: m_per_second_per_day='m/s/day'
  character(len=5), parameter :: m_per_second_per_second='m/s/s'
  character(len=7), parameter :: kg_per_kg_per_second='kg/kg/s'
  character(len=9), parameter :: kg_per_kg_per_day='kg/kg/day'
  character(len=8), parameter :: g_per_kg_per_day='g/kg/day'
  character(len=6), parameter :: g_per_kg_per_second='g/kg/s'
  character(len=7), parameter :: degC='celsius'

  public k_per_day, k_per_second, m_per_second_per_day, m_per_second_per_second, &
     kg_per_kg_per_second, kg_per_kg_per_day, g_per_kg_per_day, g_per_kg_per_second, degC

end module naming_conventions_mod
