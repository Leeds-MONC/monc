!> Wrapper around the FFTW3 bindings. This will privatise by default all of FFTW3 apart from the calls that we explicitly need
!! in our pencil solver
module fftw_mod
  use, intrinsic :: iso_c_binding
  implicit none

#ifndef TEST_MODE
  private
#endif

  include 'fftw3.f03'

  public C_DOUBLE_COMPLEX, C_PTR, FFTW_BACKWARD, FFTW_FORWARD, FFTW_ESTIMATE, fftw_plan_many_dft_r2c, &
       fftw_plan_many_dft_c2r, fftw_plan_many_dft, fftw_execute_dft, fftw_execute_dft_c2r, fftw_execute_dft_r2c, fftw_destroy_plan
end module fftw_mod
