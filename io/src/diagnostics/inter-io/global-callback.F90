!> Global callback inter IO, which registers the callback with identifiers and then the procedure is actually called in a 
!! thread once all other IO servers have issued a matching call
module global_callback_inter_io_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use configuration_parser_mod, only : io_configuration_type
  use inter_io_specifics_mod, only : handle_completion
  use allreduction_inter_io_mod, only : init_allreduction_inter_io, finalise_allreduction_inter_io, perform_inter_io_allreduction
  implicit none

#ifndef TEST_MODE
  private
#endif

  public init_global_callback_inter_io, finalise_global_callback_inter_io, perform_global_callback
contains

  !> Initialises the global callback
  !! @param io_configuration The IO server configuration
  subroutine init_global_callback_inter_io(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    call init_allreduction_inter_io(io_configuration)
  end subroutine init_global_callback_inter_io
  
  !> Finalises the global callback
  !! @param io_configuration The IO server configuration
  subroutine finalise_global_callback_inter_io(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    call finalise_allreduction_inter_io(io_configuration)
  end subroutine finalise_global_callback_inter_io

  !> Performs a global callback
  !! @param io_configuration The IO server configuration
  !! @param field_name The field name identifier
  !! @param timestep The timestep identifier
  !! @param completion_procedure The completion procedure which is the callback itself
  subroutine perform_global_callback(io_configuration, field_name, timestep, completion_procedure)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: timestep
    character(len=*), intent(in) :: field_name
    procedure(handle_completion) :: completion_procedure

    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: v
    integer :: i

    allocate(v(1))
    v=0.0
    do i=1, io_configuration%number_of_moncs
      call perform_inter_io_allreduction(io_configuration, v, 1, field_name, 2, 0, timestep, completion_procedure)
    end do
    deallocate(v)
  end subroutine perform_global_callback    
end module global_callback_inter_io_mod
