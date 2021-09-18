module shallow_convection_mod

  ! This module provides an example of tracer initialization for the BOMEX case
       
  use datadefn_mod, only : DEFAULT_PRECISION
  use logging_mod, only : LOG_INFO, LOG_ERROR, log_master_log
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type
  use grids_mod, only : vertical_grid_configuration_type, X_INDEX, Y_INDEX, Z_INDEX
  use conversions_mod, only : conv_to_string
  use optionsdatabase_mod, only : options_get_real, options_get_logical

  implicit none

#ifndef TEST_MODE
  private
#endif
  
  real(kind=DEFAULT_PRECISION) :: bl_height
  real(kind=DEFAULT_PRECISION) :: bl_tracer_value
  real(kind=DEFAULT_PRECISION) :: cl_height
  real(kind=DEFAULT_PRECISION) :: cl_tracer_value
  logical :: include_blob

  public shallow_convection_get_descriptor
contains

  type(component_descriptor_type) function shallow_convection_get_descriptor()
    shallow_convection_get_descriptor%name="shallow_convection"
    shallow_convection_get_descriptor%version=0.1
    shallow_convection_get_descriptor%initialisation=>initialisation_callback
    shallow_convection_get_descriptor%timestep=>timestep_callback
  end function shallow_convection_get_descriptor

  !> Initialise radioactive tracers: this is an example.
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state
    type(vertical_grid_configuration_type) :: vertical

    integer :: n, i_tracer, i, j, k, x_offset, y_offset
    REAL :: x_factor, y_factor ! x_pos = x_factor * (x_offset +i) etc.
    REAL :: x_pos, y_pos ! x and y in m.
    REAL :: x_size, y_size
    REAL, PARAMETER :: sz = 100.0, sy = 1000.0, sx = 1000.0 
    REAL :: rsq
    REAL :: init_time
    
    vertical=current_state%global_grid%configuration%vertical
    
    if (.not. current_state%continuation_run) then

    ! Read in parameters from options database
      bl_height=options_get_real(current_state%options_database, "bl_height") 
      bl_tracer_value=options_get_real(current_state%options_database, "bl_tracer_value") 
      cl_height=options_get_real(current_state%options_database, "cl_height") 
      cl_tracer_value=options_get_real(current_state%options_database, "cl_tracer_value") 
      include_blob=options_get_logical(current_state%options_database, "include_blob")
    
      if (current_state%timestep == 1) then
        init_time = current_state%time
      else
        init_time = current_state%time + current_state%dtm
      end if    

      call log_master_log(LOG_INFO, "[SHALLOW CONVECTION] Initialise tracers"//                         &
         " timestep: "//trim(conv_to_string(current_state%timestep))//                                  &
         " time: "//trim(conv_to_string(init_time,5)) )
    
      x_factor = current_state%global_grid%configuration%horizontal%dx  
      y_factor = current_state%global_grid%configuration%horizontal%dy  

      x_size = current_state%global_grid%size(X_INDEX) * x_factor
      y_size = current_state%global_grid%size(Y_INDEX) * y_factor
    
      do n = 1, current_state%n_radioactive_tracers
    
        i_tracer = n + current_state%radioactive_tracer_index - 1
        current_state%tracer(i_tracer)%data(:,:,:) = 0.0_DEFAULT_PRECISION
        current_state%ztracer(i_tracer)%data(:,:,:) = 0.0_DEFAULT_PRECISION
      
        ! Start at 0
        x_offset = current_state%local_grid%start(X_INDEX) - 1 - (1 + current_state%local_grid%halo_size(X_INDEX))
        y_offset = current_state%local_grid%start(Y_INDEX) - 1 - (1 + current_state%local_grid%halo_size(Y_INDEX))
        
        do i = current_state%local_grid%local_domain_start_index(X_INDEX) - current_state%local_grid%halo_size(X_INDEX), &
               current_state%local_grid%local_domain_end_index(X_INDEX)   + current_state%local_grid%halo_size(X_INDEX)
            
          x_pos = x_factor * (x_offset + i)
            
          do j = current_state%local_grid%local_domain_start_index(Y_INDEX) - current_state%local_grid%halo_size(Y_INDEX), &
                 current_state%local_grid%local_domain_end_index(Y_INDEX)   + current_state%local_grid%halo_size(Y_INDEX)

            y_pos = y_factor * (y_offset + j)
            
            do k = 1, current_state%global_grid%size(Z_INDEX)
            
              if (n .eq. 1) then
              
                ! Set tracer 1 to 1 in BL
                if (vertical%zn(k) .le. bl_height) then
                  current_state%tracer(i_tracer)%data(k, j, i) = bl_tracer_value
                  current_state%ztracer(i_tracer)%data(k, j, i) = bl_tracer_value
                end if
                               
              else if (n .eq. 2) then

                ! Set tracer 2 to 1 in cloud layer
                if (vertical%zn(k) .gt. bl_height .and. &
                    vertical%zn(k) .le. cl_height) then
                  current_state%tracer(i_tracer)%data(k, j, i) = bl_tracer_value
                  current_state%ztracer(i_tracer)%data(k, j, i) = bl_tracer_value
                end if
              
              else if (include_blob .and. (n .eq. 3)) then

                ! Set tracer 3 to blob in inversion layer
                rsq = ((vertical%zn(k) - 0.5 * (bl_height+cl_height))/sz)**2 + &
                      ((y_pos - 0.5 * y_size)/sy)**2 + ((x_pos - 0.5 * x_size)/sx)**2
                      
                current_state%tracer(i_tracer)%data(k, j, i) = EXP(-rsq/2.0)
                current_state%ztracer(i_tracer)%data(k, j, i) = EXP(-rsq/2.0)
              
              end if ! case n
              
            end do ! k
            
          end do ! j
          
        end do ! i
      
      end do ! n

    endif

  end subroutine initialisation_callback

  subroutine timestep_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

! Do nothing!
  end subroutine timestep_callback

end module shallow_convection_mod
