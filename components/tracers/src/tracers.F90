!> Description of what the component does.
!! Add as much information as is sensible.
module tracers_mod

  !> tracers source code for a new component
       
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH     
  use logging_mod, only : LOG_INFO, LOG_ERROR, log_master_log
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
       options_get_integer_array, options_get_real_array, options_get_string
  use conversions_mod, only : conv_to_string
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX

  implicit none

#ifndef TEST_MODE
  private
#endif
  !> Parameter values which determine type of surface bc for tracers
  integer, parameter, public :: TRACER_SURFACE_FLUX_FROM_DECAY=0, TRACER_SURFACE_FLUX_SPECIFIED=1, TRACER_SURFACE_VALUE_SPECIFIED=2
  logical, public :: tracers_enabled=.false., trajectories_enabled=.false., radioactive_tracers_enabled=.false.
  integer, public :: traj_interval 

  public tracers_get_descriptor, reinitialise_trajectories, get_tracer_name, get_tracer_options
  
contains

  !> Provides a description of this component for the core to register
  !! @returns The descriptor containing registration information for this component
  type(component_descriptor_type) function tracers_get_descriptor()
    tracers_get_descriptor%name="tracers"
    tracers_get_descriptor%version=0.1
    tracers_get_descriptor%initialisation=>initialisation_callback
    tracers_get_descriptor%timestep=>timestep_callback
    tracers_get_descriptor%finalisation=>finalisation_callback
  end function tracers_get_descriptor

  !> Initialisation callback hook which will do nothing.
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state
                      
  end subroutine initialisation_callback

  !> Timestep callback hook which will <DO WHAT?>
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state
    if (current_state%n_radioactive_tracers >0) then
      call tracer_decay(current_state)
    end if
  end subroutine timestep_callback

  !> Finalisation callback hook which will <DO WHAT?>
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

  end subroutine finalisation_callback

  !> Compute decay of radioactive tracers.
  !! @param current_state The current model state
  subroutine tracer_decay(current_state)
    type(model_state_type), intent(inout), target :: current_state

    REAL(kind=DEFAULT_PRECISION) :: kdecay
    REAL(kind=DEFAULT_PRECISION), PARAMETER :: sec_in_hour=3600.0

    integer :: current_y_index, current_x_index, n

    current_y_index=current_state%column_local_y
    current_x_index=current_state%column_local_x

    ! Decay of radioactive tracers
    if (current_state%n_radioactive_tracers .gt. 0) then
      do n=current_state%radioactive_tracer_index, current_state%radioactive_tracer_index + current_state%n_radioactive_tracers - 1 
        kdecay = current_state%tracer_decay_rate(n-current_state%radioactive_tracer_index+1) / sec_in_hour
        current_state%stracer(n)%data(:, current_y_index, current_x_index) = &
          current_state%stracer(n)%data(:, current_y_index, current_x_index) - &
          kdecay * current_state%tracer(n)%data(:, current_y_index, current_x_index)
      end do
    end if

  end subroutine tracer_decay
  

  !> Reinitialise trajectory tracers to position in model.
  !! @param current_state The current model state
  subroutine reinitialise_trajectories(current_state)
    use science_constants_mod, only : pi
    type(model_state_type), intent(inout), target :: current_state
    integer :: i, j, k, offset
    REAL(kind=DEFAULT_PRECISION) :: phase_factor
    REAL(kind=DEFAULT_PRECISION),dimension(current_state%local_grid%size(Z_INDEX)+current_state%local_grid%halo_size(Z_INDEX)*2, &
                                           current_state%local_grid%size(Y_INDEX)+current_state%local_grid%halo_size(Y_INDEX)*2, &
                                           current_state%local_grid%size(X_INDEX)+current_state%local_grid%halo_size(X_INDEX)*2) &
                                           :: diffr_copy, diffi_copy
      
    LOGICAL, PARAMETER :: usecyclic_x=.TRUE., usecyclic_y=.TRUE., usezn=.FALSE.
    INTEGER :: ixr, ixi, iyr, iyi, iz
    REAL :: reinit_time
    
    if (current_state%timestep == 1) then
      reinit_time = current_state%time
    else
      reinit_time = current_state%time + current_state%dtm
    end if    

    call log_master_log(LOG_INFO, "[TRACERS] Reinitialise trajectories"//                  &
       " timestep: "//trim(conv_to_string(current_state%timestep))//                       &
       " time: "//trim(conv_to_string(reinit_time,5)) )
      
    ixr = current_state%traj_tracer_index
    ixi = ixr+1 
    iyr = ixr+2 
    iyi = ixr+3 
    iz  = ixr+4
              
    if (current_state%tracer(ixr)%active) then
      diffr_copy = current_state%ztracer(ixr)%data(:,:,:)-current_state%tracer(ixr)%data(:,:,:)
      diffi_copy = current_state%ztracer(ixi)%data(:,:,:)-current_state%tracer(ixi)%data(:,:,:)
      ! Start at 0
      offset = current_state%local_grid%start(X_INDEX)-1-(1 + current_state%local_grid%halo_size(X_INDEX))
      if (usecyclic_x) then
        phase_factor = 2.0_DEFAULT_PRECISION * pi / current_state%global_grid%size(X_INDEX)
        do i=current_state%local_grid%local_domain_start_index(X_INDEX) - current_state%local_grid%halo_size(X_INDEX),&
            current_state%local_grid%local_domain_end_index(X_INDEX) + current_state%local_grid%halo_size(X_INDEX)
          current_state%tracer(ixr)%data(:,:,i) = COS( phase_factor * ( offset + i ) )
          current_state%ztracer(ixr)%data(:,:,i) = COS( phase_factor * ( offset + i ) )
          current_state%tracer(ixi)%data(:,:,i) = SIN( phase_factor * ( offset + i ) )
          current_state%ztracer(ixi)%data(:,:,i) = SIN( phase_factor * ( offset + i ) )
         end do
!        end if
      else
        ! Start at 0
        offset = current_state%local_grid%start(X_INDEX) - 1  - current_state%local_grid%halo_size(X_INDEX)- 1 
        do i = current_state%local_grid%local_domain_start_index(X_INDEX) - current_state%local_grid%halo_size(X_INDEX),&
               current_state%local_grid%local_domain_end_index(X_INDEX) + current_state%local_grid%halo_size(X_INDEX)
          current_state%tracer(ixr)%data(:,:,i) = offset + i 
          current_state%ztracer(ixr)%data(:,:,i) = offset + i 
          current_state%tracer(ixi)%data(:,:,i) = offset + i 
          current_state%ztracer(ixi)%data(:,:,i) = offset + i 
        end do
      end if
      current_state%ztracer(ixr)%data(:,:,:) = current_state%tracer(ixr)%data(:,:,:) + diffr_copy
      current_state%ztracer(ixi)%data(:,:,:) = current_state%tracer(ixi)%data(:,:,:) + diffi_copy
!        end if
    end if
    
    if (current_state%tracer(iyr)%active) then
!        if (.not. current_state%continuation_run) then
      diffr_copy = current_state%ztracer(iyr)%data(:,:,:)-current_state%tracer(iyr)%data(:,:,:)
      diffi_copy = current_state%ztracer(iyi)%data(:,:,:)-current_state%tracer(iyi)%data(:,:,:)
      ! Start at 0
      offset = current_state%local_grid%start(Y_INDEX) -1 - (1 + current_state%local_grid%halo_size(Y_INDEX))
      if (usecyclic_y) then
!        offset = -current_state%local_grid%start(Y_INDEX) 
        phase_factor = 2.0_DEFAULT_PRECISION * pi / current_state%global_grid%size(Y_INDEX)
        do j=current_state%local_grid%local_domain_start_index(Y_INDEX) - current_state%local_grid%halo_size(Y_INDEX), &
             current_state%local_grid%local_domain_end_index(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX)
          current_state%tracer(iyr)%data(:,j,:) = COS( phase_factor * ( offset + j ) )
          current_state%ztracer(iyr)%data(:,j,:) = COS( phase_factor * ( offset + j ) )
          current_state%tracer(iyi)%data(:,j,:) = SIN( phase_factor * ( offset + j ) )
          current_state%ztracer(iyi)%data(:,j,:) = SIN( phase_factor * ( offset + j ) )
        end do
      else
        do j=current_state%local_grid%local_domain_start_index(Y_INDEX) - current_state%local_grid%halo_size(Y_INDEX), &
             current_state%local_grid%local_domain_end_index(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX)
          current_state%tracer(iyr)%data(:,j,:) = offset + j
          current_state%ztracer(iyr)%data(:,j,:) = offset + j
          current_state%tracer(iyi)%data(:,j,:) = offset + j 
          current_state%ztracer(iyi)%data(:,j,:) = offset + j 
        end do
      end if
      current_state%ztracer(iyr)%data(:,:,:) = current_state%tracer(iyr)%data(:,:,:) + diffr_copy
      current_state%ztracer(iyi)%data(:,:,:) = current_state%tracer(iyi)%data(:,:,:) + diffi_copy
!        end if
    end if
    
    if (current_state%tracer(iz)%active)then
!        if (.not. current_state%continuation_run) then
      diffr_copy = current_state%ztracer(iz)%data(:,:,:)-current_state%tracer(iz)%data(:,:,:)
      if (usezn) then      
        do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
          do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)  
            current_state%tracer(iz)%data(:,j,i) = current_state%global_grid%configuration%vertical%zn(:)
            current_state%ztracer(iz)%data(:,j,i) = current_state%global_grid%configuration%vertical%zn(:)
          end do
        end do
      else
        ! start at 0 
        offset = -1
        do k = 1, current_state%local_grid%size(Z_INDEX)           
          current_state%tracer(iz)%data(k,:,:) = offset + k
          current_state%ztracer(iz)%data(k,:,:) = offset + k
        end do
      end if
      current_state%ztracer(iz)%data(:,:,:) = current_state%tracer(iz)%data(:,:,:) + diffr_copy
!        end if
    end if 
    
  end subroutine reinitialise_trajectories
  
  !> Get tracer options.
  !! @param current_state The current model state
  subroutine get_tracer_options(current_state)
    type(model_state_type), intent(inout), target :: current_state
    integer :: i_tracer
    current_state%n_tracers = 0
    current_state%traj_tracer_index = 0
    current_state%n_radioactive_tracers = 0    
    tracers_enabled = options_get_logical(current_state%options_database, "tracers_enabled")
    
    if (tracers_enabled) then
    
      call log_master_log(LOG_INFO, "Tracers enabled.")
    
      trajectories_enabled  = options_get_logical(current_state%options_database, "trajectories_enabled")
      if (trajectories_enabled) then
        call log_master_log(LOG_INFO, "Trajectories enabled.")    
        current_state%n_tracers = 5
        current_state%traj_tracer_index = 1
        ! Obtain the output traj_interval value
        traj_interval = nint(options_get_real(current_state%options_database, &
                             options_get_string(current_state%options_database,"traj_interval")))
      end if ! trajectories_enabled
      
      radioactive_tracers_enabled   = options_get_logical(current_state%options_database, "radioactive_tracers_enabled")

      if (radioactive_tracers_enabled) then
      
        call log_master_log(LOG_INFO, "Radioactive Tracers enabled.")    
        current_state%n_radioactive_tracers = options_get_integer(current_state%options_database, "n_radioactive_tracers")

        if (current_state%n_radioactive_tracers > 0) then
        
          current_state%radioactive_tracer_index = current_state%n_tracers + 1
          current_state%n_tracers = current_state%n_tracers + current_state%n_radioactive_tracers    
        
          allocate(current_state%tracer_decay_rate(current_state%n_radioactive_tracers))
          call options_get_real_array(current_state%options_database,"tracer_decay_rate",current_state%tracer_decay_rate)
          allocate(current_state%tracer_surf_bc_opt(current_state%n_radioactive_tracers))
          current_state%tracer_surf_bc_opt(:) = 0
          call options_get_integer_array(current_state%options_database, "tracer_surface_bc_option", &
            current_state%tracer_surf_bc_opt)          
          allocate(current_state%tracer_surf_bc(current_state%n_radioactive_tracers))
          call options_get_real_array(current_state%options_database,"tracer_surface_bc",current_state%tracer_surf_bc)
          if (current_state%parallel%my_rank == 0 ) then
            do i_tracer = 1,current_state%n_radioactive_tracers
              call log_master_log(LOG_INFO, "Tracer "//trim(conv_to_string(i_tracer))// &
                " decay rate = "//trim(conv_to_string(current_state%tracer_decay_rate(i_tracer),5))// &
                " surface bc type = "//trim(conv_to_string(current_state%tracer_surf_bc_opt(i_tracer)))// &
                " surface bc value = "//trim(conv_to_string(current_state%tracer_surf_bc(i_tracer),5)))
              if (current_state%tracer_surf_bc_opt(i_tracer) < 0 .or. current_state%tracer_surf_bc_opt(i_tracer) > 2) then
                call log_master_log(LOG_ERROR, "Radioactive tracer with illegal surface bc option.")          
              end if
            end do
          end if
        else
        
            call log_master_log(LOG_ERROR, "Cannot run with less than 1 radioactive tracer with radioactive tracers enabled")          
   
        end if ! n_radioactive_tracers > 0
                  
      end if ! radioactive_tracers_enabled
      
      if (current_state%n_tracers <= 0) then
        call log_master_log(LOG_ERROR, "Cannot run with less than 1 tracer with Tracers enabled")
      end if      
                    
    end if ! tracers_enabled   
  end subroutine get_tracer_options 
  
  function get_tracer_name(i, traj_tracer_index, radioactive_tracer_index, n_radioactive_tracers, n_tracers)
    integer, intent(in) :: i, traj_tracer_index, radioactive_tracer_index, n_radioactive_tracers, n_tracers
    character(len=STRING_LENGTH) :: get_tracer_name
    character(len=2), dimension(5), parameter :: traj_name = ["xr","xi","yr","yi","zr"] 
    if (traj_tracer_index .gt. 0 .and. i .lt. (traj_tracer_index+5)) then
      get_tracer_name = "traj_"//traj_name(i - traj_tracer_index + 1)
    else if (n_radioactive_tracers .gt. 0 .and. i .ge. radioactive_tracer_index .and. &
             i .lt. (radioactive_tracer_index + n_radioactive_tracers)) then
      get_tracer_name = "rad"//trim(conv_to_string(i-radioactive_tracer_index+1))
    else
      get_tracer_name = "_"//trim(conv_to_string(i))
    end if
  end function get_tracer_name
  
end module tracers_mod
