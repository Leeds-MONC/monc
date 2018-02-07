module setfluxlook_mod
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : PRESCRIBED_SURFACE_FLUXES, PRESCRIBED_SURFACE_VALUES, model_state_type
  use grids_mod, only : Z_INDEX
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use optionsdatabase_mod, only : options_get_real, options_get_integer, options_get_array_size, &
     options_get_real_array, options_get_string, options_get_logical
  use saturation_mod, only : qsaturation
  use collections_mod, only : map_type
  use science_constants_mod
  use netcdf, only : nf90_noerr, nf90_global, nf90_nowrite, nf90_inquire_attribute, nf90_open, nf90_strerror, &
       nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, nf90_get_var, nf90_inquire, nf90_close, nf90_get_att
  use logging_mod, only : LOG_INFO, LOG_WARN, LOG_ERROR, LOG_DEBUG, log_master_log, log_log, log_get_logging_level
  use q_indices_mod, only: get_q_index, standard_q_names
  use interpolation_mod, only: interpolate_point_linear_1d
  use naming_conventions_mod

  implicit none

#ifndef TEST_MODE
  private
#endif

  character(len=*), parameter ::                              &
     TIME_KEY =                  "time",                      & !<  NetCDF data time key
     SURFACE_TEMPERATURES_KEY =  "surface_temperature",       & !<  NetCDF data surface_temperatures
     SURFACE_HUMIDITIES_KEY   =  "surface_humidity",          & !<  NetCDF data surface_humidities
     SURFACE_SHF_KEY          =  "surface_sensible_heat_flux",& !<  NetCDF data surface_sensible_heat_flux
     SURFACE_LHF_KEY          =  "surface_latent_heat_flux"     !<  NetCDF data surface_latent_heat_flux
     
  integer, parameter :: LOOKUP_ENTRIES = 80    !< Number of entries for MO lookup tables
  integer, parameter :: MAX_FILE_LEN=200       !< Maximum length of surface condition input filename
  integer, parameter :: MAX_SURFACE_INPUTS=750  !< Specifies the maximum number of surface inputs through configuration file
                                               !! Inputs through netcdf files are not limitted by this.
  
  character(MAX_FILE_LEN) :: input_file

  character(len=STRING_LENGTH) :: units_surface_temp='unset'  ! units of theta variable forcing

  real(kind=DEFAULT_PRECISION) :: max_change_buoyancy_flux
  real(kind=DEFAULT_PRECISION), allocatable :: surface_boundary_input_times(:)
  real(kind=DEFAULT_PRECISION), allocatable :: surface_sensible_heat_flux(:)
  real(kind=DEFAULT_PRECISION), allocatable :: surface_latent_heat_flux(:)
  real(kind=DEFAULT_PRECISION), allocatable :: surface_temperatures(:)
  real(kind=DEFAULT_PRECISION), allocatable :: surface_humidities(:)

  integer :: iqv  ! index for vapour

  public setfluxlook_get_descriptor

contains

  !> Descriptor of this component for registration
  !! @returns The diverr component descriptor
  type(component_descriptor_type) function setfluxlook_get_descriptor()
    setfluxlook_get_descriptor%name="setfluxlook"
    setfluxlook_get_descriptor%version=0.1
    setfluxlook_get_descriptor%initialisation=>initialisation_callback
    setfluxlook_get_descriptor%timestep=>timestep_callback
  end function setfluxlook_get_descriptor

  subroutine initialisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state
    
    current_state%lookup_table_entries=LOOKUP_ENTRIES
    call read_configuration(current_state)

    if (.not. current_state%passive_q .and. current_state%number_q_fields > 0)then
      if (.not. allocated(current_state%cq))then
        allocate(current_state%cq(current_state%number_q_fields))
        current_state%cq=0.0_DEFAULT_PRECISION
      end if
      iqv = get_q_index(standard_q_names%VAPOUR, 'setfluxlook')
      current_state%cq(iqv) = ratio_mol_wts-1.0
    end if

    if (current_state%use_surface_boundary_conditions)then
    if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_FLUXES) then
       if (current_state%use_time_varying_surface_values) then
          call set_flux(current_state)
       else
          current_state%surface_temperature_flux=surface_sensible_heat_flux(1) &
               /(current_state%global_grid%configuration%vertical%rho(1)*cp)
          current_state%surface_vapour_flux = surface_latent_heat_flux(1) &
               /(current_state%global_grid%configuration%vertical%rho(1)*rlvap)
       end if
       
       current_state%fbuoy=0.                                                                 
       if(.not. current_state%passive_th) current_state%fbuoy=&
            current_state%global_grid%configuration%vertical%buoy_co(1)*current_state%surface_temperature_flux
       if(.not. current_state%passive_q .and. current_state%number_q_fields > 0)then                                                      
          current_state%fbuoy=current_state%fbuoy+current_state%cq(iqv)*current_state%surface_vapour_flux*G
       end if
       call set_look(current_state)  ! _set M-O lookup table
       current_state%theta_surf=0.0_DEFAULT_PRECISION
       current_state%surface_vapour_mixing_ratio=0.0_DEFAULT_PRECISION
       !
    else if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_VALUES) then
       !
       ! Prescribed surface temperatures
       if (current_state%use_time_varying_surface_values) then
          call set_flux(current_state)
       else
       ! If surface_values are constant then surface_temperatures prescribed in 
       ! config and read in read_configuration but if humidity not set then
       ! surface vapour (surface_vapour_mixing_ratio) set to saturated value (see read_config)
          if (current_state%saturated_surface)then
             current_state%surface_vapour_mixing_ratio = qsaturation(surface_temperatures(1),current_state%surface_pressure*0.01)
          else 
             current_state%surface_vapour_mixing_ratio =  &
                  options_get_real(current_state%options_database, "surface_vapour_mixing_ratio")
          endif
          
       ! The code below copied from set_flux as these values need to be 
       ! set for both time varying and constant surface values
       ! Set theta_v
          current_state%theta_surf = surface_temperatures(1)*&
               (current_state%surface_reference_pressure/current_state%surface_pressure)**r_over_cp
          current_state%theta_virtual_surf = current_state%theta_surf
          if (.not. current_state%passive_q .and. current_state%number_q_fields > 0) then
             current_state%theta_virtual_surf = current_state%theta_surf +  &
                  current_state%global_grid%configuration%vertical%thref(2)*  &
                  current_state%cq(iqv)*current_state%surface_vapour_mixing_ratio
          end if          

          ! Finally set up new values of THVSURF dependent constants
          current_state%cmbc=betam*current_state%global_grid%configuration%vertical%zn(2)*G*&
               von_karman_constant/current_state%theta_virtual_surf
          current_state%rcmbc=1.0_DEFAULT_PRECISION/current_state%cmbc
          current_state%ellmocon=current_state%theta_virtual_surf/(G*von_karman_constant)
       end if
    end if
    end if

  end subroutine initialisation_callback

  subroutine timestep_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    if (current_state%use_surface_boundary_conditions)then
      if (current_state%use_time_varying_surface_values)then
        call set_flux(current_state)
        call change_look(current_state)
      end if
    end if
  end subroutine timestep_callback


  subroutine set_look(current_state)
    type(model_state_type), intent(inout), target :: current_state

    integer I, ik, iters(current_state%lookup_table_entries)    ! Loop counters
    real(kind=DEFAULT_PRECISION) :: smth, &       ! Relaxation parameter for unstable case
         ob, x1, x0

    if (current_state%fbuoy .le. 0.0_DEFAULT_PRECISION) return
    current_state%velmax=100.0_DEFAULT_PRECISION
    current_state%velmin=0.1_DEFAULT_PRECISION
    current_state%aloginv=1.0_DEFAULT_PRECISION/log(current_state%velmax/current_state%velmin)
    smth=0.1_DEFAULT_PRECISION  ! _relaxation parameter for unstable case
    do ik=1, current_state%lookup_table_entries
      current_state%lookup_table_velocity(ik)=current_state%velmin*(current_state%velmax/current_state%velmin)**&
           (real(ik-1)/real(current_state%lookup_table_entries-1))            
      current_state%lookup_table_ustr(ik)=current_state%lookup_table_velocity(ik)*&
           current_state%global_grid%configuration%vertical%vk_on_zlogm
      if (current_state%fbuoy .gt. 0.0_DEFAULT_PRECISION) then         ! _unstable                                  
        iters(ik)=0
        do i=1, 30      ! @check how many iterations needed!!
          iters(ik)=iters(ik)+1
          ob=-current_state%lookup_table_ustr(ik)**3/(von_karman_constant*current_state%fbuoy)
          x1=sqrt(sqrt(1.-gammam*(current_state%global_grid%configuration%vertical%zn(2)+z0)/ob))
          x0=sqrt(sqrt(1.-gammam*z0/ob))
          current_state%lookup_table_ustr(ik)=(1.-smth)*current_state%lookup_table_ustr(ik) + smth*&
               (von_karman_constant*current_state%lookup_table_velocity(ik)) / &
               (current_state%global_grid%configuration%vertical%zlogm-(2.0_DEFAULT_PRECISION*&
               log((x1+1.0_DEFAULT_PRECISION)/(x0+1.0_DEFAULT_PRECISION)) + log((x1*x1+1.0_DEFAULT_PRECISION)/&
               (x0*x0+1.0_DEFAULT_PRECISION)) + 2.0_DEFAULT_PRECISION*atan(x0) -2.0_DEFAULT_PRECISION*atan(x1)))
        end do
      end if
    end do
    current_state%cneut=current_state%lookup_table_ustr(current_state%lookup_table_entries)/&
         current_state%lookup_table_velocity(current_state%lookup_table_entries)
    current_state%cfc=current_state%lookup_table_ustr(1)*current_state%lookup_table_velocity(1)**convective_limit  ! _Businger-Dyer
  end subroutine set_look

  subroutine change_look(current_state)
    type(model_state_type), intent(inout), target :: current_state

    real(kind=DEFAULT_PRECISION):: crelax,&     ! Relaxation parameter to move to new USTR
         dfb,&        ! Iterative difference in buoyancy flux
         rnit,&    ! real number of iterations required
         velnew,&
         dvelustr(current_state%lookup_table_entries), ob, x1, x0
    integer n, ik, &   ! Loop counters
         nit     ! Number of iterations
    if (current_state%fbuoynew .le. 0.0_DEFAULT_PRECISION) then
      current_state%fbuoy=current_state%fbuoynew
      return
    end if
    if (max_change_buoyancy_flux .gt. 0.0_DEFAULT_PRECISION) then
      rnit = abs(current_state%fbuoynew-current_state%fbuoy)/max_change_buoyancy_flux
      nit = 1+int(rnit)
    else
      nit=1
    end if              
    if (nit .gt. 10 .or. (current_state%fbuoynew .gt. 0.0_DEFAULT_PRECISION .and. &
         current_state%fbuoy .le. 0.0_DEFAULT_PRECISION)) then
      current_state%fbuoy=current_state%fbuoynew
      call set_look(current_state)
      return                                                                   
    end if
    crelax=1./sqrt(real(nit))    ! maybe better with crelax=1.
    dfb=(current_state%fbuoynew-current_state%fbuoy)/real(nit)
    do n=1, nit
      current_state%fbuoy=current_state%fbuoy+dfb                                                          
      do ik=2, current_state%lookup_table_entries-1
        dvelustr(ik)=log(current_state%lookup_table_velocity(ik+1)/current_state%lookup_table_velocity(ik-1))&
             /log(current_state%lookup_table_ustr(ik+1)/current_state%lookup_table_ustr(ik-1))
      end do
      dvelustr(1)=log(current_state%lookup_table_velocity(2)/current_state%lookup_table_velocity(1))/&
           log(current_state%lookup_table_ustr(2)/current_state%lookup_table_ustr(1))
      dvelustr(current_state%lookup_table_entries)=log(current_state%lookup_table_velocity(current_state%lookup_table_entries)/&
           current_state%lookup_table_velocity(current_state%lookup_table_entries-1))/&
           log(current_state%lookup_table_ustr(current_state%lookup_table_entries)/&
           current_state%lookup_table_ustr(current_state%lookup_table_entries-1))
      do ik=1, current_state%lookup_table_entries        
        ! compute new mean vel. based on old ustar and new fbuoy        
        ob=-current_state%lookup_table_ustr(ik)**3/(von_karman_constant*current_state%fbuoy)
        x1=sqrt(sqrt(1.-gammam*(current_state%global_grid%configuration%vertical%zn(2)+z0)/ob))
        x0=sqrt(sqrt(1.-gammam*z0/ob))
        velnew=(current_state%lookup_table_ustr(ik)/von_karman_constant)*(current_state%global_grid%configuration%vertical%zlogm-&
             (2.0_DEFAULT_PRECISION*log((x1+1.0_DEFAULT_PRECISION)/(x0+1.0_DEFAULT_PRECISION)) + &
             log((x1*x1+1.0_DEFAULT_PRECISION)/(x0*x0+1.0_DEFAULT_PRECISION)) + 2.0_DEFAULT_PRECISION*atan(x0) &
             -2.0_DEFAULT_PRECISION*atan(x1)))        
        ! relax to new ustar        
        current_state%lookup_table_ustr(ik)=current_state%lookup_table_ustr(ik)/&
             (velnew/current_state%lookup_table_velocity(ik))**(crelax/dvelustr(ik))
      end do
    end do
    current_state%cneut=current_state%lookup_table_ustr(current_state%lookup_table_entries)/&
         current_state%lookup_table_velocity(current_state%lookup_table_entries)
    current_state%cfc=current_state%lookup_table_ustr(1)*current_state%lookup_table_velocity(1)**convective_limit  ! _Businger-Dyer
    current_state%fbuoy=current_state%fbuoynew
  end subroutine change_look

  subroutine set_flux(current_state)
    type(model_state_type), intent(inout), target :: current_state

    real(kind=DEFAULT_PRECISION) :: surface_temp   ! Surface temperature
  
    if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_FLUXES) then  ! Prescribed surface fluxes
      
      ! Linear interpolation of input data...
      call interpolate_point_linear_1d(surface_boundary_input_times, & 
         surface_sensible_heat_flux/(current_state%global_grid%configuration%vertical%rho(1)*cp), &
         current_state%time, current_state%surface_temperature_flux, &
         extrapolate='constant') 
      
      call interpolate_point_linear_1d(surface_boundary_input_times, & 
         surface_latent_heat_flux/(current_state%global_grid%configuration%vertical%rho(1)*rlvap), &
         current_state%time, current_state%surface_vapour_flux, &
         extrapolate='constant')

      ! Update buoyancy flux...
      current_state%fbuoynew=0.0_DEFAULT_PRECISION
      if (.not. current_state%passive_th) current_state%fbuoynew=&
         current_state%global_grid%configuration%vertical%buoy_co(1)*current_state%surface_temperature_flux
      if (.not. current_state%passive_q) then                                                      
        current_state%fbuoynew=current_state%fbuoynew+current_state%cq(iqv)*current_state%surface_vapour_flux*G
      end if

    else if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_VALUES) then   ! Prescribed surface temperatures

      ! Linear interpolation of input data...
      call interpolate_point_linear_1d(surface_boundary_input_times, & 
         surface_temperatures, &
         current_state%time, surface_temp, &
         extrapolate='constant') 
      
      select case(trim(units_surface_temp))
      case(degC) ! degrees C
         surface_temp = surface_temp + 273.15_DEFAULT_PRECISION
      case default ! kelvin
      end select

      if (current_state%saturated_surface)then
        current_state%surface_vapour_mixing_ratio = qsaturation(surface_temp,current_state%surface_pressure*0.01)
      else
        call interpolate_point_linear_1d(surface_boundary_input_times, & 
           surface_humidities, &
           current_state%time, current_state%surface_vapour_mixing_ratio, &
           extrapolate='constant') 
      end if

      ! Set theta_v
      current_state%theta_surf = surface_temp*&
         (current_state%surface_reference_pressure/current_state%surface_pressure)**r_over_cp
      current_state%theta_virtual_surf = current_state%theta_surf
      if (.not. current_state%passive_q) then
          current_state%theta_virtual_surf = current_state%theta_surf +  &
             current_state%global_grid%configuration%vertical%thref(2)*  &
             current_state%cq(iqv)*current_state%surface_vapour_mixing_ratio
      end if

      ! Finally set up new values of THVSURF dependent constants
      current_state%cmbc=betam*current_state%global_grid%configuration%vertical%zn(2)*G*&
         von_karman_constant/current_state%theta_virtual_surf
      current_state%rcmbc=1.0_DEFAULT_PRECISION/current_state%cmbc
      current_state%ellmocon=current_state%theta_virtual_surf/(G*von_karman_constant)

    end if

  end subroutine set_flux  

  subroutine read_configuration(current_state)
    type(model_state_type), intent(inout), target :: current_state



    integer :: ncid, time_dim
    integer :: number_input_humidities


    current_state%use_surface_boundary_conditions= & 
       options_get_logical(current_state%options_database, "use_surface_boundary_conditions")

    if (current_state%use_surface_boundary_conditions)then
    current_state%type_of_surface_boundary_conditions=options_get_integer(current_state%options_database, &
       "type_of_surface_boundary_conditions")
    current_state%use_time_varying_surface_values=options_get_logical(current_state%options_database, &
       "use_time_varying_surface_values")
  
    current_state%saturated_surface = .true. ! We will change this if we find some humidity data

    input_file=options_get_string(current_state%options_database, "surface_conditions_file")
    ! Read in the input_file
    if (trim(input_file)=='' .or. trim(input_file)=='None')then
      if (current_state%use_time_varying_surface_values)then 
        allocate(surface_boundary_input_times(MAX_SURFACE_INPUTS))
        surface_boundary_input_times=0.0
        call options_get_real_array(current_state%options_database, "surface_boundary_input_times", surface_boundary_input_times)
      end if
      if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_FLUXES)then
        allocate(surface_sensible_heat_flux(MAX_SURFACE_INPUTS),   & 
           surface_latent_heat_flux(MAX_SURFACE_INPUTS)            &
           )
        surface_sensible_heat_flux=0.0
        surface_latent_heat_flux=0.0
        call options_get_real_array(current_state%options_database, "surface_sensible_heat_flux", surface_sensible_heat_flux)
        call options_get_real_array(current_state%options_database, "surface_latent_heat_flux", surface_latent_heat_flux)
        number_input_humidities=0
      else if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_VALUES) then 
        allocate(surface_temperatures(MAX_SURFACE_INPUTS),         &
           surface_humidities(MAX_SURFACE_INPUTS)                  &
           )
        surface_temperatures=0.0
        surface_humidities=0.0
        call options_get_real_array(current_state%options_database, "surface_temperatures", surface_temperatures)
        units_surface_temp=options_get_string(current_state%options_database, "units_surface_temp")
        
        call options_get_real_array(current_state%options_database, "surface_humidities", surface_humidities)
        number_input_humidities=options_get_array_size(current_state%options_database, "surface_humidities")
      end if
    else
      call check_status(nf90_open(path = trim(input_file), mode = nf90_nowrite, ncid = ncid))
      if (log_get_logging_level() .ge. LOG_DEBUG) then
        call log_master_log(LOG_DEBUG, "Reading in surface boundary conditions from:"//trim(input_file) )
      end if
      call read_dimensions(ncid, time_dim)
      call read_variables(trim(input_file), ncid, time_dim, &
         surface_boundary_input_times, surface_temperatures, surface_humidities, &
         surface_latent_heat_flux, surface_sensible_heat_flux)
      if (.not. allocated(surface_humidities))then
        number_input_humidities=0
      else
        number_input_humidities=size(surface_humidities)
      end if
      call check_status(nf90_close(ncid))
    end if
    if (number_input_humidities>0)current_state%saturated_surface=.false.

    max_change_buoyancy_flux=options_get_real(current_state%options_database, "max_change_buoyancy_flux")
  end if

    ! Allocate lookup_tables
    allocate(current_state%lookup_table_velocity(current_state%lookup_table_entries), &
         current_state%lookup_table_ustr(current_state%lookup_table_entries))

  end subroutine read_configuration

!  BELOW COPIED FROM KIDREADER - WOULD BE GOOD TO MAKE THESE UTILITY FUNCTIONS...

  !> Reads the dimensions from the NetCDF file
  !! @param ncid The NetCDF file id
  !! @param time_dim Number of elements in the time dimension
  subroutine read_dimensions(ncid, time_dim)
    integer, intent(in) :: ncid
    integer, intent(out) ::  time_dim
    integer ::  time_dimid

    call check_status(nf90_inq_dimid(ncid, TIME_KEY, time_dimid))

    call check_status(nf90_inquire_dimension(ncid, time_dimid, len=time_dim))

  end subroutine read_dimensions

  !> Reads the variables from the NetCDF KiD model file
  !! @param ncid The id of the NetCDF file
  !! @param time_dim The number of elements in the time dimension
  !! @param time The time data field that is to be read
  !! @param surface_temperatures surface temperatures (K)
  subroutine read_variables(filename, ncid, time_dim, time, surface_temperatures, surface_humidities, &
     surface_latent_heat_flux, surface_sensible_heat_flux )
    character(*), intent(in) :: filename
    integer, intent(in) :: ncid, time_dim
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable, intent(inout) :: time
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable, intent(inout) :: surface_temperatures
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable, intent(inout) :: surface_humidities
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable, intent(inout) :: surface_latent_heat_flux
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable, intent(inout) :: surface_sensible_heat_flux

    integer :: status, variable_id

    ! Do some checking on the variable contents so that we can deal with different 
    ! variable names or missing variables
    
    ! time...
    status=nf90_inq_varid(ncid, TIME_KEY, variable_id)
    if (status==nf90_noerr)then
      allocate(time(time_dim))
      call read_single_variable(ncid, TIME_KEY, data1d=time)
    else
      call log_log(LOG_ERROR, "No recognized time variable found in"//trim(filename))
    end if
    
    status=nf90_inq_varid(ncid, SURFACE_TEMPERATURES_KEY, variable_id)
    if (status==nf90_noerr)then
      allocate(surface_temperatures(time_dim))
      call read_single_variable(ncid, SURFACE_TEMPERATURES_KEY, data1d=surface_temperatures)
    end if
    
    status=nf90_inq_varid(ncid, SURFACE_HUMIDITIES_KEY, variable_id)
    if (status==nf90_noerr)then
      allocate(surface_humidities(time_dim))
      call read_single_variable(ncid, SURFACE_HUMIDITIES_KEY, data1d=surface_humidities)
    end if
    
    status=nf90_inq_varid(ncid, SURFACE_LHF_KEY, variable_id)
    if (status==nf90_noerr)then
      allocate(surface_latent_heat_flux(time_dim))
      call read_single_variable(ncid, SURFACE_LHF_KEY, data1d=surface_latent_heat_flux)
    end if
    
    status=nf90_inq_varid(ncid, SURFACE_SHF_KEY, variable_id)
    if (status==nf90_noerr)then
      allocate(surface_sensible_heat_flux(time_dim))
      call read_single_variable(ncid, SURFACE_SHF_KEY, data1d=surface_sensible_heat_flux)
    end if

  end subroutine read_variables

  !> Reads a single variable out of a NetCDF file
  !! @param ncid The NetCDF file id
  !! @param key The variable key (name) to access
  !! @param data1d Optional one dimensional data to read into
  !! @param data3d Optional three dimensional data to read into
  subroutine read_single_variable(ncid, key, data1d, data3d)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: key
    real(kind=DEFAULT_PRECISION), dimension(:), intent(inout), optional :: data1d
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout), optional :: data3d

    integer :: variable_id
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: sdata

    call check_status(nf90_inq_varid(ncid, key, variable_id))

    if (.not. present(data1d) .and. .not. present(data3d)) return

    if (present(data1d)) then
      call check_status(nf90_get_var(ncid, variable_id, data1d))
    else
      ! 3D will reshape the data to take account of the column-row major C-F transposition
      allocate(sdata(size(data3d,1),size(data3d,3), size(data3d,2)))
      call check_status(nf90_get_var(ncid, variable_id, sdata))
      data3d(:,:,:)=reshape(sdata(:,:,:),(/size(data3d,1),size(data3d,2),size(data3d,3)/))
      deallocate(sdata)
    end if
  end subroutine read_single_variable


  !> Will check a NetCDF status and write to log_log error any decoded statuses
  !! @param status The NetCDF status flag
  subroutine check_status(status)
    integer, intent(in) :: status

    if (status /= nf90_noerr) then
      call log_log(LOG_ERROR, "NetCDF returned error code of "//trim(nf90_strerror(status)))
    end if
  end subroutine check_status

end module setfluxlook_mod
