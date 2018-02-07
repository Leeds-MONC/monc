module mcclatchey_profiles

  use datadefn_mod, only : DEFAULT_PRECISION
  use optionsdatabase_mod, only : options_get_string, options_get_integer
  use state_mod, only : model_state_type
  use netcdf, only : nf90_noerr, nf90_global, nf90_nowrite,    &
       nf90_inquire_attribute, nf90_open, nf90_strerror,       &
       nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, &
       nf90_get_var, nf90_inquire, nf90_close, nf90_get_att
  use logging_mod, only : LOG_INFO, LOG_WARN, LOG_ERROR, &
       LOG_DEBUG, log_master_log, log_log, log_get_logging_level, &
       log_master_log
  use grids_mod, only : Z_INDEX
  use def_mcc_profiles, only: str_mcc_profiles

  character(len=*), parameter ::                     &
       PRESSURE_KEY =                  "plev" ,       &  !<  NetCDF data plev key in McClatchey profile
       TEMPERATURE_PROFILE_KEY =       "t"    ,       &  !<  NetCDF data for temperature
       VAPOUR_PROFILE_KEY =            "q"    ,       &  !<  NetCDF data for vapour
       OZONE_PROFILE_KEY =             "o3"              !<  NetCDF data for ozone
  
contains
  
  subroutine read_mcclatchey_profiles(current_state, mcc)
    ! This routine is the top-level routine to read the
    ! McClatchey (mcc) data from the file. It is assumed the data is
    ! in netcdf format. In Socrates, the mcc data is in CDL format
    ! (netcdf ASCII), so it needs to be converted using ncgen
    !
    ! This routine also derives centre grid values for the McClatchey profiles
    ! It is assumed that the staggering of the grid is the same as MONC, but flipped
    ! 
    type(model_state_type), target, intent(inout) :: current_state

    type(str_mcc_profiles), intent(inout) :: mcc

    real(kind=DEFAULT_PRECISION), allocatable :: weight_upper, weight_lower
    
    real(kind=DEFAULT_PRECISION), allocatable :: temp_mcc_3d(:,:,:)

    integer, parameter :: MAX_FILE_LEN=200
    character(MAX_FILE_LEN) :: mcc_input_file

    integer :: ncid, plev_dim, k
    ! Read McClatchey profile from the global or user config. Seperate
    ! profile files required from temperature, vapour and ozone. Pressure
    ! is in all files. After each get_string, there is check status for each
    ! file.

    ! read pressure and temperature
    mcc_input_file =  &
         options_get_string(current_state%options_database, "mcc_temperature_profile")
    if (log_get_logging_level() .ge. LOG_DEBUG) then
       call log_master_log(LOG_DEBUG, &
            "Reading in McClatchey pressure and temperature from:"//trim(mcc_input_file) )
    endif
    call check_mcc_status(nf90_open(path = trim(mcc_input_file), mode = nf90_nowrite, ncid = ncid))
    call read_mcc_dimensions(ncid, plev_dim)
    ! set number of McClatchey levels
    mcc%levs = plev_dim
    ! check if temp_mcc_3d allocated and allocate if needed
    if(.not. allocated(temp_mcc_3d)) then
       allocate(temp_mcc_3d(1,1,plev_dim))
    endif
    call read_mcc_variables(trim(mcc_input_file), ncid, plev_dim, &
         TEMPERATURE_PROFILE_KEY, temp_mcc_3d, mcc%p_level)
    allocate(mcc%t_level(plev_dim))
    mcc%t_level(:)=temp_mcc_3d(1,1,:)   
    call check_mcc_status(nf90_close(ncid))
    
    ! work out the centre grid values of pressure (linear interpolation)
    allocate(mcc%p_n(plev_dim))
    do k = 2, mcc%levs
       mcc%p_n(k)=(0.5*(mcc%p_level(k-1)+(mcc%p_level(k))))
    enddo
    mcc%p_n(1)=mcc%p_level(1)+(0.5*(mcc%p_n(2)-mcc%p_level(1)))
    ! work out the centre grid values of using (log(p) weighting but linear interpolation t)
    allocate(mcc%t_n(plev_dim))
    do k = 2, mcc%levs
       weight_upper = log(mcc%p_level(k))-log(mcc%p_n(k-1))
       weight_lower = log(mcc%p_level(k))-log(mcc%p_level(k-1))
       mcc%t_n(k) = (((weight_upper*mcc%t_level(k))) + &
            (weight_lower*mcc%t_level(k-1)))/(weight_upper+weight_lower)
    enddo
    mcc%t_n(1)=mcc%t_level(1)+(0.5*(mcc%t_n(2)-mcc%t_level(1)))

    ! read mcc vapour field
    mcc_input_file =  &
         options_get_string(current_state%options_database, "mcc_vapour_profile")
    if (log_get_logging_level() .ge. LOG_DEBUG) then
       call log_master_log(LOG_DEBUG, &
            "Reading in McClatchey vapour profile from:"//trim(mcc_input_file) )
    endif
    call check_mcc_status(nf90_open(path = trim(mcc_input_file), mode = nf90_nowrite, ncid = ncid))
    call read_mcc_variables(trim(mcc_input_file), ncid, plev_dim, VAPOUR_PROFILE_KEY, temp_mcc_3d)
    allocate(mcc%q_level(plev_dim))
    mcc%q_level(:)=temp_mcc_3d(1,1,:)
    call check_mcc_status(nf90_close(ncid))

    ! work out the centre grid values of using (log(p) weighting but log interpolation q and o3)
    ! change to log of mixing ratio
    allocate(mcc%q_n(plev_dim))
    do k = 2, mcc%levs
       weight_upper = log(mcc%p_level(k))-log(mcc%p_n(k-1))
       weight_lower = log(mcc%p_level(k))-log(mcc%p_level(k-1))
       mcc%q_n(k) = (((weight_upper*log(mcc%q_level(k)))) + &
            (weight_lower*log(mcc%q_level(k-1))))/(weight_upper+weight_lower)
       mcc%q_n(k) = exp(mcc%q_n(k))
    enddo
    mcc%q_n(1)=mcc%q_level(1)+(0.5*((mcc%q_n(2)-mcc%q_level(1))))

    ! read mcc ozone field
    mcc_input_file =  &
         options_get_string(current_state%options_database, "mcc_ozone_profile")
    if (log_get_logging_level() .ge. LOG_DEBUG) then
       call log_master_log(LOG_DEBUG, &
            "Reading in McClatchey ozone profile from:"//trim(mcc_input_file) )
    endif
    call check_mcc_status(nf90_open(path = trim(mcc_input_file), mode = nf90_nowrite, ncid = ncid))
    call read_mcc_variables(trim(mcc_input_file), ncid, plev_dim, OZONE_PROFILE_KEY, temp_mcc_3d)
    allocate(mcc%o3_level(plev_dim))
    mcc%o3_level(:)=temp_mcc_3d(1,1,:)
    call check_mcc_status(nf90_close(ncid))

    ! work out the centre grid values of using (log(p) weighting but log interpolation q and o3)
    ! change to log of mixing ratio
    allocate(mcc%o3_n(plev_dim))
    do k = 2, mcc%levs
       weight_upper = log(mcc%p_level(k))-log(mcc%p_n(k-1))
       weight_lower = log(mcc%p_level(k))-log(mcc%p_level(k-1))
       mcc%o3_n(k) = (((weight_upper*log(mcc%o3_level(k)))) + &
            (weight_lower*log(mcc%o3_level(k-1))))/(weight_upper+weight_lower)
       mcc%o3_n(k) = exp(mcc%o3_n(k))
    enddo
    mcc%o3_n(1)=mcc%o3_level(1)+(0.5*(mcc%o3_n(2)-mcc%o3_level(1)))
    
    deallocate(temp_mcc_3d)
    
  end subroutine read_mcclatchey_profiles

    !> Will check for a NetCDF status of the McClatchey profile file and
  !> write to log_log error any decoded statuses
  !! @param status The NetCDF status flag
  subroutine check_mcc_status(status)
    integer, intent(in) :: status
    
    if (status /= nf90_noerr) then
       call log_log(LOG_ERROR, "NetCDF returned error code of "//trim(nf90_strerror(status)))
    end if
  end subroutine check_mcc_status

  !  BELOW COPIED FROM SETFLUXLOOK, WHICH WAS COPIED FROM KIDREADER -
  !  WOULD DEFINITELY BE GOOD TO MAKE THESE GENERIC UTILITY FUNCTIONS...
  
  !> Reads the dimensions from the NetCDF file
  !! @param ncid The NetCDF file id
  !! @param time_dim Number of elements in the time dimension
  subroutine read_mcc_dimensions(ncid, plev_dim)
    integer, intent(in) :: ncid
    integer, intent(out) ::  plev_dim
    integer ::  plev_dimid

    call check_mcc_status(nf90_inq_dimid(ncid, PRESSURE_KEY, plev_dimid))

    call check_mcc_status(nf90_inquire_dimension(ncid, plev_dimid, len=plev_dim))

  end subroutine read_mcc_dimensions

  subroutine read_mcc_variables(filename, ncid, plev_dim, MCC_PROFILE_KEY, &
       temp_mcc_3d, plev)

    character(*), intent(in) :: filename
    character(*), intent(in) :: MCC_PROFILE_KEY
    
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable, intent(inout), optional :: plev
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable, intent(inout) :: temp_mcc_3d

    integer, intent(in) :: ncid, plev_dim
    integer :: status, variable_id
    
    ! Do some checking on the variable contents so that we can deal with different 
    ! variable names or missing variables
    
    ! pressure levels, only done once with the first call...
    if (present(plev)) then
       status=nf90_inq_varid(ncid, PRESSURE_KEY, variable_id)
       if (status==nf90_noerr)then
          allocate(plev(plev_dim))
          call read_single_mcc_variable(ncid, PRESSURE_KEY, data1d=plev)
       else
          call log_log(LOG_ERROR, "No recognized pressure variable found in"//trim(filename))
       end if
    endif
    
    status=nf90_inq_varid(ncid, MCC_PROFILE_KEY, variable_id)
    if (status==nf90_noerr)then
       call read_single_mcc_variable(ncid, MCC_PROFILE_KEY, data3d=temp_mcc_3d)
    else
       call log_log(LOG_ERROR, "No recognized pressure variable found in"//trim(filename))
    end if

  end subroutine read_mcc_variables

  ! AGAIN COPIED FROM SETFLUXLOOK - REALLY NEEDS TO BE A UTILITY
  
  !> Reads a single variable out of a NetCDF file
  !! @param ncid The NetCDF file id
  !! @param key The variable key (name) to access
  !! @param data1d Optional one dimensional data to read into
  !! @param data3d Optional three dimensional data to read into
  subroutine read_single_mcc_variable(ncid, key, data1d, data3d)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: key
    real(kind=DEFAULT_PRECISION), dimension(:), intent(inout), optional :: data1d
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout), optional :: data3d

    integer :: variable_id
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: sdata

    call check_mcc_status(nf90_inq_varid(ncid, key, variable_id))

    if (.not. present(data1d) .and. .not. present(data3d)) return

    if (present(data1d)) then
      call check_mcc_status(nf90_get_var(ncid, variable_id, data1d))
    else
      ! 3D will reshape the data to take account of the column-row major C-F transposition
      allocate(sdata(size(data3d,1),size(data3d,2), size(data3d,3)))
      call check_mcc_status(nf90_get_var(ncid, variable_id, sdata))
      data3d(:,:,:)=reshape(sdata(:,:,:),(/size(data3d,1),size(data3d,2),size(data3d,3)/))
      deallocate(sdata)
    end if
  end subroutine read_single_mcc_variable

  subroutine calculate_radiation_levels(current_state, mcc)
    ! This routine calculates the total number of levels to be
    ! used in the radiation calcs (socrates). This includes the
    ! number of levels in the MONC testcase plus all the
    ! McClattchey levels which fill in the top of atmosphere.
    ! The number of levels is dependent on the McClattchey profile
    ! selected in the MONC configuration
    type(model_state_type), target, intent(inout) :: current_state

    type(str_mcc_profiles), intent(inout) :: mcc
    
    real(kind=DEFAULT_PRECISION), parameter ::       &
         fracdp = 0.75,        &
         ! fraction of to dp which the 1st rad leevl can be from top
         mindp = 2000.,        &
         ! Minimum values dp in Pa
         maxdp = 10000.
         ! Maximum dp we will allow from top
    real(kind=DEFAULT_PRECISION) :: pmindiff, pcutoff
    real(kind=DEFAULT_PRECISION), allocatable :: prefn_loc(:) !(current_state%local_grid%size(Z_INDEX))
    
    integer :: k

    allocate(prefn_loc(current_state%local_grid%size(Z_INDEX)))
    
    prefn_loc(:) = current_state%global_grid%configuration%vertical%prefn(:) 
    
    pmindiff = min(max(fracdp* &
         (prefn_loc(current_state%local_grid%size(Z_INDEX)-1) - &
         prefn_loc(current_state%local_grid%size(Z_INDEX))),    &
         mindp), maxdp)
    pcutoff = prefn_loc(current_state%local_grid%size(Z_INDEX)) - pmindiff
    
    mcc%cut = 0
    do k=1, mcc%levs
       if (mcc%p_n(k) .lt. pcutoff) mcc%cut = k
    enddo

    mcc%irad_levs = current_state%local_grid%size(Z_INDEX) + mcc%cut
    mcc%n_levels = mcc%irad_levs - 1

    deallocate(prefn_loc)
    
  end subroutine calculate_radiation_levels
  
end module mcclatchey_profiles
