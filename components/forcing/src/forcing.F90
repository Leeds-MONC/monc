!> Forcing, both subsidence and large scale
module forcing_mod
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, options_get_array_size, &
     options_get_logical_array, options_get_real_array, options_get_string_array, options_get_string
  use interpolation_mod, only: piecewise_linear_1d, piecewise_linear_2d, interpolate_point_linear_2d
  use q_indices_mod, only: get_q_index, standard_q_names
  use science_constants_mod, only: seconds_in_a_day
  use naming_conventions_mod
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, log_master_log, LOG_DEBUG, log_get_logging_level, log_log, LOG_INFO
  use conversions_mod, only : conv_to_string


  ! In order to set forcing from a netcdf file, need the following netcdf modules
  use netcdf, only : nf90_noerr, nf90_global, nf90_nowrite,    &
       nf90_inquire_attribute, nf90_open, nf90_strerror,       &
       nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, &
       nf90_get_var, nf90_inquire, nf90_close, nf90_get_att

  implicit none

#ifndef TEST_MODE
  private
#endif

  character(len=*), parameter ::                              &
       TIME_KEY =                  "time",                    &  !<  NetCDF data time key
       Z_KEY =                     "z",                       &  !<  NetCDF data height(z) key
       LEV_KEY =                   "lev",                     &  !<  NetCDF data pressure level key
       TH_KEY =                    "theta_tendency",          &  !<  NetCDF data theta tendency key
       Q_KEY =                     "q_tendency",              &  !<  NetCDF data water vapour tendency key
       WSUBS_KEY =                 "wsubs"                       !<  NetCDF data subsidence velocity key 
  
  integer, parameter :: MAX_FILE_LEN=200       !< Maximum length of surface condition input filename
  character(MAX_FILE_LEN) :: input_file
        
  integer, parameter :: DIVERGENCE=0 ! Input for subsidence forcing is a divergence profile
  integer, parameter :: SUBSIDENCE=1 ! Input for subsidence forcing is the subsidence velocity profile

  integer, parameter :: TENDENCY=0   ! Input for large-scale forcing: values are tendencies (time derivatives)
  integer, parameter :: RELAXATION=1 ! Input for large-scale forcing: values are target values to relax to over timescale
  integer, parameter :: INCREMENT=2  ! Input for large-scale forcing: values are increments (deltas) over timescale

  real(kind=DEFAULT_PRECISION), allocatable :: theta_profile(:) ! Local profile to be used in the subsidence calculation
  real(kind=DEFAULT_PRECISION), allocatable :: q_profile(:) ! Local profile to be used in the subsidence calculation
  real(kind=DEFAULT_PRECISION), allocatable :: u_profile(:) ! Local profile to be used in the subsidence calculation
  real(kind=DEFAULT_PRECISION), allocatable :: v_profile(:) ! Local profile to be used in the subsidence calculation

  real(kind=DEFAULT_PRECISION), allocatable :: dtheta_profile(:) ! Local profile to be used in time-indpendent forcing
  real(kind=DEFAULT_PRECISION), allocatable :: dq_profile(:) ! Local profile to be used in the time-indpendent forcing
  real(kind=DEFAULT_PRECISION), allocatable :: du_profile(:) ! Local profile to be used in the time-indpendent forcing
  real(kind=DEFAULT_PRECISION), allocatable :: dv_profile(:) ! Local profile to be used in the time-indpendent forcing
  ! profile_diag arrays used to store the change in field due to forcing
  real(kind=DEFAULT_PRECISION), allocatable :: du_profile_diag(:), dv_profile_diag(:), dtheta_profile_diag(:), &
       dq_profile_diag(:,:)
  ! subs_profile_diag arrays used to store the change in field due to subsidence
  real(kind=DEFAULT_PRECISION), allocatable :: du_subs_profile_diag(:), dv_subs_profile_diag(:), & 
       dtheta_subs_profile_diag(:), dq_subs_profile_diag(:,:)

  ! time dependent subsidence array (from netcdf file)
  real(kind=DEFAULT_PRECISION), allocatable :: w_subs_varies_with_time(:,:)
  real(kind=DEFAULT_PRECISION), allocatable :: forcing_input_times(:)
  
  real(kind=DEFAULT_PRECISION) :: forcing_timescale_theta ! Timescale for forcing of theta
  real(kind=DEFAULT_PRECISION) :: forcing_timescale_q     ! Timescale for forcing of q
  real(kind=DEFAULT_PRECISION) :: forcing_timescale_u     ! Timescale for forcing of u
  real(kind=DEFAULT_PRECISION) :: forcing_timescale_v     ! Timescale for forcing of v

  logical :: l_constant_forcing_theta ! Use a time-independent forcing for theta
  logical :: l_constant_forcing_q     ! Use a time-independent forcing for q
  logical :: l_constant_forcing_u     ! Use a time-independent forcing for u
  logical :: l_constant_forcing_v     ! Use a time-independent forcing for v

  integer :: constant_forcing_type_theta=TENDENCY ! Method for large-scale forcing of theta
  integer :: constant_forcing_type_q=TENDENCY     ! Method for large-scale forcing of q
  integer :: constant_forcing_type_u=RELAXATION   ! Method for large-scale forcing of u
  integer :: constant_forcing_type_v=RELAXATION   ! Method for large-scale forcing of v

  logical :: l_constant_forcing_theta_z2pressure  ! profile is a function of pressure not height

  logical :: relax_to_initial_u_profile ! For relaxation, use initial profile as the target 
  logical :: relax_to_initial_v_profile ! For relaxation, use initial profile as the target 
  logical :: relax_to_initial_theta_profile ! For relaxation, use initial profile as the target 

  logical :: use_time_varying_subsidence ! Use time dependent subsidence veocity (read from file)
  logical :: use_time_varying_theta      ! Use time dependent theta forcing (read from file)
  logical :: use_time_varying_q          ! Use time dependent water vapour forcing (read from file)

  logical :: l_subs_pl_theta ! if .true. then subsidence applied to theta field
  logical :: l_subs_pl_q     ! if .true. then subsidence applied to q fields
  
  logical :: l_subs_local_theta ! if .true. then subsidence applied locally (i.e. not with mean fields) to theta field
  logical :: l_subs_local_q     ! if .true. then subsidence applied locally (i.e. not with mean fields) to q fields

  character(len=STRING_LENGTH), dimension(:), allocatable :: names_force_pl_q  ! names of q variables to force

  ! Local tendency diagnostic variables for this component
  ! 3D tendency fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
       tend_3d_u, tend_3d_v, tend_3d_th,tend_3d_qv,                  &
       tend_3d_ql,tend_3d_qi,tend_3d_qr,tend_3d_qs,tend_3d_qg,       &
       tend_3d_tabs
  logical :: l_tend_3d_u, l_tend_3d_v, l_tend_3d_th,l_tend_3d_qv,               &
             l_tend_3d_ql,l_tend_3d_qi,l_tend_3d_qr,l_tend_3d_qs,l_tend_3d_qg,  &
             l_tend_3d_tabs
  ! Local mean tendency profile fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::                         &
       tend_pr_tot_u, tend_pr_tot_v, tend_pr_tot_th,tend_pr_tot_qv,                  &
       tend_pr_tot_ql,tend_pr_tot_qi,tend_pr_tot_qr,tend_pr_tot_qs,tend_pr_tot_qg,   &
       tend_pr_tot_tabs
  logical :: l_tend_pr_tot_u, l_tend_pr_tot_v,l_tend_pr_tot_th,l_tend_pr_tot_qv,                        &
             l_tend_pr_tot_ql,l_tend_pr_tot_qi,l_tend_pr_tot_qr,l_tend_pr_tot_qs,l_tend_pr_tot_qg,      &
             l_tend_pr_tot_tabs
  ! q indices
  integer :: iqv=0, iql=0, iqr=0, iqi=0, iqs=0, iqg=0

  integer :: diagnostic_generation_frequency

  ! Contains time varying forcing profile information
  type time_varying_forcing_profile
    real(kind=DEFAULT_PRECISION), allocatable :: forcing_times(:)    ! input forcing times
    real(kind=DEFAULT_PRECISION), allocatable :: forcing_values(:,:) ! input forcing values, interpolated to MONC heights
  end type time_varying_forcing_profile

  type(time_varying_forcing_profile), allocatable :: time_varying_subsidence, time_varying_theta, time_varying_q

  logical :: convert_input_theta_from_temperature=.false. ! If .true. input forcing data is for temperature and should
                                                          ! be converted to theta (potential temerature).

  public forcing_get_descriptor

contains

  !> Provides the component descriptor for the core to register
  !! @returns The descriptor describing this component
  type(component_descriptor_type) function forcing_get_descriptor()
    forcing_get_descriptor%name="forcing"
    forcing_get_descriptor%version=0.1
    forcing_get_descriptor%initialisation=>init_callback
    forcing_get_descriptor%timestep=>timestep_callback
    forcing_get_descriptor%finalisation=>finalisation_callback

    forcing_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    forcing_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(forcing_get_descriptor%published_fields(18+10+10))

    forcing_get_descriptor%published_fields(1)="u_subsidence"
    forcing_get_descriptor%published_fields(2)="v_subsidence"
    forcing_get_descriptor%published_fields(3)="th_subsidence"
    forcing_get_descriptor%published_fields(4)="vapour_mmr_subsidence"
    forcing_get_descriptor%published_fields(5)="cloud_mmr_subsidence"
    forcing_get_descriptor%published_fields(6)="rain_mmr_subsidence"
    forcing_get_descriptor%published_fields(7)="ice_mmr_subsidence"
    forcing_get_descriptor%published_fields(8)="snow_mmr_subsidence"
    forcing_get_descriptor%published_fields(9)="graupel_mmr_subsidence"
    forcing_get_descriptor%published_fields(10)="u_large_scale"
    forcing_get_descriptor%published_fields(11)="v_large_scale"
    forcing_get_descriptor%published_fields(12)="th_large_scale"
    forcing_get_descriptor%published_fields(13)="vapour_mmr_large_scale"
    forcing_get_descriptor%published_fields(14)="cloud_mmr_large_scale"
    forcing_get_descriptor%published_fields(15)="rain_mmr_large_scale"
    forcing_get_descriptor%published_fields(16)="ice_mmr_large_scale"
    forcing_get_descriptor%published_fields(17)="snow_mmr_large_scale"
    forcing_get_descriptor%published_fields(18)="graupel_mmr_large_scale"

    forcing_get_descriptor%published_fields(18+1)= "tend_u_forcing_3d_local"
    forcing_get_descriptor%published_fields(18+2)= "tend_v_forcing_3d_local"
    forcing_get_descriptor%published_fields(18+3)= "tend_th_forcing_3d_local"
    forcing_get_descriptor%published_fields(18+4)= "tend_qv_forcing_3d_local"
    forcing_get_descriptor%published_fields(18+5)= "tend_ql_forcing_3d_local"
    forcing_get_descriptor%published_fields(18+6)= "tend_qi_forcing_3d_local"
    forcing_get_descriptor%published_fields(18+7)= "tend_qr_forcing_3d_local"
    forcing_get_descriptor%published_fields(18+8)= "tend_qs_forcing_3d_local"
    forcing_get_descriptor%published_fields(18+9)= "tend_qg_forcing_3d_local"
    forcing_get_descriptor%published_fields(18+10)="tend_tabs_forcing_3d_local"

    forcing_get_descriptor%published_fields(18+10+1)= "tend_u_forcing_profile_total_local"
    forcing_get_descriptor%published_fields(18+10+2)= "tend_v_forcing_profile_total_local"
    forcing_get_descriptor%published_fields(18+10+3)= "tend_th_forcing_profile_total_local"
    forcing_get_descriptor%published_fields(18+10+4)= "tend_qv_forcing_profile_total_local"
    forcing_get_descriptor%published_fields(18+10+5)= "tend_ql_forcing_profile_total_local"
    forcing_get_descriptor%published_fields(18+10+6)= "tend_qi_forcing_profile_total_local"
    forcing_get_descriptor%published_fields(18+10+7)= "tend_qr_forcing_profile_total_local"
    forcing_get_descriptor%published_fields(18+10+8)= "tend_qs_forcing_profile_total_local"
    forcing_get_descriptor%published_fields(18+10+9)= "tend_qg_forcing_profile_total_local"
    forcing_get_descriptor%published_fields(18+10+10)="tend_tabs_forcing_profile_total_local"

  end function forcing_get_descriptor

  !> Field information retrieval callback, this returns information for a specific components published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve information for
  !! @param field_information Populated with information about the field
  subroutine field_information_retrieval_callback(current_state, name, field_information)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_information
    integer :: strcomp

    field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
    field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
    field_information%number_dimensions=1
    field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)

    if (name .eq. "u_subsidence") then
      field_information%enabled=current_state%u%active .and. l_subs_pl_theta .and. &
           allocated(current_state%global_grid%configuration%vertical%olzubar)
    else if (name .eq. "v_subsidence") then
      field_information%enabled=current_state%v%active .and. l_subs_pl_theta .and. &
           allocated(current_state%global_grid%configuration%vertical%olzvbar)
    else if (name .eq. "th_subsidence") then
      field_information%enabled=current_state%th%active .and. l_subs_pl_theta .and. &
           allocated(current_state%global_grid%configuration%vertical%olzthbar)
    else if (l_subs_pl_q) then
      if (name .eq. "vapour_mmr_subsidence" .or. name .eq. "cloud_mmr_subsidence" ) then
          field_information%enabled=.not. current_state%passive_q .and. & 
               current_state%number_q_fields .gt. 0               .and. &
               allocated(current_state%global_grid%configuration%vertical%olzqbar)
      else if (name .eq. "rain_mmr_subsidence" ) then
          field_information%enabled=current_state%rain_water_mixing_ratio_index .gt. 0 .and. &
               allocated(current_state%global_grid%configuration%vertical%olzqbar)   
      else if (name .eq. "ice_mmr_subsidence" ) then
          field_information%enabled=  current_state%ice_water_mixing_ratio_index .gt. 0 .and. &
               allocated(current_state%global_grid%configuration%vertical%olzqbar)     
      else if (name .eq. "snow_mmr_subsidence" ) then
          field_information%enabled=  current_state%snow_water_mixing_ratio_index .gt. 0 .and. &
               allocated(current_state%global_grid%configuration%vertical%olzqbar) 
      else if (name .eq. "graupel_mmr_subsidence" ) then
          field_information%enabled=  current_state%graupel_water_mixing_ratio_index .gt. 0 .and. &
               allocated(current_state%global_grid%configuration%vertical%olzqbar)            
      end if

    else if (name .eq. "u_large_scale") then
      field_information%enabled=current_state%u%active .and. l_constant_forcing_u
    else if (name .eq. "v_large_scale") then
      field_information%enabled=current_state%v%active .and. l_constant_forcing_v
    else if (name .eq. "th_large_scale") then
       field_information%enabled=current_state%th%active .and. l_constant_forcing_theta

    else if (l_constant_forcing_q) then 
      if (name .eq. "vapour_mmr_large_scale" .or. name .eq. "cloud_mmr_large_scale" ) then
         field_information%enabled=.not. current_state%passive_q .and. & 
              current_state%number_q_fields .gt. 0               .and. &
              allocated(current_state%global_grid%configuration%vertical%olzqbar)
      else if (name .eq. "rain_mmr_large_scale" ) then
         field_information%enabled=current_state%rain_water_mixing_ratio_index .gt. 0 .and. &
              allocated(current_state%global_grid%configuration%vertical%olzqbar)   
      else if (name .eq. "ice_mmr_large_scale" ) then
         field_information%enabled=  current_state%ice_water_mixing_ratio_index .gt. 0 .and. &
              allocated(current_state%global_grid%configuration%vertical%olzqbar)     
      else if (name .eq. "snow_mmr_large_scale" ) then
         field_information%enabled=  current_state%snow_water_mixing_ratio_index .gt. 0 .and. &
              allocated(current_state%global_grid%configuration%vertical%olzqbar) 
      else if (name .eq. "graupel_mmr_large_scale" ) then
         field_information%enabled=  current_state%graupel_water_mixing_ratio_index .gt. 0 .and. &
              allocated(current_state%global_grid%configuration%vertical%olzqbar)            
      end if   
   end if

    ! Field information for 3d
    strcomp=INDEX(name, "forcing_3d_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=3
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_u_forcing_3d_local") then
        field_information%enabled=l_tend_3d_u
      else if (name .eq. "tend_v_forcing_3d_local") then
        field_information%enabled=l_tend_3d_v
      else if (name .eq. "tend_th_forcing_3d_local") then
        field_information%enabled=l_tend_3d_th
      else if (name .eq. "tend_qv_forcing_3d_local") then
        field_information%enabled=l_tend_3d_qv
      else if (name .eq. "tend_ql_forcing_3d_local") then
        field_information%enabled=l_tend_3d_ql
      else if (name .eq. "tend_qi_forcing_3d_local") then
        field_information%enabled=l_tend_3d_qi
      else if (name .eq. "tend_qr_forcing_3d_local") then
        field_information%enabled=l_tend_3d_qr
      else if (name .eq. "tend_qs_forcing_3d_local") then
        field_information%enabled=l_tend_3d_qs
      else if (name .eq. "tend_qg_forcing_3d_local") then
        field_information%enabled=l_tend_3d_qg
      else if (name .eq. "tend_tabs_forcing_3d_local") then
        field_information%enabled=l_tend_3d_tabs
      else
        field_information%enabled=.true.
      end if

    end if !end 3d check

    ! Field information for profiles
    strcomp=INDEX(name, "forcing_profile_total_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=1
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_u_forcing_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_u
      else if (name .eq. "tend_v_forcing_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_v
      else if (name .eq. "tend_th_forcing_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_th
      else if (name .eq. "tend_qv_forcing_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qv
      else if (name .eq. "tend_ql_forcing_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_ql
      else if (name .eq. "tend_qi_forcing_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qi
      else if (name .eq. "tend_qr_forcing_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qr
      else if (name .eq. "tend_qs_forcing_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qs
      else if (name .eq. "tend_qg_forcing_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qg
      else if (name .eq. "tend_tabs_forcing_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_tabs
      else
        field_information%enabled=.true.
      end if

    end if !end profile check

  end subroutine field_information_retrieval_callback

  !> Field value retrieval callback, this returns the value of a specific published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve the value for
  !! @param field_value Populated with the value of the field
  subroutine field_value_retrieval_callback(current_state, name, field_value)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_value_type), intent(out) :: field_value

    integer :: k, n, column_size

    column_size=current_state%local_grid%size(Z_INDEX)

    ! subsidence diagnostics
    if (name .eq. "u_subsidence") then
       allocate(field_value%real_1d_array(column_size))
       field_value%real_1d_array(:)=du_subs_profile_diag(:)
    else if (name .eq. "v_subsidence") then
       allocate(field_value%real_1d_array(column_size))
       field_value%real_1d_array(:)=dv_subs_profile_diag(:)
    else if (name .eq. "th_subsidence") then
       allocate(field_value%real_1d_array(column_size))
       field_value%real_1d_array(:)= dtheta_subs_profile_diag(:)
    else if (name .eq. "vapour_mmr_subsidence") then
       allocate(field_value%real_1d_array(column_size))
       field_value%real_1d_array(:)=dq_subs_profile_diag(:,iqv)
    else if (name .eq. "cloud_mmr_subsidence") then
       allocate(field_value%real_1d_array(column_size))
       field_value%real_1d_array(:)=dq_subs_profile_diag(:,iql)
    else if (name .eq. "rain_mmr_subsidence") then
       allocate(field_value%real_1d_array(column_size))
       field_value%real_1d_array(:)=dq_subs_profile_diag(:,iqr)
    else if (name .eq. "ice_mmr_subsidence") then
       allocate(field_value%real_1d_array(column_size))
       field_value%real_1d_array(:)=dq_subs_profile_diag(:,iqi)
    else if (name .eq. "snow_mmr_subsidence") then
       allocate(field_value%real_1d_array(column_size))
       field_value%real_1d_array(:)=dq_subs_profile_diag(:,iqs)
    else if (name .eq. "graupel_mmr_subsidence") then
       allocate(field_value%real_1d_array(column_size))
       field_value%real_1d_array(:)=dq_subs_profile_diag(:,iqg)
    ! Large-scale forcing diagnostics   
    else if (name .eq. "u_large_scale") then
      allocate(field_value%real_1d_array(column_size))
      field_value%real_1d_array=get_averaged_diagnostics(current_state, du_profile_diag)
    else if (name .eq. "v_large_scale") then
      allocate(field_value%real_1d_array(column_size))
      field_value%real_1d_array=get_averaged_diagnostics(current_state, dv_profile_diag)
    else if (name .eq. "th_large_scale") then
      allocate(field_value%real_1d_array(column_size))
      field_value%real_1d_array=dtheta_profile_diag
    else if (name .eq. "vapour_mmr_large_scale") then
       allocate(field_value%real_1d_array(column_size))
       field_value%real_1d_array(:)=get_averaged_diagnostics(current_state, dq_profile_diag(:,iqv))
    else if (name .eq. "cloud_mmr_large_scale") then
       allocate(field_value%real_1d_array(column_size))
       field_value%real_1d_array(:)=get_averaged_diagnostics(current_state, dq_profile_diag(:,iql))
    else if (name .eq. "rain_mmr_large_scale") then
       allocate(field_value%real_1d_array(column_size))
       field_value%real_1d_array(:)=get_averaged_diagnostics(current_state, dq_profile_diag(:,iqr))
    else if (name .eq. "ice_mmr_large_scale") then
       allocate(field_value%real_1d_array(column_size))
       field_value%real_1d_array(:)=get_averaged_diagnostics(current_state, dq_profile_diag(:,iqi))
    else if (name .eq. "snow_mmr_large_scale") then
       allocate(field_value%real_1d_array(column_size))
       field_value%real_1d_array(:)=get_averaged_diagnostics(current_state, dq_profile_diag(:,iqs))
    else if (name .eq. "graupel_mmr_large_scale") then
       allocate(field_value%real_1d_array(column_size))
       field_value%real_1d_array(:)=get_averaged_diagnostics(current_state, dq_profile_diag(:,iqg))  

    ! 3d Tendency Fields
    else if (name .eq. "tend_u_forcing_3d_local" .and. allocated(tend_3d_u)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_u)
    else if (name .eq. "tend_v_forcing_3d_local" .and. allocated(tend_3d_v)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_v)
    else if (name .eq. "tend_th_forcing_3d_local" .and. allocated(tend_3d_th)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_th)
    else if (name .eq. "tend_qv_forcing_3d_local" .and. allocated(tend_3d_qv)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qv)
    else if (name .eq. "tend_ql_forcing_3d_local" .and. allocated(tend_3d_ql)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_ql)
    else if (name .eq. "tend_qi_forcing_3d_local" .and. allocated(tend_3d_qi)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qi)
    else if (name .eq. "tend_qr_forcing_3d_local" .and. allocated(tend_3d_qr)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qr)
    else if (name .eq. "tend_qs_forcing_3d_local" .and. allocated(tend_3d_qs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qs)
    else if (name .eq. "tend_qg_forcing_3d_local" .and. allocated(tend_3d_qg)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qg)
    else if (name .eq. "tend_tabs_forcing_3d_local" .and. allocated(tend_3d_tabs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_tabs)

    ! Profile Tendency Fields
    else if (name .eq. "tend_u_forcing_profile_total_local" .and. allocated(tend_pr_tot_u)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_u)
    else if (name .eq. "tend_v_forcing_profile_total_local" .and. allocated(tend_pr_tot_v)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_v)
    else if (name .eq. "tend_th_forcing_profile_total_local" .and. allocated(tend_pr_tot_th)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_th)
    else if (name .eq. "tend_qv_forcing_profile_total_local" .and. allocated(tend_pr_tot_qv)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qv)
    else if (name .eq. "tend_ql_forcing_profile_total_local" .and. allocated(tend_pr_tot_ql)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_ql)
    else if (name .eq. "tend_qi_forcing_profile_total_local" .and. allocated(tend_pr_tot_qi)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qi)
    else if (name .eq. "tend_qr_forcing_profile_total_local" .and. allocated(tend_pr_tot_qr)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qr)
    else if (name .eq. "tend_qs_forcing_profile_total_local" .and. allocated(tend_pr_tot_qs)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qs)
    else if (name .eq. "tend_qg_forcing_profile_total_local" .and. allocated(tend_pr_tot_qg)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qg)
    else if (name .eq. "tend_tabs_forcing_profile_total_local" .and. allocated(tend_pr_tot_tabs)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_tabs)
    end if

  end subroutine field_value_retrieval_callback  

  !> Initialises the forcing data structures
  subroutine init_callback(current_state)

    type(model_state_type), target, intent(inout) :: current_state

    integer :: nq_force ! The number of q fields apply large-scale time-independent forcing
    integer :: nzq      ! The number of input levels for subsidence/divergence profile
    integer :: i,n  ! loop counters
    integer :: iq   ! temporary q varible index

    integer :: ncid ! id for the netcdf file  
    integer :: time_dim ! number of elements in time variable, read from input file
    integer :: z_dim ! number of elements in height variable, read from input file
 
    ! Input arrays from config (always 1D) -  subsidence profile
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_subs_pl  ! subsidence node for q variables
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_subs_pl  ! subsidence node height values for q variables

    ! Input arrays from config (always 1D) - time-independent forcing
    real(kind=DEFAULT_PRECISION), dimension(:, :), allocatable :: f_force_pl_q   ! Forcing values for q variables
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_force_pl_q      ! Forcing height values for q variables
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_force_pl_theta  ! Forcing values for theta variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_force_pl_theta  ! Forcing height values for theta variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_force_pl_u      ! Forcing values for u variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_force_pl_u      ! Forcing height values for u variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_force_pl_v      ! Forcing values for v variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_force_pl_v      ! Forcing height values for v variable

    ! Read from netcdf file if used - always 2D
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: f_subs_2d  ! subsidence node for q variables
    
    integer :: subsidence_input_type=DIVERGENCE  ! Determines if we're reading in a subsidence velocity or divergence
    
    real(kind=DEFAULT_PRECISION), allocatable :: f_force_pl_q_tmp(:) !temporary 1D storage of forcing for q fields
    real(kind=DEFAULT_PRECISION), allocatable :: zgrid(:)            ! z grid to use in interpolation   

    character(len=STRING_LENGTH), dimension(:), allocatable :: units_q_force  ! units of q variable forcing
    character(len=STRING_LENGTH) :: units_theta_force='unset'  ! units of theta variable forcing
    character(len=STRING_LENGTH) :: units_u_force='unset'  ! units of theta variable forcing
    character(len=STRING_LENGTH) :: units_v_force='unset'  ! units of theta variable forcing

    integer :: k
    logical :: l_qdiag

    allocate(u_profile(current_state%local_grid%size(Z_INDEX)),      &
       v_profile(current_state%local_grid%size(Z_INDEX)),            &
       theta_profile(current_state%local_grid%size(Z_INDEX)), &
       q_profile(current_state%local_grid%size(Z_INDEX)))

    allocate(dtheta_profile(current_state%local_grid%size(Z_INDEX)), &
       dq_profile(current_state%local_grid%size(Z_INDEX)),           &
       du_profile(current_state%local_grid%size(Z_INDEX)),           &
       dv_profile(current_state%local_grid%size(Z_INDEX)))

    allocate(du_profile_diag(current_state%local_grid%size(Z_INDEX)), &
         dv_profile_diag(current_state%local_grid%size(Z_INDEX)),     &
         dtheta_profile_diag(current_state%local_grid%size(Z_INDEX)), &
         dq_profile_diag(current_state%local_grid%size(Z_INDEX), current_state%number_q_fields))
    
    allocate(du_subs_profile_diag(current_state%local_grid%size(Z_INDEX)), &
         dv_subs_profile_diag(current_state%local_grid%size(Z_INDEX)),     &
         dtheta_subs_profile_diag(current_state%local_grid%size(Z_INDEX)), &
         dq_subs_profile_diag(current_state%local_grid%size(Z_INDEX), current_state%number_q_fields))

    allocate(zgrid(current_state%local_grid%size(Z_INDEX)))

    ! assign microphysics indexes, needed for the diagnostic output
    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then 
       iqv=get_q_index(standard_q_names%VAPOUR, 'forcing')                         
       iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'forcing')
    endif
    if (current_state%rain_water_mixing_ratio_index > 0) &
         iqr = current_state%rain_water_mixing_ratio_index
    if (current_state%ice_water_mixing_ratio_index > 0) &
         iqi = current_state%ice_water_mixing_ratio_index 
    if (current_state%snow_water_mixing_ratio_index > 0) &
         iqs = current_state%snow_water_mixing_ratio_index
    if (current_state%graupel_water_mixing_ratio_index > 0) &
         iqg = current_state%graupel_water_mixing_ratio_index

    ! time_varying forcing initialization
    use_time_varying_subsidence= &
         options_get_logical(current_state%options_database, "use_time_varying_subsidence")
    if ((l_subs_pl_theta .or. l_subs_pl_q) .and. use_time_varying_subsidence) then
        allocate(time_varying_subsidence)
        call init_time_varying_forcing(current_state, time_varying_subsidence, WSUBS_KEY,                  &
             options_get_string(current_state%options_database, "varying_subsidence_file"), &
             options_get_string(current_state%options_database, "varying_subsidence_coordinate"))
    end if

    use_time_varying_theta= &
         options_get_logical(current_state%options_database, "use_time_varying_theta")
    if (use_time_varying_theta) then
        convert_input_theta_from_temperature=options_get_logical(current_state%options_database, &
                                                                       "convert_input_theta_from_temperature")
        allocate(time_varying_theta)
        call init_time_varying_forcing(current_state, time_varying_theta, TH_KEY,                     &
             options_get_string(current_state%options_database, "varying_theta_file"), &
             options_get_string(current_state%options_database, "varying_theta_coordinate"))
    end if

    use_time_varying_q= &
         options_get_logical(current_state%options_database, "use_time_varying_q")
    if (use_time_varying_q) then
        allocate(time_varying_q)
        call init_time_varying_forcing(current_state, time_varying_q, Q_KEY,                      &
             options_get_string(current_state%options_database, "varying_q_file"), &
             options_get_string(current_state%options_database, "varying_q_coordinate"))
    end if

    ! Subsidence forcing initialization
    
    l_subs_pl_theta=options_get_logical(current_state%options_database, "l_subs_pl_theta")
    l_subs_pl_q=options_get_logical(current_state%options_database, "l_subs_pl_q")
    subsidence_input_type=options_get_integer(current_state%options_database, "subsidence_input_type")
    l_subs_local_theta=options_get_logical(current_state%options_database, "subsidence_local_theta")
    l_subs_local_q=options_get_logical(current_state%options_database, "subsidence_local_q")

    if ((l_subs_pl_theta .and. .not. l_subs_local_theta) .or. &
       (l_subs_pl_q .and. .not. l_subs_local_q))then
      if (.not. is_component_enabled(current_state%options_database, "mean_profiles")) then
        call log_master_log(LOG_ERROR, "subsidence requires the mean profiles component to be enabled")
      end if
    end if

    if ((l_subs_pl_theta .or. l_subs_pl_q) .and. .not. use_time_varying_subsidence) then
      allocate(z_subs_pl(options_get_array_size(current_state%options_database, "z_subs_pl")), &
           f_subs_pl(options_get_array_size(current_state%options_database, "f_subs_pl")))
      call options_get_real_array(current_state%options_database, "z_subs_pl", z_subs_pl)
      call options_get_real_array(current_state%options_database, "f_subs_pl", f_subs_pl)      
      ! Get profiles
      zgrid=current_state%global_grid%configuration%vertical%z(:)
      call piecewise_linear_1d(z_subs_pl(1:size(z_subs_pl)), f_subs_pl(1:size(f_subs_pl)), zgrid, &
           current_state%global_grid%configuration%vertical%w_subs)
      if (subsidence_input_type==DIVERGENCE) then
        current_state%global_grid%configuration%vertical%w_subs(:) = &
            -1.0*current_state%global_grid%configuration%vertical%w_subs(:)*zgrid(:)
      end if
      deallocate(z_subs_pl, f_subs_pl)  
    end if

   ! Time independent large-scale forcing (proxy for e.g. advection/radiation)
   ! This probably isn't the right place to be doing this
   if (.not. allocated(current_state%l_forceq))then
      allocate(current_state%l_forceq(current_state%number_q_fields))
      current_state%l_forceq=.false.
   end if

    l_constant_forcing_theta=options_get_logical(current_state%options_database, "l_constant_forcing_theta")
    l_constant_forcing_q=options_get_logical(current_state%options_database, "l_constant_forcing_q")
    l_constant_forcing_u=options_get_logical(current_state%options_database, "l_constant_forcing_u")
    l_constant_forcing_v=options_get_logical(current_state%options_database, "l_constant_forcing_v")

    if (l_constant_forcing_q) then
      allocate(names_force_pl_q(options_get_array_size(current_state%options_database, "names_constant_forcing_q")))
      call options_get_string_array(current_state%options_database, "names_constant_forcing_q", names_force_pl_q)
    end if
    
    if (l_constant_forcing_theta)then
      constant_forcing_type_theta=options_get_integer(current_state%options_database, "constant_forcing_type_theta")
      forcing_timescale_theta=options_get_real(current_state%options_database, "forcing_timescale_theta")
      l_constant_forcing_theta_z2pressure=options_get_logical(current_state%options_database,"l_constant_forcing_theta_z2pressure")

      allocate(z_force_pl_theta(options_get_array_size(current_state%options_database, "z_force_pl_theta")), &
           f_force_pl_theta(options_get_array_size(current_state%options_database, "f_force_pl_theta")))
      call options_get_real_array(current_state%options_database, "z_force_pl_theta", z_force_pl_theta)
      call options_get_real_array(current_state%options_database, "f_force_pl_theta", f_force_pl_theta)
      ! Get profiles
      relax_to_initial_theta_profile=options_get_logical(current_state%options_database, "relax_to_initial_theta_profile")
      if (relax_to_initial_theta_profile)then
        current_state%global_grid%configuration%vertical%theta_force(:) = &
           current_state%global_grid%configuration%vertical%theta_init(:)
      else
        if (l_constant_forcing_theta_z2pressure)then
          zgrid=current_state%global_grid%configuration%vertical%zn(:)
        else
          zgrid=current_state%global_grid%configuration%vertical%prefn(:)
        end if
        call piecewise_linear_1d(z_force_pl_theta(1:size(z_force_pl_theta)), f_force_pl_theta(1:size(f_force_pl_theta)), zgrid, &
           current_state%global_grid%configuration%vertical%theta_force)
      end if
      
      ! Unit conversions...
      convert_input_theta_from_temperature=options_get_logical(current_state%options_database, &
                                                                       "convert_input_theta_from_temperature")
      if (convert_input_theta_from_temperature)then ! Input is temperature not theta
        current_state%global_grid%configuration%vertical%theta_force(:) =   &
           current_state%global_grid%configuration%vertical%theta_force(:)* &
           current_state%global_grid%configuration%vertical%prefrcp(:)
      end if
      
      if (constant_forcing_type_theta==TENDENCY)then
        units_theta_force=options_get_string(current_state%options_database, "units_theta_force")
        select case(trim(units_theta_force))
        case(k_per_day)
          current_state%global_grid%configuration%vertical%theta_force(:) = &
             current_state%global_grid%configuration%vertical%theta_force(:)/seconds_in_a_day
        case default !(k_per_second)
        end select
      end if
      deallocate(z_force_pl_theta, f_force_pl_theta)
   end if
                                                           
#ifdef U_ACTIVE
    if (l_constant_forcing_u)then
      constant_forcing_type_u=options_get_integer(current_state%options_database, "constant_forcing_type_u")
      forcing_timescale_u=options_get_real(current_state%options_database, "forcing_timescale_u")
      relax_to_initial_u_profile=options_get_logical(current_state%options_database, "relax_to_initial_u_profile")
      if (relax_to_initial_u_profile)then
        current_state%global_grid%configuration%vertical%u_force(:) = &
           current_state%global_grid%configuration%vertical%u_init(:)
      else
        allocate(z_force_pl_u(options_get_array_size(current_state%options_database, "z_force_pl_u")), &
           f_force_pl_u(options_get_array_size(current_state%options_database, "f_force_pl_u")))
        call options_get_real_array(current_state%options_database, "z_force_pl_u", z_force_pl_u)
        call options_get_real_array(current_state%options_database, "f_force_pl_u", f_force_pl_u)
        ! Get profiles
        zgrid=current_state%global_grid%configuration%vertical%zn(:)
        call piecewise_linear_1d(z_force_pl_u(1:size(z_force_pl_u)), f_force_pl_u(1:size(f_force_pl_u)), zgrid, &
           current_state%global_grid%configuration%vertical%u_force)
        deallocate(z_force_pl_u, f_force_pl_u)
      end if


      if (constant_forcing_type_u==TENDENCY)then
        ! Unit conversions...
        units_u_force=options_get_string(current_state%options_database, "units_u_force")
        select case(trim(units_u_force))
        case(m_per_second_per_day)
          current_state%global_grid%configuration%vertical%u_force(:) = &
             current_state%global_grid%configuration%vertical%u_force(:)/seconds_in_a_day
        case default  !(m_per_second_per_second)
        end select
      end if
    end if
#endif
                                                           
#ifdef V_ACTIVE
    if (l_constant_forcing_v)then
      constant_forcing_type_v=options_get_integer(current_state%options_database, "constant_forcing_type_v")
      forcing_timescale_v=options_get_real(current_state%options_database, "forcing_timescale_v")
      relax_to_initial_v_profile=options_get_logical(current_state%options_database, "relax_to_initial_v_profile")
      if (relax_to_initial_v_profile)then
        current_state%global_grid%configuration%vertical%v_force(:) = &
           current_state%global_grid%configuration%vertical%v_init(:)
      else
        allocate(z_force_pl_v(options_get_array_size(current_state%options_database, "z_force_pl_v")), &
           f_force_pl_v(options_get_array_size(current_state%options_database, "f_force_pl_v")))
        call options_get_real_array(current_state%options_database, "z_force_pl_v", z_force_pl_v)
        call options_get_real_array(current_state%options_database, "f_force_pl_v", f_force_pl_v)
        ! Get profiles
        zgrid=current_state%global_grid%configuration%vertical%zn(:)
        call piecewise_linear_1d(z_force_pl_v(1:size(z_force_pl_v)), f_force_pl_v(1:size(f_force_pl_v)), zgrid, &
           current_state%global_grid%configuration%vertical%v_force)
        deallocate(z_force_pl_v, f_force_pl_v)
      end if


      if (constant_forcing_type_v==TENDENCY)then
        ! Unit conversions...
        units_v_force=options_get_string(current_state%options_database, "units_v_force")
        select case(trim(units_v_force))
        case(m_per_second_per_day)
          current_state%global_grid%configuration%vertical%v_force(:) = &
             current_state%global_grid%configuration%vertical%v_force(:)/seconds_in_a_day
        case default !(m_per_second_per_second)
        end select
      end if
    end if 
#endif   
        
    if (l_constant_forcing_q) then
      constant_forcing_type_q=options_get_integer(current_state%options_database, "constant_forcing_type_q")
      forcing_timescale_q=options_get_real(current_state%options_database, "forcing_timescale_q")
      nq_force=size(names_force_pl_q)
      allocate(z_force_pl_q(options_get_array_size(current_state%options_database, "z_force_pl_q")))
      call options_get_real_array(current_state%options_database, "z_force_pl_q", z_force_pl_q)
      nzq=size(z_force_pl_q)
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      allocate(f_force_pl_q_tmp(nq_force*nzq))
      call options_get_real_array(current_state%options_database, "f_force_pl_q", f_force_pl_q_tmp)
      allocate(f_force_pl_q(nzq, nq_force))
      f_force_pl_q(1:nzq, 1:nq_force)=reshape(f_force_pl_q_tmp, (/nzq, nq_force/))

      allocate(units_q_force(options_get_array_size(current_state%options_database, "units_q_force")))
      call options_get_string_array(current_state%options_database, "units_q_force", units_q_force)
      do n=1, nq_force
        iq=get_q_index(trim(names_force_pl_q(n)), 'forcing:time-independent')
        call piecewise_linear_1d(z_force_pl_q(1:nzq), f_force_pl_q(1:nzq,n), zgrid, &
           current_state%global_grid%configuration%vertical%q_force(:,iq))
        
        current_state%l_forceq(iq)=.true.

        ! Unit conversions...
        if (constant_forcing_type_q==TENDENCY)then
          select case(trim(units_q_force(n)))
          case(kg_per_kg_per_day)
            current_state%global_grid%configuration%vertical%q_force(:,iq) = &
               current_state%global_grid%configuration%vertical%q_force(:,iq)/seconds_in_a_day
          case(g_per_kg_per_day)
            current_state%global_grid%configuration%vertical%q_force(:,iq) = &
               0.001*current_state%global_grid%configuration%vertical%q_force(:,iq)/seconds_in_a_day
          case(g_per_kg_per_second)
            current_state%global_grid%configuration%vertical%q_force(:,iq) = &
               0.001*current_state%global_grid%configuration%vertical%q_force(:,iq)
          case default !(kg_per_kg_per_second)
          end select
        end if
      end do
      deallocate(f_force_pl_q_tmp, units_q_force, f_force_pl_q, z_force_pl_q)  
    end if

    deallocate(zgrid)

    ! Set tendency diagnostic logicals based on availability
    ! Need to use 3d tendencies to compute the profiles, so they will be allocated
    !      in the case where profiles are available
    l_qdiag =  (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0)

    l_tend_pr_tot_u   = current_state%u%active
    l_tend_pr_tot_v   = current_state%v%active
    l_tend_pr_tot_th  = current_state%th%active
    l_tend_pr_tot_qv  = l_qdiag .and. current_state%number_q_fields .ge. 1
    l_tend_pr_tot_ql  = l_qdiag .and. current_state%number_q_fields .ge. 2
    l_tend_pr_tot_qi  = l_qdiag .and. current_state%number_q_fields .ge. 11
    l_tend_pr_tot_qr  = l_qdiag .and. current_state%number_q_fields .ge. 11
    l_tend_pr_tot_qs  = l_qdiag .and. current_state%number_q_fields .ge. 11
    l_tend_pr_tot_qg  = l_qdiag .and. current_state%number_q_fields .ge. 11
    l_tend_pr_tot_tabs  = l_tend_pr_tot_th

    l_tend_3d_u   = current_state%u%active .or. l_tend_pr_tot_u
    l_tend_3d_v   = current_state%v%active .or. l_tend_pr_tot_v
    l_tend_3d_th  = current_state%th%active .or. l_tend_pr_tot_th
    l_tend_3d_qv  = (l_qdiag .and. current_state%number_q_fields .ge. 1)  .or. l_tend_pr_tot_qv
    l_tend_3d_ql  = (l_qdiag .and. current_state%number_q_fields .ge. 2)  .or. l_tend_pr_tot_ql
    l_tend_3d_qi  = (l_qdiag .and. current_state%number_q_fields .ge. 11) .or. l_tend_pr_tot_qi
    l_tend_3d_qr  = (l_qdiag .and. current_state%number_q_fields .ge. 11) .or. l_tend_pr_tot_qr
    l_tend_3d_qs  = (l_qdiag .and. current_state%number_q_fields .ge. 11) .or. l_tend_pr_tot_qs
    l_tend_3d_qg  = (l_qdiag .and. current_state%number_q_fields .ge. 11) .or. l_tend_pr_tot_qg
    l_tend_3d_tabs = l_tend_3d_th

    ! Allocate 3d tendency fields upon availability
    if (l_tend_3d_u) then
      allocate( tend_3d_u(current_state%local_grid%size(Z_INDEX),  &
                          current_state%local_grid%size(Y_INDEX),  &
                          current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_v) then
      allocate( tend_3d_v(current_state%local_grid%size(Z_INDEX),  &
                          current_state%local_grid%size(Y_INDEX),  &
                          current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_th) then
      allocate( tend_3d_th(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qv) then
      allocate( tend_3d_qv(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_ql) then
      allocate( tend_3d_ql(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qi) then
      allocate( tend_3d_qi(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qr) then
      allocate( tend_3d_qr(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qs) then
      allocate( tend_3d_qs(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qg) then
      allocate( tend_3d_qg(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_tabs) then
      allocate( tend_3d_tabs(current_state%local_grid%size(Z_INDEX),  &
                             current_state%local_grid%size(Y_INDEX),  &
                             current_state%local_grid%size(X_INDEX)   )    )
    endif

    ! Allocate profile tendency fields upon availability
    if (l_tend_pr_tot_u) then
      allocate( tend_pr_tot_u(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_v) then
      allocate( tend_pr_tot_v(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_th) then
      allocate( tend_pr_tot_th(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_qv) then
      allocate( tend_pr_tot_qv(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_ql) then
      allocate( tend_pr_tot_ql(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_qi) then
      allocate( tend_pr_tot_qi(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_qr) then
      allocate( tend_pr_tot_qr(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_qs) then
      allocate( tend_pr_tot_qs(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_qg) then
      allocate( tend_pr_tot_qg(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_tabs) then
      allocate( tend_pr_tot_tabs(current_state%local_grid%size(Z_INDEX)) )
    endif

    ! Save the sampling_frequency to force diagnostic calculation on select time steps
    diagnostic_generation_frequency=options_get_integer(current_state%options_database, "sampling_frequency")


  end subroutine init_callBack

  !> Called for each data column and will determine the forcing values in x and y which are then applied to the field
  !! source terms
  !! @param current_state The current model state_mod
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: current_x_index, current_y_index, target_x_index, target_y_index, k
    logical :: calculate_diagnostics
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: temp_prof

    calculate_diagnostics = ((current_state%time_basis .and. current_state%timestep == current_state%sample_timestep) .or.      &
                      (.not. current_state%time_basis .and. mod(current_state%timestep, diagnostic_generation_frequency) == 0))

    current_x_index=current_state%column_local_x
    current_y_index=current_state%column_local_y
    target_y_index=current_y_index-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=current_x_index-current_state%local_grid%halo_size(X_INDEX)

    if (current_state%first_timestep_column) then
      du_profile_diag=0.0_DEFAULT_PRECISION
      dv_profile_diag=0.0_DEFAULT_PRECISION
      dtheta_profile_diag=0.0_DEFAULT_PRECISION
      dq_profile_diag=0.0_DEFAULT_PRECISION
      du_subs_profile_diag=0.0_DEFAULT_PRECISION
      dv_subs_profile_diag=0.0_DEFAULT_PRECISION
      dtheta_subs_profile_diag=0.0_DEFAULT_PRECISION
      dq_subs_profile_diag=0.0_DEFAULT_PRECISION
      if (l_tend_pr_tot_u) then
        tend_pr_tot_u(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_v) then
        tend_pr_tot_v(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_th) then
        tend_pr_tot_th(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qv) then
        tend_pr_tot_qv(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_ql) then
        tend_pr_tot_ql(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qi) then
        tend_pr_tot_qi(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qr) then
        tend_pr_tot_qr(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qs) then
        tend_pr_tot_qs(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qg) then
        tend_pr_tot_qg(:)=0.0_DEFAULT_PRECISION
      end if
      if (l_tend_pr_tot_tabs) then
        tend_pr_tot_tabs(:)=0.0_DEFAULT_PRECISION
      endif
    end if    

    if (current_state%halo_column .or. current_state%timestep<3) return

    if (calculate_diagnostics) &
        call save_precomponent_tendencies(current_state, current_x_index, current_y_index, target_x_index, target_y_index)

    ! AH: perform subsidence calculation but first determine if time varying or constant
    !     If timevarying then work out the profile of subsidence for the given time and 
    !     assign to w_subs, which is used in apply_subsidence_to...
    !     
    if (use_time_varying_subsidence) then 
       call interpolate_point_linear_2d(time_varying_subsidence%forcing_times,         &
            time_varying_subsidence%forcing_values,                                    &
            current_state%time, current_state%global_grid%configuration%vertical%w_subs)
    endif ! if not w_subs is constant and set in the init_callback


    ! Apply time-varying theta and q (vapour only).
    ! This functionality permits the user to apply a constant forcing separately as long as the 
    !   theta forcing is consistently in theta or absolute temperature units in both cases because
    !   they share the convert_input_theta_from_temperature logical.
    if (use_time_varying_theta) then
       ! Obtain the profile, interpolated to the current time
       call interpolate_point_linear_2d(time_varying_theta%forcing_times,         &
            time_varying_theta%forcing_values,                                    &
            current_state%time, temp_prof)

      ! Unit conversions
      if (convert_input_theta_from_temperature)then ! Input is temperature not theta
        temp_prof = temp_prof * current_state%global_grid%configuration%vertical%prefrcp(:)  
      end if

      ! Record the diagnostic and apply the forcing
      dtheta_profile_diag = dtheta_profile_diag + temp_prof
      do k=2,current_state%local_grid%size(Z_INDEX)-1
        current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x) = &
             current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x) &
             + temp_prof(k)
      end do
    endif ! use_time_varying_theta

    if (use_time_varying_q) then
       ! Obtain the profile, interpolated to the current time
       call interpolate_point_linear_2d(time_varying_q%forcing_times,         &
            time_varying_q%forcing_values,                                    &
            current_state%time, temp_prof)

      ! Record the diagnostic and apply the forcing
      dq_profile_diag(:,iqv) = dq_profile_diag(:,iqv) + temp_prof
      do k=2,current_state%local_grid%size(Z_INDEX)-1
        current_state%sq(iqv)%data(k,current_state%column_local_y,current_state%column_local_x) = &
             current_state%sq(iqv)%data(k,current_state%column_local_y,current_state%column_local_x) &
             + temp_prof(k)
      end do
    endif ! use_time_varying_q


    if (l_subs_pl_theta) then
      call apply_subsidence_to_flow_fields(current_state)
      call apply_subsidence_to_theta(current_state)
    end if
    if (l_subs_pl_q) call apply_subsidence_to_q_fields(current_state)

    if (l_constant_forcing_theta) call apply_time_independent_forcing_to_theta(current_state)                                                         
#ifdef U_ACTIVE
    if (l_constant_forcing_u) call apply_time_independent_forcing_to_u(current_state)
#endif                                                       
#ifdef V_ACTIVE
    if (l_constant_forcing_v) call apply_time_independent_forcing_to_v(current_state)
#endif
    if (l_constant_forcing_q) call apply_time_independent_forcing_to_q(current_state)

    if (calculate_diagnostics) &
        call compute_component_tendencies(current_state, current_x_index, current_y_index, target_x_index, target_y_index)

  end subroutine timestep_callback

  subroutine apply_subsidence_to_flow_fields(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: k
    real(kind=DEFAULT_PRECISION) :: usub, vsub    


    if (l_subs_local_theta)then ! Use local gradients not global means
      u_profile(:)=current_state%zu%data(:,current_state%column_local_y,current_state%column_local_x)
      v_profile(:)=current_state%zv%data(:,current_state%column_local_y,current_state%column_local_x)
    else
      u_profile(:)=current_state%global_grid%configuration%vertical%olzubar(:)
      v_profile(:)=current_state%global_grid%configuration%vertical%olzvbar(:)
    end if

    do k=2,current_state%local_grid%size(Z_INDEX)-1                                                             
#ifdef U_ACTIVE
      usub =  2.0 * (current_state%global_grid%configuration%vertical%w_subs(k-1)*              &
         current_state%global_grid%configuration%vertical%tzc1(k)*(u_profile(k)-u_profile(k-1)) &
           + current_state%global_grid%configuration%vertical%w_subs(k)*                        &
           current_state%global_grid%configuration%vertical%tzc2(k)*                            &
           (u_profile(k+1) - u_profile(k)))
      current_state%su%data(k,current_state%column_local_y,current_state%column_local_x) =      &
           current_state%su%data(k,current_state%column_local_y,current_state%column_local_x) - usub
      du_subs_profile_diag(k) = du_subs_profile_diag(k) - usub
#endif
#ifdef V_ACTIVE     
      vsub =  2.0 * (current_state%global_grid%configuration%vertical%w_subs(k-1)*              &
         current_state%global_grid%configuration%vertical%tzc1(k)*(v_profile(k)-v_profile(k-1)) &
         + current_state%global_grid%configuration%vertical%w_subs(k)*                          &
         current_state%global_grid%configuration%vertical%tzc2(k)*                              &
         (v_profile(k+1) - v_profile(k)))
      current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x) =      &
           current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x) - vsub
      dv_subs_profile_diag(k) = dv_subs_profile_diag(k) - vsub 
#endif
    end do
    k=current_state%local_grid%size(Z_INDEX)
#ifdef U_ACTIVE
    usub =  2.0 * (current_state%global_grid%configuration%vertical%w_subs(k-1)*                &
         current_state%global_grid%configuration%vertical%tzc1(k)*(u_profile(k)-u_profile(k-1)))
    current_state%su%data(k,current_state%column_local_y,current_state%column_local_x) =        &
         current_state%su%data(k,current_state%column_local_y,current_state%column_local_x) - usub
    du_subs_profile_diag(k) = du_subs_profile_diag(k) - usub
#endif
#ifdef V_ACTIVE
    vsub =  2.0 * (current_state%global_grid%configuration%vertical%w_subs(k-1)*                &
         current_state%global_grid%configuration%vertical%tzc1(k)*(v_profile(k)-v_profile(k-1)))
    current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x) =        &
         current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x) - vsub
    dv_subs_profile_diag(k) = dv_subs_profile_diag(k) - vsub
#endif
  end subroutine apply_subsidence_to_flow_fields

  subroutine apply_subsidence_to_theta(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: k
    real(kind=DEFAULT_PRECISION) :: thsub

    if (l_subs_local_theta)then ! Use local gradients not global means
      theta_profile(:)=current_state%zth%data(:,current_state%column_local_y,current_state%column_local_x) &
         + current_state%global_grid%configuration%vertical%thref(:)
    else
      theta_profile(:)=current_state%global_grid%configuration%vertical%olzthbar(:) &
         + current_state%global_grid%configuration%vertical%thref(:)
    end if

    do k=2,current_state%local_grid%size(Z_INDEX)-1
      thsub = current_state%global_grid%configuration%vertical%w_subs(k)*   &
         (theta_profile(k+1) - theta_profile(k))*                           &
         current_state%global_grid%configuration%vertical%rdzn(k+1)
      
      current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x) = &
           current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x) - thsub
      dtheta_subs_profile_diag(k) = dtheta_subs_profile_diag(k) - thsub
    end do
    k=current_state%local_grid%size(Z_INDEX)
    thsub = current_state%global_grid%configuration%vertical%w_subs(k)* &
       (theta_profile(k) - theta_profile(k-1))*                         &
         current_state%global_grid%configuration%vertical%rdzn(k)

    current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x) = &
         current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x) - thsub
    dtheta_subs_profile_diag(k) = dtheta_subs_profile_diag(k) - thsub
  end subroutine apply_subsidence_to_theta

  subroutine apply_subsidence_to_q_fields(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: k, n
    real(kind=DEFAULT_PRECISION) :: qsub


    do n=1,current_state%number_q_fields
      if (l_subs_local_q)then ! Use local gradients not global means
        q_profile(:)=current_state%zq(n)%data(:,current_state%column_local_y,current_state%column_local_x)
      else
        q_profile(:)=current_state%global_grid%configuration%vertical%olzqbar(:,n)
      end if
      do k=2,current_state%local_grid%size(Z_INDEX)-1
        qsub = current_state%global_grid%configuration%vertical%w_subs(k)*  &
           (q_profile(k+1) - q_profile(k))*                               &
           current_state%global_grid%configuration%vertical%rdzn(k+1)
        current_state%sq(n)%data(k,current_state%column_local_y,current_state%column_local_x) = &
           current_state%sq(n)%data(k,current_state%column_local_y,current_state%column_local_x) - qsub
        dq_subs_profile_diag(k,n) = dq_subs_profile_diag(k,n) - qsub
      end do
      k=current_state%local_grid%size(Z_INDEX)
      qsub = current_state%global_grid%configuration%vertical%w_subs(k)*    &
         (q_profile(k) - q_profile(k-1))*                                   &
         current_state%global_grid%configuration%vertical%rdzn(k)
      current_state%sq(n)%data(k,current_state%column_local_y,current_state%column_local_x) = &
         current_state%sq(n)%data(k,current_state%column_local_y,current_state%column_local_x) - qsub
      dq_subs_profile_diag(k,n) = dq_subs_profile_diag(k,n) - qsub
    end do
  end subroutine apply_subsidence_to_q_fields

  subroutine apply_time_independent_forcing_to_theta(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: k
    real(kind=DEFAULT_PRECISION) :: dtm_scale

    if (constant_forcing_type_theta==TENDENCY)then
      dtm_scale=1.0_DEFAULT_PRECISION
    else !  constant_forcing_type_theta==(RELAXATION or INCREMENT)
      dtm_scale=1.0_DEFAULT_PRECISION/forcing_timescale_theta
    end if
    
    if (constant_forcing_type_theta==RELAXATION)then
      dtheta_profile(:)=dtm_scale * (current_state%global_grid%configuration%vertical%theta_force(:) - &
         current_state%zth%data(:,current_state%column_local_y,current_state%column_local_x) -       &
         current_state%global_grid%configuration%vertical%thref(:))
    else !  constant_forcing_type_theta==(TENDENCY or INCREMENT)
      dtheta_profile(:)=dtm_scale * current_state%global_grid%configuration%vertical%theta_force(:)
   end if


    dtheta_profile_diag=dtheta_profile_diag+(dtheta_profile)

    do k=2,current_state%local_grid%size(Z_INDEX)-1
      current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x) = &
           current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x) &
           + dtheta_profile(k)
    end do 

  end subroutine apply_time_independent_forcing_to_theta

  subroutine apply_time_independent_forcing_to_q(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: n, k
    real(kind=DEFAULT_PRECISION) :: dtm_scale

    do n=1,current_state%number_q_fields
      if (current_state%l_forceq(n))then
        if (constant_forcing_type_q==TENDENCY)then
          dtm_scale=1.0_DEFAULT_PRECISION
        else !  constant_forcing_type_q==(RELAXATION or INCREMENT)
          dtm_scale=1.0_DEFAULT_PRECISION/forcing_timescale_q
        end if
        
        if (constant_forcing_type_q==RELAXATION)then
          dq_profile(:)=dtm_scale * (current_state%global_grid%configuration%vertical%q_force(:,n) - &
             current_state%zq(n)%data(:,current_state%column_local_y,current_state%column_local_x))
        else !  constant_forcing_type_q==(TENDENCY or INCREMENT)
          dq_profile(:)=dtm_scale * current_state%global_grid%configuration%vertical%q_force(:,n)
        end if

        dq_profile_diag(:,n)=dq_profile_diag(:,n)+dq_profile

        do k=2,current_state%local_grid%size(Z_INDEX)-1
          current_state%sq(n)%data(k,current_state%column_local_y,current_state%column_local_x) = &
             current_state%sq(n)%data(k,current_state%column_local_y,current_state%column_local_x) &
             + dq_profile(k)
        end do
      end if
    end do

  end subroutine apply_time_independent_forcing_to_q

  subroutine apply_time_independent_forcing_to_u(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: k
    real(kind=DEFAULT_PRECISION) :: dtm_scale

    if (constant_forcing_type_u==TENDENCY)then
      dtm_scale=1.0_DEFAULT_PRECISION
    else !  constant_forcing_type_u==(RELAXATION or INCREMENT)
      dtm_scale=1.0_DEFAULT_PRECISION/forcing_timescale_u
    end if
    
    if (constant_forcing_type_u==RELAXATION)then
      du_profile(:)=dtm_scale * (current_state%global_grid%configuration%vertical%u_force(:) - &
           current_state%global_grid%configuration%vertical%olzubar(:))
    else !  constant_forcing_type_u==(TENDENCY or INCREMENT)
      du_profile(:)=dtm_scale * current_state%global_grid%configuration%vertical%u_force(:)
    end if

    du_profile_diag=du_profile_diag+du_profile

    do k=2,current_state%local_grid%size(Z_INDEX)-1
      current_state%su%data(k,current_state%column_local_y,current_state%column_local_x) = &
           current_state%su%data(k,current_state%column_local_y,current_state%column_local_x) &
           + du_profile(k)
    end do

  end subroutine apply_time_independent_forcing_to_u

  subroutine apply_time_independent_forcing_to_v(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: k
    real(kind=DEFAULT_PRECISION) :: dtm_scale

    if (constant_forcing_type_v==TENDENCY)then
      dtm_scale=1.0_DEFAULT_PRECISION
    else !  constant_forcing_type_v==(RELAXATION or INCREMENT)
      dtm_scale=1.0_DEFAULT_PRECISION/forcing_timescale_v
    end if
    
    if (constant_forcing_type_v==RELAXATION)then
      dv_profile(:)=dtm_scale * (current_state%global_grid%configuration%vertical%v_force(:) - &
           current_state%global_grid%configuration%vertical%olzvbar(:) )
    else !  constant_forcing_type_v==(TENDENCY or INCREMENT)
      dv_profile(:)=dtm_scale * current_state%global_grid%configuration%vertical%v_force(:)
    end if

    dv_profile_diag=dv_profile_diag+dv_profile

    do k=2,current_state%local_grid%size(Z_INDEX)-1
      current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x) = &
           current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x) &
           + dv_profile(k)
    end do
  end subroutine apply_time_independent_forcing_to_v

  !> Finalises the component 
  !! @current_state Current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    deallocate(theta_profile, q_profile, u_profile, v_profile)
    deallocate(dtheta_profile, dq_profile, du_profile, dv_profile)
    if (allocated(du_profile_diag)) deallocate(du_profile_diag)
    if (allocated(dv_profile_diag)) deallocate(dv_profile_diag)
    if (allocated(dtheta_profile_diag)) deallocate(dtheta_profile_diag)
    if (allocated(dq_profile_diag)) deallocate(dq_profile_diag)

    if (allocated(tend_3d_u)) deallocate(tend_3d_u)
    if (allocated(tend_3d_v)) deallocate(tend_3d_v)
    if (allocated(tend_3d_th)) deallocate(tend_3d_th)
    if (allocated(tend_3d_qv)) deallocate(tend_3d_qv)
    if (allocated(tend_3d_ql)) deallocate(tend_3d_ql)
    if (allocated(tend_3d_qi)) deallocate(tend_3d_qi)
    if (allocated(tend_3d_qr)) deallocate(tend_3d_qr)
    if (allocated(tend_3d_qs)) deallocate(tend_3d_qs)
    if (allocated(tend_3d_qg)) deallocate(tend_3d_qg)
    if (allocated(tend_3d_tabs)) deallocate(tend_3d_tabs)

    if (allocated(tend_pr_tot_u)) deallocate(tend_pr_tot_u)
    if (allocated(tend_pr_tot_v)) deallocate(tend_pr_tot_v)
    if (allocated(tend_pr_tot_th)) deallocate(tend_pr_tot_th)
    if (allocated(tend_pr_tot_qv)) deallocate(tend_pr_tot_qv)
    if (allocated(tend_pr_tot_ql)) deallocate(tend_pr_tot_ql)
    if (allocated(tend_pr_tot_qi)) deallocate(tend_pr_tot_qi)
    if (allocated(tend_pr_tot_qr)) deallocate(tend_pr_tot_qr)
    if (allocated(tend_pr_tot_qs)) deallocate(tend_pr_tot_qs)
    if (allocated(tend_pr_tot_qg)) deallocate(tend_pr_tot_qg)
    if (allocated(tend_pr_tot_tabs)) deallocate(tend_pr_tot_tabs)

  end subroutine finalisation_callback

  !> Averages some diagnostic values across all local horizontal points
  !! @param current_state Current model state
  !! @param diagnostics_summed The diagnostics, which are summed up values from each local point, which need to be averaged
  !! @returns The averaged diagnostics across all local horizontal points
  function get_averaged_diagnostics(current_state, diagnostics_summed)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: diagnostics_summed
    real(kind=DEFAULT_PRECISION), dimension(size(diagnostics_summed)) :: get_averaged_diagnostics

    integer :: horizontal_points

    horizontal_points=current_state%local_grid%size(X_INDEX) * current_state%local_grid%size(Y_INDEX)

    get_averaged_diagnostics(:)=diagnostics_summed(:)/horizontal_points
  end function get_averaged_diagnostics

   !> Save the 3d tendencies coming into the component.
  !! @param current_state Current model state
  !! @param cxn The current slice, x, index
  !! @param cyn The current column, y, index.
  !! @param txn target_x_index
  !! @param tyn target_y_index
  subroutine save_precomponent_tendencies(current_state, cxn, cyn, txn, tyn)
    type(model_state_type), target, intent(in) :: current_state
    integer, intent(in) ::  cxn, cyn, txn, tyn

    ! Save 3d tendency fields upon request (of 3d or profiles) and availability
    if (l_tend_3d_u) then
      tend_3d_u(:,tyn,txn)=current_state%su%data(:,cyn,cxn)
    endif
    if (l_tend_3d_v) then
      tend_3d_v(:,tyn,txn)=current_state%sv%data(:,cyn,cxn)
    endif
    if (l_tend_3d_th) then
      tend_3d_th(:,tyn,txn)=current_state%sth%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qv) then
      tend_3d_qv(:,tyn,txn)=current_state%sq(iqv)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_ql) then
      tend_3d_ql(:,tyn,txn)=current_state%sq(iql)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qi) then
      tend_3d_qi(:,tyn,txn)=current_state%sq(iqi)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qr) then
      tend_3d_qr(:,tyn,txn)=current_state%sq(iqr)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qs) then
      tend_3d_qs(:,tyn,txn)=current_state%sq(iqs)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qg) then
      tend_3d_qg(:,tyn,txn)=current_state%sq(iqg)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_tabs) then
      tend_3d_tabs(:,tyn,txn)=current_state%sth%data(:,cyn,cxn) * current_state%global_grid%configuration%vertical%rprefrcp(:)
    endif

  end subroutine save_precomponent_tendencies

   !> Computation of component tendencies
  !! @param current_state Current model state
  !! @param cxn The current slice, x, index
  !! @param cyn The current column, y, index.
  !! @param txn target_x_index
  !! @param tyn target_y_index
  subroutine compute_component_tendencies(current_state, cxn, cyn, txn, tyn)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  cxn, cyn, txn, tyn

    ! Calculate change in tendency due to component
    if (l_tend_3d_u) then
      tend_3d_u(:,tyn,txn)=current_state%su%data(:,cyn,cxn)       - tend_3d_u(:,tyn,txn)
    endif
    if (l_tend_3d_v) then
      tend_3d_v(:,tyn,txn)=current_state%sv%data(:,cyn,cxn)       - tend_3d_v(:,tyn,txn)
    endif
    if (l_tend_3d_th) then
      tend_3d_th(:,tyn,txn)=current_state%sth%data(:,cyn,cxn)     - tend_3d_th(:,tyn,txn)
    endif
    if (l_tend_3d_qv) then
      tend_3d_qv(:,tyn,txn)=current_state%sq(iqv)%data(:,cyn,cxn) - tend_3d_qv(:,tyn,txn)
    endif
    if (l_tend_3d_ql) then
      tend_3d_ql(:,tyn,txn)=current_state%sq(iql)%data(:,cyn,cxn) - tend_3d_ql(:,tyn,txn)
    endif
    if (l_tend_3d_qi) then
      tend_3d_qi(:,tyn,txn)=current_state%sq(iqi)%data(:,cyn,cxn) - tend_3d_qi(:,tyn,txn)
    endif
    if (l_tend_3d_qr) then
      tend_3d_qr(:,tyn,txn)=current_state%sq(iqr)%data(:,cyn,cxn) - tend_3d_qr(:,tyn,txn)
    endif
    if (l_tend_3d_qs) then
      tend_3d_qs(:,tyn,txn)=current_state%sq(iqs)%data(:,cyn,cxn) - tend_3d_qs(:,tyn,txn)
    endif
    if (l_tend_3d_qg) then
      tend_3d_qg(:,tyn,txn)=current_state%sq(iqg)%data(:,cyn,cxn) - tend_3d_qg(:,tyn,txn)
    endif
    if (l_tend_3d_tabs) then
      tend_3d_tabs(:,tyn,txn)=   &
       current_state%sth%data(:,cyn,cxn) * current_state%global_grid%configuration%vertical%rprefrcp(:)   &
        - tend_3d_tabs(:,tyn,txn)
    endif

   ! Add local tendency fields to the profile total
    if (l_tend_pr_tot_u) then
      tend_pr_tot_u(:)=tend_pr_tot_u(:)   + tend_3d_u(:,tyn,txn)
    endif
    if (l_tend_pr_tot_v) then
      tend_pr_tot_v(:)=tend_pr_tot_v(:)   + tend_3d_v(:,tyn,txn)
    endif
    if (l_tend_pr_tot_th) then
      tend_pr_tot_th(:)=tend_pr_tot_th(:) + tend_3d_th(:,tyn,txn)
    endif
    if (l_tend_pr_tot_qv) then
      tend_pr_tot_qv(:)=tend_pr_tot_qv(:) + tend_3d_qv(:,tyn,txn)
    endif
    if (l_tend_pr_tot_ql) then
      tend_pr_tot_ql(:)=tend_pr_tot_ql(:) + tend_3d_ql(:,tyn,txn)
    endif
    if (l_tend_pr_tot_qi) then
      tend_pr_tot_qi(:)=tend_pr_tot_qi(:) + tend_3d_qi(:,tyn,txn)
    endif
    if (l_tend_pr_tot_qr) then
      tend_pr_tot_qr(:)=tend_pr_tot_qr(:) + tend_3d_qr(:,tyn,txn)
    endif
    if (l_tend_pr_tot_qs) then
      tend_pr_tot_qs(:)=tend_pr_tot_qs(:) + tend_3d_qs(:,tyn,txn)
    endif
    if (l_tend_pr_tot_qg) then
      tend_pr_tot_qg(:)=tend_pr_tot_qg(:) + tend_3d_qg(:,tyn,txn)
    endif
    if (l_tend_pr_tot_tabs) then
      tend_pr_tot_tabs(:)=tend_pr_tot_tabs(:) + tend_3d_tabs(:,tyn,txn)
    endif

  end subroutine compute_component_tendencies


  !> Sets the published field value from the temporary diagnostic values held by this component.
  !! @param field_value Populated with the value of the field
  !! @param real_1d_field Optional one dimensional real of values to publish
  !! @param real_2d_field Optional two dimensional real of values to publish
  subroutine set_published_field_value(field_value, real_1d_field, real_2d_field, real_3d_field)
    type(component_field_value_type), intent(inout) :: field_value
    real(kind=DEFAULT_PRECISION), dimension(:), optional :: real_1d_field
    real(kind=DEFAULT_PRECISION), dimension(:,:), optional :: real_2d_field
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), optional :: real_3d_field

    if (present(real_1d_field)) then
      allocate(field_value%real_1d_array(size(real_1d_field)), source=real_1d_field)
    else if (present(real_2d_field)) then
      allocate(field_value%real_2d_array(size(real_2d_field, 1), size(real_2d_field, 2)), source=real_2d_field)
    else if (present(real_3d_field)) then
      allocate(field_value%real_3d_array(size(real_3d_field, 1), size(real_3d_field, 2), size(real_3d_field, 3)), &
               source=real_3d_field)
    end if
  end subroutine set_published_field_value

  !> Will check a NetCDF status and write to log_log error any decoded statuses
  !! @param status The NetCDF status flag. This is a copy of the routine called in setfluxlook
  subroutine check_forcing_status(status)
    integer, intent(in) :: status

    if (status /= nf90_noerr) then
      call log_log(LOG_ERROR, "NetCDF returned error code of "//trim(nf90_strerror(status)))
    end if
  end subroutine check_forcing_status

  
  !> Reads the dimensions for forcing from the NetCDF file. This routine assumes the forcing  uses only time and height.
  !! @param ncid The NetCDF file id
  !! @param vert_key The vertical coordinate key of the input data 
  !! @param time_dim Number of elements in the time dimension
  !! @param z_dim Number of elements in the vertical dimension
  subroutine read_2d_forcing_dimensions(ncid, vert_key, time_dim, z_dim)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: vert_key
    integer, intent(out) ::  time_dim
    integer, intent(out) ::  z_dim
    integer ::  time_dimid, z_dimid

    call check_forcing_status(nf90_inq_dimid(ncid, TIME_KEY, time_dimid))
    call check_forcing_status(nf90_inquire_dimension(ncid, time_dimid, len=time_dim))

    call check_forcing_status(nf90_inq_dimid(ncid, vert_key, z_dimid))
    call check_forcing_status(nf90_inquire_dimension(ncid, z_dimid, len=z_dim))

  end subroutine read_2d_forcing_dimensions

  !> Reads the variables from the NetCDF forcing file. The 2d variables are assumed to be time and height
  !! @param ncid The id of the NetCDF file
  !! @param time_dim The number of elements in the time dimension
  !! @param time The time data field that is to be read
  !! @param vert_key The vertical coordinate key of the input data
  !! @param z_dim is the number of elements in the height dimension
  !! @param force_2d_key is the string that defines the forcing variable in teh NetCDF file
  !! @param force_2d_var is the forcing data field that is read with dimension (t_dim, z_dim)
  subroutine read_2d_forcing_variables(filename, ncid, time_dim, time, vert_key, z_dim, z_profile, &
       force_2d_key, force_2d_var )

    character(*), intent(in) :: filename
    character(len=*), intent(in) :: force_2d_key, vert_key
    integer, intent(in) :: ncid, time_dim, z_dim
    real(kind=DEFAULT_PRECISION), intent(inout) :: time(:), z_profile(:)
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable, intent(inout) :: force_2d_var
 
    integer :: status, variable_id

    ! Do some checking on the variable contents so that we can deal with different 
    ! variable names or missing variables
    
    ! time and height 
    status=nf90_inq_varid(ncid, TIME_KEY, variable_id)
    if (status==nf90_noerr)then
       call read_single_forcing_variable(ncid, TIME_KEY, data1d=time)
    else
      call log_log(LOG_ERROR, "No recognized time variable found in"//trim(filename))
    end if

    status=nf90_inq_varid(ncid, vert_key, variable_id)
    if (status==nf90_noerr)then
       call read_single_forcing_variable(ncid, vert_key, data1d=z_profile)
    else
      call log_log(LOG_ERROR, "No recognized '"//trim(vert_key)//"' vertical coordinate variable found in"//trim(filename))
    end if
    
    status=nf90_inq_varid(ncid, force_2d_key, variable_id)
    if (status==nf90_noerr)then
       call read_single_forcing_variable(ncid, force_2d_key, data2d=force_2d_var)
    else
       call log_log(LOG_ERROR, "No recognized forcing variable found in"//trim(filename))
   end if

  end subroutine read_2d_forcing_variables

  !> Reads a single variable out of a NetCDF file
  !! @param ncid The NetCDF file id
  !! @param key The variable key (name) to access
  !! @param data1d Optional one dimensional data to read into
  !! @param data3d Optional three dimensional data to read into
  subroutine read_single_forcing_variable(ncid, key, data1d, data2d)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: key
    real(kind=DEFAULT_PRECISION), dimension(:), intent(inout), optional :: data1d
    real(kind=DEFAULT_PRECISION), dimension(:,:), intent(inout), optional :: data2d

    integer :: variable_id
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: sdata

    call check_forcing_status(nf90_inq_varid(ncid, key, variable_id))

    if (.not. present(data1d) .and. .not. present(data2d)) return

    if (present(data1d)) then
      call check_forcing_status(nf90_get_var(ncid, variable_id, data1d))
    else
      ! 2D 
      allocate(sdata(size(data2d,1),size(data2d,2)))
      call check_forcing_status(nf90_get_var(ncid, variable_id, sdata))
      data2d(:,:)=reshape(sdata(:,:),(/size(data2d,1),size(data2d,2)/))
      !deallocate(sdata)
    end if
  end subroutine read_single_forcing_variable

  !> Sets up time-varying forcing profiles
  !! @param current_state Current model state
  !! @param tvdata The time-varying data structure
  !! @param key The variable key (name) to access
  !! @param filename The input NetCDF file name
  !! @param coordinate The vertical coordinate of the input data [ height | pressure ]
  subroutine init_time_varying_forcing(current_state, tvdata, key, filename, coordinate)
    type(model_state_type), target, intent(inout) :: current_state
    type(time_varying_forcing_profile), intent(inout) :: tvdata
    character(len=*), intent(in) :: key, filename, coordinate

    character(STRING_LENGTH) :: vert_key
    integer :: ncid, time_dim_len, vert_dim_len  ! Input file parameters
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: vert_coords   ! contains input vertical coordinates
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: input_forcing
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: vert_grid

    ! Check for valid vertical coordinate specification
    if (trim(coordinate) .eq. 'pressure') then
      vert_key = LEV_KEY
      if (key .eq. WSUBS_KEY) then
        call piecewise_linear_1d(current_state%global_grid%configuration%vertical%zn(:),     &
                                 current_state%global_grid%configuration%vertical%prefn(:),  &
                                 current_state%global_grid%configuration%vertical%z(:),      &
                                 vert_grid)    ! get pressure values on w-levels (z)
      else
        vert_grid = current_state%global_grid%configuration%vertical%prefn(:)
      end if
    else if (trim(coordinate) .eq. 'height') then
      vert_key = Z_KEY
      if (key .eq. WSUBS_KEY) then
        vert_grid = current_state%global_grid%configuration%vertical%z(:)
      else
        vert_grid = current_state%global_grid%configuration%vertical%zn(:)
      end if
    else
      call log_log(LOG_ERROR, "Must specify vertical coordinate for forcing file as 'height' [m] or 'pressure' [Pa] "// &
                              "for '"//trim(key)//"' of file: "//trim(filename))
    end if

    ! open forcing file
    call check_forcing_status(nf90_open(path = trim(filename), mode = nf90_nowrite, ncid = ncid))

    ! read the dimension sizes and allocate space to receive the data
    call read_2d_forcing_dimensions(ncid, trim(vert_key), time_dim_len, vert_dim_len)
    allocate(tvdata%forcing_times(time_dim_len), vert_coords(vert_dim_len), input_forcing(vert_dim_len, time_dim_len), &
             tvdata%forcing_values(size(vert_grid), time_dim_len) )

    ! read the forcing coordinates and data, then close file
    call read_2d_forcing_variables(trim(filename), ncid, time_dim_len, tvdata%forcing_times, &
                                   trim(vert_key), vert_dim_len, vert_coords, &
                                   key, input_forcing)
    call check_forcing_status(nf90_close(ncid))


    ! interpolate forcing levels onto the MONC vertical grid (vert_grid), for all forcing times
    !   Linear gradient extrapolation beyond the input height bounds is implicitly performed.
    !     This behaviour is different from that in piecewise_linear_1d, which does not extrapolate.
    call piecewise_linear_2d(vert_coords, tvdata%forcing_times, input_forcing, &  ! input coordinates and data
                             vert_grid, tvdata%forcing_values)                    ! output vertical coords and data
    ! clean up
    deallocate(vert_coords, input_forcing)

  end subroutine init_time_varying_forcing


end module forcing_mod
