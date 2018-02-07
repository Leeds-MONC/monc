!> Common checkpoint functionality which is used by reader and writers to NetCDF checkpoints
module checkpointer_common_mod
#ifndef TEST_MODE
  use netcdf, only : nf90_ebaddim, nf90_enotatt, nf90_enotvar, nf90_noerr, nf90_strerror
#else
  use dummy_netcdf_mod, only : nf90_ebaddim, nf90_enotatt, nf90_enotvar, nf90_noerr, nf90_strerror
#endif
  use logging_mod, only : LOG_ERROR, log_log
  implicit none

  character(len=*), parameter ::  X_DIM_KEY = "x", &                 !< X dimension/variable key
                                  Y_DIM_KEY="y", &                   !< Y dimension/variable key
                                  Z_DIM_KEY="z", &                   !< Z dimension/variable key                                  
                                  ZN_DIM_KEY="zn", &                  
                                  Q_DIM_KEY="q", &
                                  U_KEY = "u_nogal", &                          !< U variable NetCDF key
                                  V_KEY = "v_nogal", &                          !< V variable NetCDF key
                                  W_KEY = "w", &                          !< W variable NetCDF key
                                  Q_KEY = "q", &                          !< Q variable NetCDF key
                                  ZU_KEY = "zu", &
                                  ZV_KEY = "zv", &
                                  ZW_KEY = "zw" , &
                                  ZQ_KEY = "zq", &
                                  X_KEY = "x", &
                                  Y_KEY = "y", &
                                  Z_KEY = "z", &
                                  ZN_KEY="zn", &
                                  TH_KEY = "th", &                         !< Theta variable NetCDF key
                                  ZTH_KEY = "zth", &
                                  P_KEY = "p", &                          !< Pressure variable NetCDF key
                                  TIMESTEP="timestep", &                  !< Timestep NetCDF key
                                  TIME_KEY="time",&
                                  DTM_KEY="dtm",&
                                  DTM_NEW_KEY="dtm_new",&
                                  ABSOLUTE_NEW_DTM_KEY="absolute_new_dtm",&
                                  UGAL="ugal",&
                                  VGAL="vgal",&
                                  EMPTY_DIM_KEY="empty_dim", &            !< Empty dimension key
                                  KEY_VALUE_PAIR_KEY="kvp", &  !< Key-value pair dimension key
                                  OPTIONS_DIM_KEY="number_options", &       !< Options dimension key
                                  OPTIONS_KEY="options_database", &                !< Options variable key
                                  STRING_DIM_KEY="string",&          !< String dimension key
                                  TITLE_ATTRIBUTE_KEY="title",&
                                  CREATED_ATTRIBUTE_KEY="created",&
                                  NQFIELDS="nqfields", &
                                  Q_INDICES_DIM_KEY="active_q_indicies", &
                                  Q_INDICES_KEY="q_indicies", &
                                  X_RESOLUTION="x_resolution", &
                                  Y_RESOLUTION="y_resolution", &
                                  X_TOP="x_top", &
                                  Y_TOP="y_top", &
                                  X_BOTTOM="x_bottom", &
                                  Y_BOTTOM="y_bottom", &
                                  Q_FIELD_ANONYMOUS_NAME="q_qfield", &
                                  ZQ_FIELD_ANONYMOUS_NAME="zq_qfield", &
                                  THREF="thref", &
                                  OLUBAR="olubar", &
                                  OLZUBAR="olzubar", &
                                  OLVBAR="olvbar", &
                                  OLZVBAR="olzvbar", &
                                  OLTHBAR="olthbar", &
                                  OLZTHBAR="olzthbar", &
                                  OLQBAR="olqbar", &
                                  OLQBAR_ANONYMOUS_NAME="olqbar_qfield", &
                                  OLZQBAR="olzqbar", &
                                  OLZQBAR_ANONYMOUS_NAME="olzqbar_qfield", &
                                  RAD_LAST_TIME_KEY="rad_last_time", &
                                  STH_LW_KEY="sth_lw", &
                                  STH_SW_KEY="sth_sw"

  integer, parameter :: MAX_STRING_LENGTH = 100   !< Maximum string length (stored size)

contains

  !> Will check a NetCDF status and write to log_log error any decoded statuses. Can be used to decode
  !! whether a dimension or variable exists within the NetCDF data file
  !! @param status The NetCDF status flag
  !! @param foundFlag Whether the field has been found or not
  subroutine check_status(status, found_flag)
    integer, intent(in) :: status
    logical, intent(out), optional :: found_flag

    if (present(found_flag)) then
      found_flag = status /= nf90_ebaddim .and. status /= nf90_enotatt .and. status /= nf90_enotvar
      if (.not. found_flag) return
    end if

    if (status /= nf90_noerr) then
      call log_log(LOG_ERROR, "NetCDF returned error code of "//trim(nf90_strerror(status)))
    end if
  end subroutine check_status

  !> Removes NetCDF C style null termination of string. This is placed right at the end, after any
  !! spaces so trim will not actually trim any spaces due to null terminator
  !! @param netCDFString The NetCDF string to remove the null terminator from which is modified
  subroutine remove_null_terminator_from_string(net_cdf_string)
    character(len=*), intent(inout) :: net_cdf_string
    integer :: i
    do i=1,len(net_cdf_string)
      if (iachar(net_cdf_string(i:i)) == 0) then
        net_cdf_string(i:len(net_cdf_string)) = ' '
        exit
      end if
    end do
  end subroutine remove_null_terminator_from_string
end module checkpointer_common_mod
