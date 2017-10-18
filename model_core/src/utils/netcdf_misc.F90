!> NetCDF misc functionality which can be shared between modules that work with NetCDF files
module netcdf_misc_mod
  use datadefn_mod, only : STRING_LENGTH
  use netcdf, only : nf90_ebaddim, nf90_enotatt, nf90_enotvar, nf90_noerr, nf90_strerror
  use logging_mod, only : LOG_ERROR, log_log
  implicit none

#ifndef TEST_MODE
  private
#endif

  public check_netcdf_status
contains

  !> Will check a NetCDF status and write to log_log error any decoded statuses. Can be used to decode
  !! whether a dimension or variable exists within the NetCDF data file
  !! @param status The NetCDF status flag
  !! @param foundFlag Whether the field has been found or not
  subroutine check_netcdf_status(status, found_flag, error_message_details)
    integer, intent(in) :: status
    character(len=STRING_LENGTH), intent(in), optional :: error_message_details
    logical, intent(out), optional :: found_flag

    if (present(found_flag)) then
      found_flag = status /= nf90_ebaddim .and. status /= nf90_enotatt .and. status /= nf90_enotvar
      if (.not. found_flag) return
    end if

    if (status /= nf90_noerr) then
       if (present(error_message_details)) then
          call log_log(LOG_ERROR, "NetCDF returned error code of "//trim(nf90_strerror(status))&
             //" ("//trim(error_message_details)//")")
       else
          call log_log(LOG_ERROR, "NetCDF returned error code of "//trim(nf90_strerror(status)))
       endif
    end if
  end subroutine check_netcdf_status
end module netcdf_misc_mod
