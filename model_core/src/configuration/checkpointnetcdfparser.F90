!> Loads in the configuration stored in a NetCDF checkpoint file for the model to start from
module configuration_checkpoint_netcdf_parser_mod
  use datadefn_mod, only : STRING_LENGTH
  use collections_mod, only : hashmap_type
  use netcdf, only : NF90_NOWRITE, NF90_NETCDF4, NF90_MPIIO, NF90_NOERR, nf90_strerror, nf90_open, nf90_close, &
       nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, nf90_get_var
  use logging_mod, only : LOG_ERROR, log_master_log
  use conversions_mod, only : conv_is_integer, conv_to_integer, conv_is_real, conv_to_real, conv_is_logical, conv_to_logical, &
       conv_single_real_to_double
  use optionsdatabase_mod, only : options_add
  use mpi, only : MPI_INFO_NULL
  use netcdf_misc_mod, only : check_netcdf_status
  implicit none

#ifndef TEST_MODE
  private
#endif

  character(len=*), parameter :: OPTIONS_KEY="options_database", &  !< The options key which references the configuration
       OPTIONS_DIM_KEY="number_options" !< Options dimension key

  public parse_configuration_checkpoint_netcdf
contains

  !> Will parse the NetCDF checkpoint file and loads the configuration into the options database
  !! @param options_database The options database
  !! @param checkpoint_name Name of the checkpoint file
  !! @param communicator MPI communicator for parallel IO
  subroutine parse_configuration_checkpoint_netcdf(options_database, checkpoint_name, communicator)
    type(hashmap_type), intent(inout) :: options_database
    character(*), intent(in) :: checkpoint_name
    integer, intent(in) :: communicator

    integer :: ncid

    call check_netcdf_status(nf90_open(path = checkpoint_name, mode = NF90_NOWRITE, ncid = ncid))
    call load_options(options_database, ncid)
    call check_netcdf_status(nf90_close(ncid))
  end subroutine parse_configuration_checkpoint_netcdf

  !> Will read in and initialise the options database from the contents of the checkpoint file
  !! @param options_database The options database to store the values into
  !! @param ncid The NetCDF file id
  subroutine load_options(options_database, ncid)
    type(hashmap_type), intent(inout) :: options_database
    integer, intent(in) :: ncid

    integer :: i, options_id, number_options
    character(len=STRING_LENGTH) :: key, value

    number_options=get_number_of_options(ncid)
    call check_netcdf_status(nf90_inq_varid(ncid, OPTIONS_KEY, options_id))

    do i=1, number_options
      call check_netcdf_status(nf90_get_var(ncid, options_id, key, (/ 1, 1, i /)))
      call check_netcdf_status(nf90_get_var(ncid, options_id, value, (/ 1, 2, i /)))
      ! NetCDF does C style null termination right at the end, need to remove this so can trim spaces etc
      call remove_null_terminator_from_string(key)
      call remove_null_terminator_from_string(value)
      if (conv_is_integer(trim(value))) then
        call options_add(options_database, trim(key), conv_to_integer(trim(value)))
      else if (conv_is_real(trim(value))) then
        call options_add(options_database, trim(key), conv_single_real_to_double(conv_to_real(trim(value))))
      else if (conv_is_logical(trim(value))) then
        call options_add(options_database, trim(key), conv_to_logical(trim(value)))
      else
        call options_add(options_database, trim(key), trim(value))
      end if
    end do
  end subroutine load_options

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

  !> Retrieves the number of option key-value pairs that are present in the checkpoint file
  !! @param ncid The NetCDF file id
  !! @returns The number of options that will be read in
  integer function get_number_of_options(ncid)
    integer, intent(in) :: ncid

    integer :: options_dimid, options_dim
    call check_netcdf_status(nf90_inq_dimid(ncid, OPTIONS_DIM_KEY, options_dimid))
    call check_netcdf_status(nf90_inquire_dimension(ncid, options_dimid, len=options_dim))
    get_number_of_options=options_dim
  end function get_number_of_options
end module configuration_checkpoint_netcdf_parser_mod
