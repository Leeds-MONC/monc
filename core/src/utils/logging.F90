!> Logging utility
module logging_mod
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer, parameter, private :: MASTER_PROCESS = 0

  integer, parameter, public :: LOG_ERROR = 1 !< Only log ERROR messages
  integer, parameter, public :: LOG_WARN = 2  !< Log WARNING and ERROR messages
  integer, parameter, public :: LOG_INFO = 3  !< Log INFO, WARNING and ERROR messages
  integer, parameter, public :: LOG_DEBUG = 4 !< Log DEBUG, INFO, WARNING and ERROR messages

  character(len=*), parameter ::  ERROR_PREAMBLE = "[ERROR]",&  !< Error message preamble
                                  WARN_PREAMBLE = "[WARN]",& !< Warning message preable
                                  INFO_PREAMBLE = "[INFO]",&    !< Info message preamble
                                  DEBUG_PREAMBLE = "[DEBUG]"    !< Debug message preamble

  !> Private current logging state, the default is INFO level
  integer :: current_log_level=LOG_INFO
  !> Whether an error log message is fatal and stops execution
  logical, parameter :: error_is_fatal=.true.
  !> Whether or not I am the master process
  logical, save :: i_am_master=.true.

  public initialise_logging, log_log, log_newline, log_master_log, log_master_newline, log_is_master, &
       log_set_logging_level, log_get_logging_level

contains

  !> Initialises the logging. This is done to make it easier for master logging only, so that we don't have to pass the
  !! process ID in each time
  !! @param pid The process Id
  subroutine initialise_logging(pid)
    integer, intent(in) :: pid

    i_am_master=pid==MASTER_PROCESS
  end subroutine initialise_logging  

  !> Will log just from the master process
  !! @param level The logging level
  !! @param message The logging message
  !! @param myId My process id
  subroutine log_master_log(level, message)
    integer, intent(in) :: level
    character(len=*), intent(in) :: message

    if (i_am_master) then
      call log_log(level, message)
    else
      if (error_is_fatal .and. level .eq. LOG_ERROR) call abort()
    end if
  end subroutine log_master_log
  
  !> The master process will log a new line to stdio
  subroutine log_master_newline()
    if (i_am_master) write(*,*) ""
  end subroutine log_master_newline

  !> Determines whether the process is the master logging process. This might be preferable rather than
  !! calling masterLog due to having to construct the message for masterLog regardless
  !! @param pid The process id to test
  logical function log_is_master()
    log_is_master = i_am_master
  end function log_is_master  

  !> Logs a message at the specified level. If the level is above the current level then the message is ignored.
  !! If applicable then the message will be written to stdout or a provided string
  !! @param level The logging level
  !! @param message The message to log
  !! @param str Optional string to write the message into (rather than stdout)
  subroutine log_log(level, message, str)
    integer, intent(in) :: level
    character(len=*), intent(in) :: message
    character(len=*), intent(out), optional :: str

    if (level .le. current_log_level) then
      if (present(str)) then
        write(str,"(A,A,A)") trim(get_logging_preamble(level)), " ",message
      else
        write(*,"(A,A,A)") trim(get_logging_preamble(level)), " ", message
      end if
    end if
    if (error_is_fatal .and. level .eq. LOG_ERROR) call abort()
  end subroutine log_log

  !> Will log a new line to the stdout
  subroutine log_newline()
    write(*,*) ""
  end subroutine log_newline  

  !> Returns the string preamble that applies to the specific logging level
  !! @param level The logging level to translate into a string
  !! @returns String representation of specific level
  character(len=7) function get_logging_preamble(level)
    integer, intent(in) :: level

    if (level .eq. LOG_ERROR) then
      get_logging_preamble = ERROR_PREAMBLE
    else if (level .eq. LOG_WARN) then
      get_logging_preamble = WARN_PREAMBLE
    else if (level .eq. LOG_INFO) then
      get_logging_preamble = INFO_PREAMBLE
    else
      get_logging_preamble = DEBUG_PREAMBLE
    end if
  end function get_logging_preamble

  !> Sets the logging level, messages with less priority will be ignored
  !! @param level The logging level to adopt
  subroutine log_set_logging_level(level)
    integer, intent(in) :: level

    current_log_level = level
  end subroutine log_set_logging_level

  !> Retrieves the current logging level
  !! @returns The current logging level
  integer function log_get_logging_level()
    log_get_logging_level = current_log_level
  end function log_get_logging_level
end module logging_mod
