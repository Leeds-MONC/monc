! Tests the logging_mod utility functions
module test_logging_mod
    use fruit, only : assert_equals, assert_not_equals, assert_true
    use logging_mod, only : LOG_ERROR, LOG_WARN, LOG_INFO, LOG_DEBUG, ERROR_PREAMBLE, &
         WARN_PREAMBLE, INFO_PREAMBLE, DEBUG_PREAMBLE, log_set_logging_level, &
         log_log, get_logging_preamble, log_get_logging_level, error_is_fatal
    implicit none

    contains
      
      ! Tests the logging_mod by sending messages into a string. It will run at different 
      ! levels to ensure that
      ! messages are written iff the logging_mod level is appropriate
    subroutine test_logging_level
      character(len=30) :: outpt
      
      integer i, j, space
      
      do i=1,4
         call log_set_logging_level(i)
         call assert_equals(i, log_get_logging_level(), "Logging level has been set correctly")
         do j=1,4
            call log_log(j, str(j)//" level", outpt)
            if (j .le. i) then
               space = index(outpt, "] ")
               call assert_true(space .gt. 0, "The logging_mod message terminates")
               outpt = outpt(space+2 : len(outpt))
               call assert_equals(str(j)//" level", outpt, &
                    "Message included for specific logging_mod level")
            else
               call assert_not_equals(str(j)//" level", outpt, &
                    "Message ignored for specific logging_mod level")
            end if
         end do
      end do
    end subroutine test_logging_level
    
    ! Tests that the specific logging_mod levels produce the correct preambles
    subroutine test_logging_preamble()
      character(len=7) :: preamble
      
      preamble = get_logging_preamble(LOG_ERROR)
      call assert_equals(ERROR_PREAMBLE, preamble)

      preamble = get_logging_preamble(LOG_WARN)
      call assert_equals(WARN_PREAMBLE, preamble)

      preamble = get_logging_preamble(LOG_INFO)
      call assert_equals(INFO_PREAMBLE, preamble)

      preamble = get_logging_preamble(LOG_DEBUG)
      call assert_equals(DEBUG_PREAMBLE, preamble)
    end subroutine test_logging_preamble

    ! Helper function to convert an integer into a string
    character(len=15) function str(k)
      integer, intent(in) :: k
      write (str, *) k
      str = adjustl(str)
    end function str

end module test_logging_mod

! Driver for logging_mod utility tests
program test_logging_mod_driver
  use fruit, only : init_fruit, run_test_case, fruit_summary
  use test_logging_mod, only : test_logging_level, test_logging_preamble

  implicit none
  logical :: error_fatal
  call init_fruit
  error_fatal = .false.  ! Error messages continue to allow for testing
  call run_test_case(test_logging_level, "Test logging_mod level")
  call run_test_case(test_logging_preamble, "Test logging_mod preamble")
  call fruit_summary
end program test_logging_mod_driver
