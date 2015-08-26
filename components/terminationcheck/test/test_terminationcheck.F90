! Tests the termination checking and determination functionality
module test_terminationcheck_mod
  use fruit, only : assert_equals
  use terminationcheck_mod, only : timestep_callback, consolidation_callback, modeldump_callback, options_add
  use state_mod, only : model_state_type
  implicit none

  integer, parameter :: MAX_TESTS = 1000  ! Number of tests to run for each case, some will require termination some not

contains

  ! Tests the timestep termination checking functionality
  subroutine test_timestep_callback()
    integer :: i, nn_timesteps
    type(model_state_type) :: current_state
    real :: r

    call init_random_seed

    do i=1,MAX_TESTS
      call random_number(r)
      current_state%timestep = int(r*10000)+1
      call random_number(r)
      nn_timesteps= int(r*9)+1

      current_state%last_timestep_column = .true. ! Forces the check (only once per timestep)
      current_state%continue_timestep=.true.
      call options_add(current_state%options_database, "nn_timesteps", nn_timesteps)
      call timestep_callback(current_state)
      call assert_equals(mod(current_state%timestep, nn_timesteps) /= 0, current_state%continue_timestep,&
        "Timestep completion consistent with expectations")
    end do
  end subroutine test_timestep_callback

  ! Tests the consolidation stage termination checking functionality
  subroutine test_consolidation_callback()
    integer :: i, nn_timesteps, nn_consolidation, multiplier
    type(model_state_type) :: current_state
    real :: r

    call init_random_seed

    do i=1,MAX_TESTS
      call random_number(r)
      nn_timesteps= int(r*9)+1
      call random_number(r)
      multiplier = int(r*199)+1
      call random_number(r)
      nn_consolidation = int(r*9)+1

      current_state%timestep = nn_timesteps * multiplier ! timestep already terminated hence multiplication
      current_state%continue_consolidation=.true.
      call options_add(current_state%options_database, "nn_timesteps", nn_timesteps)
      call options_add(current_state%options_database, "nn_consolidation", nn_consolidation)
      call consolidation_callback(current_state)
      call assert_equals(mod(current_state%timestep / nn_timesteps, nn_consolidation) /= 0, current_state%continue_consolidation,&
        "Consolidation completion consistent with expectations")
    end do
  end subroutine test_consolidation_callback

  ! Tests the model dump termination checking functionality
  subroutine test_modeldump_callback()
    integer :: i, nn_timesteps, nn_consolidation, nn_modeldump, multiplier
    type(model_state_type) :: current_state
    real :: r

    call init_random_seed

    do i=1,MAX_TESTS
      call random_number(r)
      nn_timesteps= int(r*9)+1
      call random_number(r)
      multiplier = int(r*39)+1
      call random_number(r)
      nn_consolidation = int(r*9)+1
      call random_number(r)
      nn_modeldump = int(r*9)+1

      current_state%timestep = nn_timesteps * nn_consolidation * multiplier  ! timestep and consolidation already terminated hence multiplication
      current_state%continue_modeldump=.true.
      call options_add(current_state%options_database, "nn_timesteps", nn_timesteps)
      call options_add(current_state%options_database, "nn_consolidation", nn_consolidation)
      call options_add(current_state%options_database, "nn_modeldump", nn_modeldump)
      call modeldump_callback(current_state)
      call assert_equals(mod(current_state%timestep / nn_timesteps / nn_consolidation, nn_modeldump) /= 0,&
       current_state%continue_modeldump, "Model dump completion consistent with expectations")
    end do
  end subroutine test_modeldump_callback

  ! Helper subroutine to initialise the random seed (based on the clock)
  subroutine init_random_seed()
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put=seed)

    deallocate(seed)
  end subroutine init_random_seed
end module test_terminationcheck_mod

! Driver for termination checking tests
program test_terminationcheck_driver
  use fruit, only : init_fruit, run_test_case, fruit_summary
  use test_terminationcheck_mod, only : test_timestep_callback, test_consolidation_callback, test_modeldump_callback

  implicit none

  call init_fruit
  call run_test_case(test_timestep_callback, "Test time-stepping completion")
  call run_test_case(test_consolidation_callback, "Test consolidation completion")
  call run_test_case(test_modeldump_callback, "Test model dump completion")
  call fruit_summary
end program test_terminationcheck_driver
