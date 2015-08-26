! Tests the component registry functionality
module test_registry_mod
  use fruit, only : assert_equals, add_fail, assert_true, assert_false, set_unit_name
  use registry_mod, only : get_component_info, get_all_registered_components, register_component, &
       execute_initialisation_callbacks, execute_timestep_callbacks, execute_consolidation_callbacks, &
       execute_modeldump_callbacks, execute_finalisation_callbacks, free_registry, deregister_component
  use state_mod, only : model_state_type
  use monc_component_mod, only : component_descriptor_type
  use collections_mod, only : map_type, c_size, c_value_at, c_key_at
  implicit none

  ! Counters to determine the number of callback calls for testing
  integer :: init_calls = 0
  integer :: timestep_calls = 0
  integer :: consolidation_calls = 0
  integer :: modeldump_calls = 0
  integer :: finalisation_calls = 0

contains
  ! Tests registering components of different names, checks that they all register and then check the info of each one to
  ! ensure that the name and version numbers correspond to the correct values
  subroutine test_register()
    type(component_descriptor_type), pointer :: descriptor
    type(map_type) :: component_info
    class(*), pointer :: data
    integer :: i

    call set_unit_name('test_register')
    do i=1,10
      allocate(descriptor)
      descriptor%version=100-i
      descriptor%name="Test "//str(i)
      call register_component(descriptor)
    end do

    component_info = get_all_registered_components()
    call assert_equals(10, c_size(component_info), "Number of registered components after registrations")

    do i=1,10
      call assert_equals("Test "//str(i), c_key_at(component_info, i))
      data => c_value_at(component_info, i)
      select type(data)
        type is (real)
        call assert_equals(real(100-i), data, "Version number of component at i correct")
        class default
        call add_fail("Unknown type")
      end select
    end do
  end subroutine test_register

  ! Tests the deregistration of components. Will register a load of components, check these are registered and then
  ! deregister each of them and ensure that they have been removed
  subroutine test_deregister()
    type(component_descriptor_type), pointer :: descriptor
    type(map_type) :: component_info
    integer :: i

    do i=1,10
      allocate(descriptor)
      descriptor%version=100-i
      descriptor%name="Test "//str(i)
      call register_component(descriptor)
    end do

    component_info = get_all_registered_components()
    call assert_equals(10, c_size(component_info), "Number of registered components after registrations")

    do i=1,10
      call deregister_component("Test "//str(i))
      component_info = get_all_registered_components()
      call assert_equals(10-i, c_size(component_info), "Component at i has de-registered")
    end do
  end subroutine test_deregister

  ! Tests the component detailed information is correct for each registered component. Also provides incorrect (non registered)
  ! names to ensure that NULL is returned and the registry handles this correctly. The name, version and call back pointers are
  ! all checked for consistency
  subroutine test_component_information()
    type(model_state_type) :: testing_state
    type(component_descriptor_type), pointer :: data
    integer :: init_callbacks, timestep_callbacks, consolidation_callbacks, modeldump_callbacks,&
      finalisation_callbacks, i

    call free_registry() ! Clear the registry to eliminate residue of previous unit tests
    call insert_component_callbacks(init_callbacks, timestep_callbacks, consolidation_callbacks, &
      modeldump_callbacks, finalisation_callbacks)

    do i=1,120
      data => get_component_info("Test "//str(i))
      if (i .le. 100) then
        call assert_true(associated(data), "Testing there is some component information if less than 100")
        call assert_equals("Test "//str(i), data%name, "Compare registered and expected name")
        call assert_equals(real(100-i), data%version, "Compare registered and expected version")
        ! Check each callback is associated if i is within range or not if i is not
        call assert_equals(merge(.true., .false., i .le. init_callbacks), associated(data%initialisation), &
          "Consistency of initialisation call-back")
        call assert_equals(merge(.true., .false., i .le. timestep_callbacks), associated(data%timestep), &
          "Consistency of timestep call-back")
        call assert_equals(merge(.true., .false., i .le. consolidation_callbacks), associated(data%consolidation), &
          "Consistency of consolidation call-back")
        call assert_equals(merge(.true., .false., i .le. modeldump_callbacks), associated(data%modeldump), &
          "Consistency of model dump call-back")
        call assert_equals(merge(.true., .false., i .le. finalisation_callbacks), associated(data%finalisation), &
          "Consistency of finalisation all-back")
      else
        call assert_false(associated(data), "No component if greater than 100") ! Not registered i>100 so should be NULL
      end if
    end do

  end subroutine test_component_information

  ! Will register a hundred components and then de-register components from 25 to 50 (inclusive, so 26
  ! deregistrations.) Will then execute the callback stages and ensure that the number of calls is
  ! consistent with the expected number taking into account the deregistrations
  subroutine test_component_removal_callbacks()
    type(model_state_type) :: testing_state
    integer :: init_callbacks, timestep_callbacks, consolidation_callbacks, modeldump_callbacks,&
      finalisation_callbacks, i

    call clear_counters()
    call free_registry() ! Clear the registry to eliminate residue of previous unit tests
    call insert_component_callbacks(init_callbacks, timestep_callbacks, consolidation_callbacks, &
      modeldump_callbacks, finalisation_callbacks)

    do i=25,50
      ! Now deregister 25-50 callbacks
      call deregister_component("Test "//str(i))
    end do

    ! In reference to the 25-50 inclusive removals, recalculate the expected number of each stage's callback
    init_callbacks = calculate_remaining_calls(init_callbacks, 25, 50)
    timestep_callbacks = calculate_remaining_calls(timestep_callbacks, 25, 50)
    consolidation_callbacks = calculate_remaining_calls(consolidation_callbacks, 25, 50)
    modeldump_callbacks = calculate_remaining_calls(modeldump_callbacks, 25, 50)
    finalisation_callbacks = calculate_remaining_calls(finalisation_callbacks, 25, 50)

    ! Call the stages
    call execute_initialisation_callbacks(testing_state)
    call execute_timestep_callbacks(testing_state)
    call execute_consolidation_callbacks(testing_state)
    call execute_modelDump_callbacks(testing_state)
    call execute_finalisation_callbacks(testing_state)

    ! Check number of calls in each stages's callbacks are appropriate
    call assert_equals(init_callbacks, init_calls, "Number of initialisation call-backs post removal")
    call assert_equals(timestep_callbacks, timestep_calls, "Number of timestep call-backs post removal")
    call assert_equals(consolidation_callbacks, consolidation_calls, "Number of consolidation call-backs post removal")
    call assert_equals(modeldump_callbacks, modeldump_calls, "Number of model dump call-backs post removal")
    call assert_equals(finalisation_callbacks, finalisation_calls, "Number of finalisation call-backs post removal")

  end subroutine test_component_removal_callbacks

  ! Will register a number of components, execute the state callbacks and ensure that the number of callbacks in
  ! each stage are consistent which the number registered
  subroutine test_component_callbacks()
    type(model_state_type) :: testing_state
    integer :: init_callbacks, timestep_callbacks, consolidation_callbacks, modeldump_callbacks, finalisation_callbacks

    call clear_counters()
    call free_registry() ! Clear the registry to remove any residue from previous unit tests
    call insert_component_callbacks(init_callbacks, timestep_callbacks, consolidation_callbacks, &
         modeldump_callbacks, finalisation_callbacks)

    ! Call the callbacks for each stage
    call execute_initialisation_callbacks(testing_state)
    call execute_timestep_callbacks(testing_state)
    call execute_consolidation_callbacks(testing_state)
    call execute_modelDump_callbacks(testing_state)
    call execute_finalisation_callbacks(testing_state)

    ! Check that the number of callback calls is consistent with what we expected for each stage
    call assert_equals(init_callbacks, init_calls, "Number of initialisation call-backs")
    call assert_equals(timestep_callbacks, timestep_calls, "Number of timestep call-backs")
    call assert_equals(consolidation_callbacks, consolidation_calls, "Number of consolidation call-backs")
    call assert_equals(modeldump_callbacks, modeldump_calls, "Number of model dump call-backs")
    call assert_equals(finalisation_callbacks, finalisation_calls, "Number of finalisation call-backs")
  end subroutine test_component_callbacks

  ! Will register a number of call backs and then reregister identical named callbacks. This tests that
  ! if someone reregisters a callback then the new registration takes the place of the previous callback
  subroutine test_component_replacement_callbacks()
    type(model_state_type), target :: testing_state
    integer :: init_callbacks, timestep_callbacks, consolidation_callbacks, modeldump_callbacks, finalisation_callbacks

    call clear_counters()
    call free_registry() ! Clear the registry to remove any residue of previous unit tests
    call insert_component_callbacks(init_callbacks, timestep_callbacks, consolidation_callbacks, &
         modeldump_callbacks, finalisation_callbacks)
    ! Recreate our callbacks
    call generate_dummy_callbacks(init_callbacks, timestep_callbacks, consolidation_callbacks, &
         modeldump_callbacks, finalisation_callbacks)

    ! Execute callbacks for each stage
    call execute_initialisation_callbacks(testing_state)
    call execute_timestep_callbacks(testing_state)
    call execute_consolidation_callbacks(testing_state)
    call execute_modelDump_callbacks(testing_state)
    call execute_finalisation_callbacks(testing_state)

    ! Check the number of calls to each callback is consistent with what we initially expected
    call assert_equals(init_callbacks, init_calls, "Number of initialisation call-backs post replacement")
    call assert_equals(timestep_callbacks, timestep_calls, "Number of timestep call-backs post replacement")
    call assert_equals(consolidation_callbacks, consolidation_calls, "Number of consolidation call-backs post replacement")
    call assert_equals(modeldump_callbacks, modeldump_calls, "Number of model dump call-backs post replacement")
    call assert_equals(finalisation_callbacks, finalisation_calls, "Number of finalisation call-backs post replacement")
  end subroutine test_component_replacement_callbacks

  ! A helper function which calculates the remaining calls from an initial number and the start and
  ! end numbers (inclusive) of components that have been deregistered
  integer function calculate_remaining_calls(orig_value, a, b)
    integer, intent(in) :: orig_value, a, b

    if (orig_value .ge. a) then
      calculate_remaining_calls = orig_value - merge(orig_value - a+1, a+1, orig_value .le. b)
    else
      calculate_remaining_calls = orig_value
    end if
  end function calculate_remaining_calls

  ! A helper subroutine to insert the component callbacks. It generates a random number of callbacks
  ! for each stage and then inserts components based upon these.
  subroutine insert_component_callbacks(init_callbacks, timestep_callbacks, consolidation_callbacks,&
    modelDump_callbacks, finalisation_callbacks)

    integer, intent(out) :: init_callbacks, timestep_callbacks, consolidation_callbacks, &
         modelDump_callbacks, finalisation_callbacks
    real :: r

    call init_random_seed()

    call random_number(r)
    init_callbacks = int(r*99)+1
    call random_number(r)
    timestep_callbacks = int(r*99)+1
    call random_number(r)
    consolidation_callbacks = int(r*99)+1
    call random_number(r)
    modelDump_callbacks = int(r*99)+1
    call random_number(r)
    finalisation_callbacks = int(r*99)+1

    call generate_dummy_callbacks(init_callbacks, timestep_callbacks, consolidation_callbacks, &
         modelDump_callbacks, finalisation_callbacks)
  end subroutine insert_component_callbacks

  ! Registers a hundred dummy components and specific callbacks depending upon the numbers provided
  subroutine generate_dummy_callbacks(init_callbacks, timestep_callbacks, consolidation_callbacks, &
       modelDump_callbacks, finalisation_callbacks)

    integer, intent(in) :: init_callbacks, timestep_callbacks, consolidation_callbacks, &
         modelDump_callbacks, finalisation_callbacks
    type(component_descriptor_type), pointer :: descriptor
    integer :: i

    do i=1,100
      allocate(descriptor)
      descriptor%name="Test "//str(i)
      descriptor%version=100-i
      if (i .le. init_callbacks) descriptor%initialisation=>internal_test_init
      if (i .le. timestep_callbacks) descriptor%timestep=>internal_test_timestep
      if (i .le. consolidation_callbacks) descriptor%consolidation=>internal_test_consolidation
      if (i .le. modelDump_callbacks) descriptor%modeldump=>internal_test_modelDump
      if (i .le. finalisation_callbacks) descriptor%finalisation=>internal_test_finalisation
      call register_component(descriptor)
    end do
  end subroutine generate_dummy_callbacks

  ! Test callback function, will increment the init integer value
  subroutine internal_test_init(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    init_calls = init_calls + 1
  end subroutine internal_test_init

  ! Test callback function, will increment the timestep integer value
  subroutine internal_test_timestep(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    timestep_calls = timestep_calls + 1
  end subroutine internal_test_timestep

  ! Test callback function, will increment the consolidation integer value
  subroutine internal_test_consolidation(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    consolidation_calls = consolidation_calls + 1
  end subroutine internal_test_consolidation

  ! Test callback function, will increment the modeldump integer value
  subroutine internal_test_modeldump(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    modeldump_calls = modeldump_calls + 1
  end subroutine internal_test_modeldump

  ! Test callback function, will increment the finalisation integer value
  subroutine internal_test_finalisation(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    finalisation_calls = finalisation_calls + 1
  end subroutine internal_test_finalisation

  ! Helper function to convert an integer into a string
  character(len=15) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
  end function str

  ! Helper function to seed the random number generator
  subroutine init_random_seed()
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)

    deallocate(seed)
  end subroutine init_random_seed

  ! Clears the call back counters
  subroutine clear_counters()
    init_calls = 0
    timestep_calls = 0
    consolidation_calls = 0
    modeldump_calls = 0
    finalisation_calls = 0
  end subroutine clear_counters

end module test_registry_mod

! The driver for testing the registry
program test_registry_driver
  use fruit, only : init_fruit, run_test_case, fruit_summary
  use test_registry_mod, only : test_register, test_deregister, test_component_information, test_component_callbacks, &
       test_component_replacement_callbacks, test_component_removal_callbacks

  implicit none

  call init_fruit
  call run_test_case(test_register, "Component registration")
  call run_test_case(test_deregister, "Component de-registration")
  call run_test_case(test_component_information, "Retrieval of component information")
  call run_test_case(test_component_callbacks, "Component call-backs")
  call run_test_case(test_component_replacement_callbacks, "Component replacement with call-backs")
  call run_test_case(test_component_removal_callbacks, "Component removal with call-backs")
  call fruit_summary
end program test_registry_driver
