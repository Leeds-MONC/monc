!> Checkpointing NetCDF functionality
!!
!! Both reads and initalises the model based upon a checkpoint and writes (dumps)
!! model status to a checkpoint
module checkpointer_mod
  use datadefn_mod, only : STRING_LENGTH
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type
  use conversions_mod, only : conv_to_string
  use optionsdatabase_mod, only : options_get_string, options_has_key, options_get_integer, options_get_logical
  use logging_mod, only : LOG_INFO, log_master_log, log_master_newline
  use checkpointer_write_checkpoint_mod, only : write_checkpoint_file
  use checkpointer_read_checkpoint_mod, only : read_checkpoint_file
  implicit none

#ifndef TEST_MODE
  private
#endif

  character(len=STRING_LENGTH), save :: checkpoint_file !< The checkpoint write file base name
  logical, save :: unique_per_dump, &     !< Whether to make each model dump a unique filename
       enable_write
  integer, save :: checkpoint_frequency

  public checkpointer_get_descriptor

contains

  !> Provides registry information for the component
  !! @returns The component descriptor that describes this component
  type(component_descriptor_type) function checkpointer_get_descriptor()
    checkpointer_get_descriptor%name="checkpointer"
    checkpointer_get_descriptor%version=0.1
    checkpointer_get_descriptor%initialisation=>initialisation_callback
    checkpointer_get_descriptor%timestep=>timestep_callback
    checkpointer_get_descriptor%finalisation=>finalisation_callback
  end function checkpointer_get_descriptor

  !> Initialisation hook, if appropriate (depends on command line arguments) then will read in
  !! an existing checkpoint file and use this as the basis for model state
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    character(len=STRING_LENGTH) :: internal_write_mode

    checkpoint_frequency=options_get_integer(current_state%options_database, "checkpoint_frequency")
    checkpoint_file=options_get_string(current_state%options_database, "checkpoint_file")
    unique_per_dump=options_get_logical(current_state%options_database, "checkpoint_unique_per_dump")
    internal_write_mode=options_get_string(current_state%options_database, "checkpoint_internal_write")
    if (trim(internal_write_mode) .eq. "always") then
      enable_write=.true.
    else if (trim(internal_write_mode) .eq. "never") then
      enable_write=.false.
    else
      ! Auto mode
      enable_write=.not. current_state%io_server_enabled
    end if

    if (options_has_key(current_state%options_database, "checkpoint")) then
      call read_checkpoint_file(current_state, options_get_string(current_state%options_database, "checkpoint"))
    end if
  end subroutine initialisation_callback

  !> The timestep hook will dump out model state_mod to a checkpoint file
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (enable_write .and. checkpoint_frequency .gt. 0) then
      if (mod(current_state%timestep, checkpoint_frequency) == 0) call perform_checkpoint_dump(current_state)
    end if
  end subroutine timestep_callback

  !> Called on termination to write out the status of the model run to checkpoint
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (enable_write) call perform_checkpoint_dump(current_state)
  end subroutine finalisation_callback

  !> Performs the checkpoint dump and timings. This can be called as part of the timestep or at the end of a
  !! model run
  !! @param current_state The current model state
  subroutine perform_checkpoint_dump(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    character(len=STRING_LENGTH) :: unique_fn
    real :: start_dump_time, end_dump_time
    integer :: ierr

    call cpu_time(start_dump_time)
    if (unique_per_dump) then
      call generate_unique_filename(current_state, unique_fn)
      call write_checkpoint_file(current_state, unique_fn)
    else
      call write_checkpoint_file(current_state, checkpoint_file)
    end if
    ! Barrier here to ensure all processes dumped before log_log stats (is there a better way?)
    call mpi_barrier(current_state%parallel%monc_communicator, ierr)
    call cpu_time(end_dump_time)
    call log_dump_stats(current_state, start_dump_time, end_dump_time)
  end subroutine perform_checkpoint_dump

  !> Will dump out the model dump statistics
  !! @param current_state The current model state
  !! @param startTime The start CPU time of the dump
  !! @param endTime The end CPU time of the dump
  subroutine log_dump_stats(current_state, start_time, end_time)
    type(model_state_type), intent(inout) :: current_state
    real :: start_time, end_time

    call log_master_newline()
    call log_master_log(LOG_INFO, "Model dump completed in "//trim(conv_to_string(int((end_time-start_time)*1000)))//"ms")
  end subroutine log_dump_stats  

  !> Generates a unique filename based upon the base one specified and the number 
  !! of completed timesteps
  !! @param current_state The current model state
  !! @param newName The new name that is produced by this subroutine
  subroutine generate_unique_filename(current_state, new_name)
    type(model_state_type), intent(inout) :: current_state
    character(len=STRING_LENGTH), intent(out) :: new_name

    integer :: dot_posn
    
    dot_posn=index(checkpoint_file, ".")
    if (dot_posn .gt. 0) then
      new_name = checkpoint_file(1:dot_posn-1)
    else
      new_name=checkpoint_file
    end if
    new_name=trim(new_name)//"_"//trim(conv_to_string(current_state%timestep))
    if (dot_posn .gt. 0) then
      new_name=trim(new_name)//checkpoint_file(dot_posn:len(checkpoint_file))
    end if    
  end subroutine generate_unique_filename  
end module checkpointer_mod
