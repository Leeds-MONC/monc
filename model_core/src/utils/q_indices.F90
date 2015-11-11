!> This manages the Q variables and specifically the mapping between names and the index that they are stored at
module q_indices_mod
  use conversions_mod, only : conv_to_string
  use logging_mod, only : LOG_INFO, LOG_ERROR, log_master_log
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer, parameter :: n_maxqs=100 !< Maximum number of Q variables to manage
  integer, parameter :: n_maxqname=100 !< Length to allocate for a Q variable name

  !< Internal storage type for tracking the Q variable states
  type :: q_metadata_type
     character(n_maxqname) :: name
     logical :: l_used=.false.
  end type q_metadata_type
  
  !< Register of the Q variables
  type(q_metadata_type) :: q_register(n_maxqs)

  public q_metadata_type, q_indices_add, get_indices_descriptor, get_max_number_q_indices, &
       get_number_active_q_indices, set_q_index
contains

  !> Sets a Q index to be active at a specific index and sets the name
  !! @param index The index to set the Q index at
  !! @param name The name to set as the Q index name
  subroutine set_q_index(index, name)
    integer, intent(in) :: index
    character(len=*), intent(in) :: name

    q_register(index)%name=adjustr(trim(name))
    q_register(index)%l_used=.true.
  end subroutine set_q_index  

  !> Gets the maximum number of Q indicies
  !! @returns The maximum number of Q indicies
  integer function get_max_number_q_indices()
    get_max_number_q_indices=n_maxqs
  end function get_max_number_q_indices

  !> Gets the number of active Q indicies (i.e. those allocated to specific uses)
  !! @returns The number of active Q indicies
  integer function get_number_active_q_indices()

    integer :: i

    get_number_active_q_indices=0
    do i=1, n_maxqs
      if (q_register(i)%l_used) get_number_active_q_indices=get_number_active_q_indices+1
    end do    
  end function get_number_active_q_indices  

  !> Retrieves the indicies descriptor at a specific location
  !! @param i The index to retrieve the descriptor at
  !! @returns The corresponding descriptor
  function get_indices_descriptor(i)
    integer :: i
    type(q_metadata_type) :: get_indices_descriptor
    
    get_indices_descriptor=q_register(i)
  end function get_indices_descriptor  

  !> Add in a new entry into the register if the name does
  !! not already exist or return the index of the pre-existing variable.
  !! @param name variable name
  !! @param assigning_component name of component which is assigning this variable
  !! @returns The variable index
  integer function q_indices_add(name, assigning_component)
    character(*), intent(in) :: name
    character(*), optional :: assigning_component

    integer :: iname, i_unused

    i_unused=0
    
    do iname=n_maxqs, 1, -1
      if (trim(name) == trim(q_register(iname)%name)) then
        ! Here we find the variable has already been added
        q_indices_add=iname
        return
      end if
      ! Counting backwards we will eventually find the first open slot
      if (.not. q_register(iname)%l_used) then
        i_unused=iname
      end if
    end do

    if (i_unused == 0) then
      call log_master_log(LOG_ERROR, 'Somehow we need to extend q_register')
    end if

    ! Not already defined, so populate the empty slot
    q_indices_add=i_unused
    q_register(q_indices_add)%name = adjustr(trim(name))
    q_register(q_indices_add)%l_used = .true.

    if (present(assigning_component)) then
      call log_master_log(LOG_INFO, 'q variable #'//trim(conv_to_string(q_indices_add))//' is assigned to '//trim(name)//&
           '. Assigned from component: '//trim(assigning_component))
    else
      call log_master_log(LOG_INFO, 'q variable #'//trim(conv_to_string(q_indices_add))//' is assigned to '//trim(name))
    end if

  end function q_indices_add
end module q_indices_mod
