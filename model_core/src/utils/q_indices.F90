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


  ! standard names
  type :: standard_q_names_type
     character(LEN=6 ) :: VAPOUR                 = 'vapour' 
     character(LEN=17) :: CLOUD_LIQUID_MASS      = 'cloud_liquid_mass'
     character(LEN=9 ) :: RAIN_MASS              = 'rain_mass'
     character(LEN=8 ) :: ICE_MASS               = 'ice_mass'
     character(LEN=9 ) :: SNOW_MASS              = 'snow_mass'
     character(LEN=12) :: GRAUPEL_MASS           = 'graupel_mass'
     character(LEN=19) :: CLOUD_LIQUID_NUMBER    = 'cloud_liquid_number'
     character(LEN=11) :: RAIN_NUMBER            = 'rain_number'
     character(LEN=10) :: ICE_NUMBER             = 'ice_number'
     character(LEN=11) :: SNOW_NUMBER            = 'snow_number'
     character(LEN=14) :: GRAUPEL_NUMBER         = 'graupel_number'
     character(LEN=17) :: RAIN_THIRD_MOMENT      = 'rain_third_moment'
     character(LEN=17) :: SNOW_THIRD_MOMENT      = 'snow_third_moment'
     character(LEN=20) :: GRAUPEL_THIRD_MOMENT   = 'graupel_third_moment'
     character(LEN=15) :: AITKEN_SOL_MASS        = 'aitken_sol_mass'
     character(LEN=17) :: AITKEN_SOL_NUMBER      = 'aitken_sol_number'
     character(LEN=14) :: ACCUM_SOL_MASS         = 'accum_sol_mass'
     character(LEN=16) :: ACCUM_SOL_NUMBER       = 'accum_sol_number'
     character(LEN=15) :: COARSE_SOL_MASS        = 'coarse_sol_mass'
     character(LEN=17) :: COARSE_SOL_NUMBER      = 'coarse_sol_number'
     character(LEN=17) :: ACTIVE_SOL_LIQUID      = 'active_sol_liquid'
     character(LEN=15) :: ACTIVE_SOL_RAIN        = 'active_sol_rain'
     character(LEN=16) :: COARSE_DUST_MASS       = 'coarse_dust_mass'
     character(LEN=18) :: COARSE_DUST_NUMBER     = 'coarse_dust_number'
     character(LEN=16) :: ACTIVE_INSOL_ICE       = 'active_insol_ice'
     character(LEN=14) :: ACTIVE_SOL_ICE         = 'active_sol_ice'
     character(LEN=19) :: ACTIVE_INSOL_LIQUID    = 'active_insol_liquid'
     character(LEN=16) :: ACCUM_INSOL_MASS       = 'accum_insol_mass'
     character(LEN=18) :: ACCUM_INSOL_NUMBER     = 'accum_insol_number'
     character(LEN=17) :: ACTIVE_SOL_NUMBER      = 'active_sol_number'
     character(LEN=19) :: ACTIVE_INSOL_NUMBER    = 'active_insol_number'
  end type standard_q_names_type

  type(standard_q_names_type) :: standard_q_names

  public q_metadata_type, get_q_index, get_indices_descriptor, get_max_number_q_indices, &
       get_number_active_q_indices, set_q_index, standard_q_names


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
  integer function get_q_index(name, assigning_component)
    character(*), intent(in) :: name
    character(*), optional :: assigning_component

    integer :: iname, i_unused

    i_unused=0
    
    do iname=n_maxqs, 1, -1
      if (trim(name) == trim(q_register(iname)%name)) then
        ! Here we find the variable has already been added
        get_q_index=iname
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
    get_q_index=i_unused
    q_register(get_q_index)%name = adjustr(trim(name))
    q_register(get_q_index)%l_used = .true.

    if (present(assigning_component)) then
      call log_master_log(LOG_INFO, 'q variable #'//trim(conv_to_string(get_q_index))//' is assigned to '//trim(name)//&
           '. Assigned from component: '//trim(assigning_component))
    else
      call log_master_log(LOG_INFO, 'q variable #'//trim(conv_to_string(get_q_index))//' is assigned to '//trim(name))
    end if

  end function get_q_index
end module q_indices_mod
