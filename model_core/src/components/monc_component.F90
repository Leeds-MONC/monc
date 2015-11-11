!> Interfaces and types that MONC components must specify
!!
!! This module should be used by all MONC components to specify a description
!! of themselves (for registration) and it also specifies the interfaces of
!! each callback procedure
module monc_component_mod
    use state_mod, only : model_state_type
    use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
    implicit none

#ifndef TEST_MODE
    private
#endif

    integer, parameter :: COMPONENT_SCALAR_FIELD_TYPE = 1, COMPONENT_ARRAY_FIELD_TYPE=2
    integer, parameter :: COMPONENT_INTEGER_DATA_TYPE = 1, COMPONENT_DOUBLE_DATA_TYPE=5

    !> Index of each priority value in the descriptor
    integer, public, parameter :: INIT_PRIORITY_INDEX=1, TIMESTEP_PRIORITY_INDEX=2, FINALISATION_PRIORITY_INDEX=5 

    !> Wrapper type for the value returned for a published field from a component
    type, public :: component_field_value_type
       real(kind=DEFAULT_PRECISION), dimension(:,:,:,:), allocatable :: real_4d_array
       real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: real_3d_array
       real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: real_2d_array
       real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: real_1d_array
       real(kind=DEFAULT_PRECISION) :: scalar_real
       integer :: scalar_int
    end type component_field_value_type

    type, public :: component_field_information_type
       integer :: field_type, data_type, number_dimensions, dimension_sizes(4)
       logical :: enabled
    end type component_field_information_type

    !> Description of a component
    !!
    !! Provided to the registry to register a component. It is also used by the registry to store information
    !! about components. Each component should initialise the name, version and appropriate procedure pointers
    !! that their component implements. Leave any other procedure pointers null as this indicated not to install
    !! a call back hook for that specific stage.
    type, public :: component_descriptor_type
      procedure(component_initialisation), nopass, pointer :: initialisation => null()  !< Initialisation callback pointer
      procedure(component_timestep), nopass, pointer :: timestep => null()              !< Timestep callback pointer
      procedure(component_finalisation), nopass, pointer :: finalisation => null()      !< Finalisation callback pointer      
      procedure(component_get_field_value), nopass, pointer :: field_value_retrieval => null()
      procedure(component_get_field_information), nopass, pointer :: field_information_retrieval => null()
      character(len=STRING_LENGTH), dimension(:), pointer :: published_fields => null()
      character(len=STRING_LENGTH) :: name !< Component name
      real(kind=DEFAULT_PRECISION) :: version !< Component version number
    end type component_descriptor_type

    !> Interface defining the signature of each callback hook that a component may specify.
    !!
    !! Each hook accepts the current model state (which may be modified.)
    abstract interface
      function component_get_description()
        import :: component_descriptor_type
        type(component_descriptor_type) :: component_get_description
      end function component_get_description
      !> Component initialisation callback hook signature
      !! @param simulationState The current model state which may be modified.
      subroutine component_initialisation(current_state)
        import :: model_state_type
        type(model_state_type), target, intent(inout) :: current_state
      end subroutine component_initialisation

      !> Component timestep callback hook signature
      !! @param simulationState The current model state which may be modified.
      subroutine component_timestep(current_state)
        import :: model_state_type
        type(model_state_type), target, intent(inout) :: current_state
      end subroutine component_timestep

      !> Component finalisation callback hook signature
      !! @param simulationState The current model state which may be modified.
      subroutine component_finalisation(current_state)
        import :: model_state_type
        type(model_state_type), target, intent(inout) :: current_state
      end subroutine component_finalisation

      subroutine component_get_field_information(current_state, name, field_sizing)
        import :: component_field_information_type, model_state_type
        type(model_state_type), target, intent(inout) :: current_state
        character(len=*), intent(in) :: name
        type(component_field_information_type), intent(out) :: field_sizing
      end subroutine component_get_field_information      

      !> Retrieves a specific published field value from a component
      !! @param name The name of the published field to access
      !! @param field_value The fields value wrapper
      subroutine component_get_field_value(current_state, name, field_value)
        import :: component_field_value_type, model_state_type
        type(model_state_type), target, intent(inout) :: current_state
        character(len=*), intent(in) :: name
        type(component_field_value_type), intent(out) :: field_value
      end subroutine component_get_field_value      
    end interface

    public COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
         component_initialisation, component_timestep, component_finalisation, component_get_description
end module monc_component_mod
