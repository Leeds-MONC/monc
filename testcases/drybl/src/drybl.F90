!> Dry boundary layer test case, which represents test case 1 in the LEM
module drybl_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : PRESCRIBED_SURFACE_FLUXES, model_state_type
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer, parameter :: MAX_SIZE_SEED_ARRAY=256, I_SEED=7, ISD=1
  logical, parameter :: USE_PSEUDO_RANDOM = .false.

  public drybl_get_descriptor
contains

  type(component_descriptor_type) function drybl_get_descriptor()
    drybl_get_descriptor%name="drybl"
    drybl_get_descriptor%version=0.1
    drybl_get_descriptor%initialisation=>initialisation_callback
  end function drybl_get_descriptor

  !> Sets up the field values for this test case
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    ! Only generate the dry boundary layer if starting from fresh, otherwise this is loaded in
    if (current_state%timestep == 1) call generate_drybl(current_state)
  end subroutine initialisation_callback

  !> Generates the dry boundary layer
  !! @param current_state The current model state
  subroutine generate_drybl(current_state)
    type(model_state_type), intent(inout), target :: current_state
    
    integer, dimension(MAX_SIZE_SEED_ARRAY) :: iranseed
    integer :: i, j, k
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(X_INDEX), current_state%local_grid%size(Y_INDEX), &
         current_state%local_grid%size(Z_INDEX)) :: randarr

    current_state%type_of_surface_boundary_conditions=PRESCRIBED_SURFACE_FLUXES
    
    if (USE_PSEUDO_RANDOM) then
      call fill_pseudo_array(randarr)
    else
      iranseed(1:ISD)=I_SEED
      call random_seed(get=iranseed)
      call random_number(randarr)
    end if

    do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
      do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
        do k=2, current_state%local_grid%local_domain_end_index(Z_INDEX)
#ifdef U_ACTIVE
          current_state%u%data(k,j,i)=current_state%surface_geostrophic_wind_x
#endif
#ifdef V_ACTIVE
          current_state%v%data(k,j,i)=current_state%surface_geostrophic_wind_y
#endif
#ifdef W_ACTIVE
          current_state%w%data(k,j,i) = 0.0_DEFAULT_PRECISION
          if (current_state%global_grid%configuration%vertical%z(k) .lt. 500.0_DEFAULT_PRECISION) then            
            current_state%w%data(k,j,i) = current_state%w%data(k,j,i)+(randarr(&
                 (i-current_state%local_grid%local_domain_start_index(X_INDEX))+1, &
                 (j-current_state%local_grid%local_domain_start_index(Y_INDEX))+1, k)-0.5_DEFAULT_PRECISION)*&
                 (1.0_DEFAULT_PRECISION- current_state%global_grid%configuration%vertical%z(k)/500.0_DEFAULT_PRECISION)
          end if
#endif
        end do
#ifdef W_ACTIVE
        current_state%w%data(current_state%local_grid%local_domain_end_index(Z_INDEX),j,i)=0.0_DEFAULT_PRECISION
        current_state%w%data(1,j,i)=0.0_DEFAULT_PRECISION
#endif
        if (current_state%use_viscosity_and_diffusion) then
#ifdef U_ACTIVE
          current_state%u%data(1,j,i)=-current_state%u%data(2,j,i)
#endif
#ifdef V_ACTIVE
          current_state%v%data(1,j,i)=-current_state%v%data(2,j,i)
#endif
        else
#ifdef U_ACTIVE
          current_state%u%data(1,j,i)=current_state%u%data(2,j,i)
#endif
#ifdef V_ACTIVE
          current_state%v%data(1,j,i)=current_state%v%data(2,j,i)
#endif
        end if
      end do      
    end do
#ifdef U_ACTIVE
    current_state%zu%data=current_state%u%data
#endif
#ifdef V_ACTIVE
    current_state%zv%data=current_state%v%data
#endif
#ifdef W_ACTIVE
    current_state%zw%data=current_state%w%data
#endif
  end subroutine generate_drybl  

  !> Generates a pseudo random array, which allows us to do comparisons between different compilers and architectures
  !! @param data The data to fill with the pseudo random values
  subroutine fill_pseudo_array(data)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: data

    integer :: i, j, k, m
    m=size(data,3)+(size(data,2)*5)+(size(data,1)*10)
    do i=1, size(data,3)
      do j=1, size(data,2)
        do k=1, size(data,1)          
          data(k,j,i)=(i+j*5+k*10) / (2.0_DEFAULT_PRECISION*m)          
        end do
      end do
    end do
  end subroutine fill_pseudo_array
end module drybl_mod
