!> Add random noise into the fields
module randomnoise_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : FORWARD_STEPPING, model_state_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
     options_get_logical_array, options_get_real_array, options_get_string_array 
  use interpolation_mod, only: piecewise_linear_1d
  use q_indices_mod, only: q_indices_add

implicit none

#ifndef TEST_MODE
  private
#endif
    
  integer, parameter :: MAX_SIZE_SEED_ARRAY=256, I_SEED=7, ISD=1

  public randomnoise_get_descriptor
contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function randomnoise_get_descriptor()
    randomnoise_get_descriptor%name="randomnoise"
    randomnoise_get_descriptor%version=0.1
    randomnoise_get_descriptor%initialisation=>initialisation_callback
  end function randomnoise_get_descriptor

  !> The initialisation callback sets up the buoyancy coefficient
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer, dimension(MAX_SIZE_SEED_ARRAY) :: iranseed
    real(kind=DEFAULT_PRECISION), dimension(current_state%global_grid%size(X_INDEX), current_state%global_grid%size(Y_INDEX), &
       current_state%global_grid%size(Z_INDEX)) :: randarr

    integer :: nq_rand ! The number of q fields to add noise to
    integer :: nzq     ! The number of input levels for noise
    integer :: i,j,k,n ! loop counters
    integer :: iq      ! temporary q varible index

    integer, parameter :: MAXQIN=10, MAXZIN=10
    integer, parameter :: UNSET_REAL=-999.0
 
    real(kind=DEFAULT_PRECISION) :: f_rand_pl_q(MAXZIN,MAXQIN)=UNSET_REAL   ! Random Noise node amplitude for q variables
    real(kind=DEFAULT_PRECISION) :: z_rand_pl_q(MAXZIN)=UNSET_REAL     ! Random Noise node height values for q variables
    real(kind=DEFAULT_PRECISION) :: f_rand_pl_theta(MAXZIN)=UNSET_REAL ! Random Noise node amplitude for theta variable
    real(kind=DEFAULT_PRECISION) :: z_rand_pl_theta(MAXZIN)=UNSET_REAL ! Random Noise node height values for theta variable

    logical :: l_rand_pl_theta ! if .true. then random noise added to theta field
    logical :: l_rand_pl_q     ! if .true. then random noise added to q fields

    character(len=STRING_LENGTH) :: names_rand_pl_q(MAXZIN)='unset'  ! names of q variables to add random noise to
    
    real(kind=DEFAULT_PRECISION), allocatable :: f_rand_pl_q_tmp(:) !temporary 1D storage of random noise for q field
    real(kind=DEFAULT_PRECISION), allocatable :: zgrid(:)  ! z grid to use in interpolation

    integer :: number_in ! Number of entries input into array

    allocate(zgrid(current_state%local_grid%local_domain_end_index(Z_INDEX)))

    l_rand_pl_theta=options_get_logical(current_state%options_database, "l_rand_pl_theta")
    l_rand_pl_q=options_get_logical(current_state%options_database, "l_rand_pl_q")
    if (l_rand_pl_q)call options_get_string_array(current_state%options_database, "names_rand_pl_q", names_rand_pl_q)

    iranseed(1:ISD)=I_SEED

    if (l_rand_pl_theta)then
      ! Get random numbers
      call random_seed(get=iranseed)
      call random_number(randarr)

      ! Get amplitude profiles
      call options_get_real_array(current_state%options_database, "z_rand_pl_theta", z_rand_pl_theta)
      call options_get_real_array(current_state%options_database, "f_rand_pl_theta", f_rand_pl_theta)
      do i=1,size(z_rand_pl_theta)
        if (z_rand_pl_theta(i) == UNSET_REAL) exit 
      end do
      number_in=i-1
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      call piecewise_linear_1d(z_rand_pl_theta(1:number_in), f_rand_pl_theta(1:number_in), zgrid, &
         current_state%global_grid%configuration%vertical%theta_rand)
      do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
        do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
          do k=2, current_state%local_grid%local_domain_end_index(Z_INDEX)          
            current_state%th%data(k,j,i) = current_state%th%data(k,j,i) + &
               current_state%global_grid%configuration%vertical%theta_rand(k) * 2.0 * (randarr( &
               i-current_state%local_grid%local_domain_start_index(X_INDEX)+current_state%local_grid%start(X_INDEX), &
               j-current_state%local_grid%local_domain_start_index(Y_INDEX)+current_state%local_grid%start(Y_INDEX), &
               k)-0.5)
          end do
        end do
      end do

    end if

    if (l_rand_pl_q)then
      do i=1,size(names_rand_pl_q)
        if (trim(names_rand_pl_q(i)) == trim('unset')) exit 
      end do
      nq_rand=i-1
      call options_get_real_array(current_state%options_database, "z_rand_pl_q", z_rand_pl_q)
      do i=1,size(z_rand_pl_q)
        if (z_rand_pl_q(i) == UNSET_REAL) exit 
      end do
      nzq=i-1
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      allocate(f_rand_pl_q_tmp(nq_rand*nzq))
      call options_get_real_array(current_state%options_database, "f_rand_pl_q", f_rand_pl_q_tmp)
      f_rand_pl_q(1:nzq, 1:nq_rand)=reshape(f_rand_pl_q_tmp, (/nzq, nq_rand/))
      do n=1,nq_rand
        ! Get random numbers
        call random_seed(get=iranseed)
        call random_number(randarr)
        
        iq=q_indices_add(trim(names_rand_pl_q(n)), 'random noise')
        zgrid=current_state%global_grid%configuration%vertical%zn(:)
        call piecewise_linear_1d(z_rand_pl_q(1:number_in), f_rand_pl_q(1:nzq,n), zgrid, &
           current_state%global_grid%configuration%vertical%q_rand(:,iq))
        do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
          do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
            do k=2, current_state%local_grid%local_domain_end_index(Z_INDEX)          
              current_state%q(iq)%data(k,j,i) = current_state%q(iq)%data(k,j,i) + &
                 current_state%global_grid%configuration%vertical%q_rand(k,iq) * 2.0 * (randarr( &
                 i-current_state%local_grid%local_domain_start_index(X_INDEX)+current_state%local_grid%start(X_INDEX), &
                 j-current_state%local_grid%local_domain_start_index(Y_INDEX)+current_state%local_grid%start(Y_INDEX), &
                 k)-0.5)
            end do
          end do
        end do
      end do
    end if

    deallocate(zgrid)
    
  end subroutine initialisation_callback
end module randomnoise_mod