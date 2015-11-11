!> This is the iterative pressure solver and uses a Jacobi preconditioned BiCGStab which we implement here
module iterativesolver_mod
  use monc_component_mod, only : component_descriptor_type
  use collections_mod, only : map_type
  use optionsdatabase_mod, only : options_get_real, options_get_integer, options_get_logical
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX, local_grid_type, grid_configuration_type
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use logging_mod, only : LOG_WARN, LOG_DEBUG, log_log, log_get_logging_level
  use conversions_mod, only : conv_to_string
  use communication_types_mod, only : halo_communication_type, neighbour_description_type, field_data_wrapper_type
  use halo_communication_mod, only : copy_buffer_to_field, copy_field_to_buffer, perform_local_data_copy_for_field, &
       init_halo_communication, finalise_halo_communication, initiate_nonblocking_halo_swap, complete_nonblocking_halo_swap, &
       blocking_halo_swap, get_single_field_per_halo_cell
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, log_master_log
  use mpi, only : MPI_MAX, MPI_SUM, MPI_COMM_WORLD, MPI_REQUEST_NULL, MPI_STATUSES_IGNORE
  implicit none

#ifndef TEST_MODE
  private
#endif

  !> A helper type to abstract the concrete details of the matrix
  type matrix_type
     real(kind=DEFAULT_PRECISION) :: n, s, e, w
     real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: u, d, p, lu_d, lu_u, vol
  end type matrix_type

  real(kind=DEFAULT_PRECISION) :: tolerance, relaxation     !< Solving tollerance
  integer :: max_iterations, &     !< Maximum number of BiCGStab iterations
       preconditioner_iterations     !< Number of preconditioner iterations to perform per call
  logical :: symm_prob
  
  real(kind=DEFAULT_PRECISION), parameter :: TINY = 1.0e-16 !< Minimum residual - if we go below this then something has gone wrong

  type(halo_communication_type), save :: halo_swap_state      !< The halo swap state as initialised by that module
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: psource, prev_p      !< Passed to BiCGStab as the RHS
  logical :: first_run=.true.
  type(matrix_type) :: A

  public iterativesolver_get_descriptor
contains

  !> Descriptor of the iterative solver component used by the registry
  !! @returns The iterative solver component descriptor
  type(component_descriptor_type) function iterativesolver_get_descriptor()
    iterativesolver_get_descriptor%name="iterativesolver"
    iterativesolver_get_descriptor%version=0.1
    iterativesolver_get_descriptor%initialisation=>initialisation_callback
    iterativesolver_get_descriptor%timestep=>timestep_callback
    iterativesolver_get_descriptor%finalisation=>finalisation_callback
  end function iterativesolver_get_descriptor

  !> Initialisation callback hook which will set up the halo swapping state and allocate some data
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (.not. is_component_enabled(current_state%options_database, "diverr")) then
      call log_master_log(LOG_ERROR, "The iterative solver component requires the diverr component to be enabled")
    end if

    tolerance=options_get_real(current_state%options_database, "tolerance")
    max_iterations=options_get_integer(current_state%options_database, "max_iterations")
    preconditioner_iterations=options_get_integer(current_state%options_database, "preconditioner_iterations")
    symm_prob=options_get_logical(current_state%options_database, "symm_prob")

    call init_halo_communication(current_state, get_single_field_per_halo_cell, halo_swap_state, 1, .false.)

    allocate(psource(current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2, &
         current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2, &
         current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2),&
         prev_p(current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2, &
         current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2, &
         current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2))

    A=create_problem_matrix(current_state%local_grid%size(Z_INDEX))
    call set_matrix_for_poisson(current_state%global_grid%configuration, A, current_state%local_grid%size(Z_INDEX))
  end subroutine initialisation_callback

  !> Timestep callback, this ignores all but the last column where it calls the solver
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: i_strt, i_end, j_strt, j_end, k_end

    i_strt = current_state%local_grid%local_domain_start_index(X_INDEX)
    i_end  = current_state%local_grid%local_domain_end_index(X_INDEX)
    j_strt = current_state%local_grid%local_domain_start_index(Y_INDEX)
    j_end  = current_state%local_grid%local_domain_end_index(Y_INDEX)
    k_end  = current_state%local_grid%size(Z_INDEX)

    call complete_psrce_calculation(current_state, current_state%local_grid%halo_size(Y_INDEX), &
         current_state%local_grid%halo_size(X_INDEX))
     
    call initiate_nonblocking_halo_swap(current_state, halo_swap_state, copy_p_to_halo_buffer)
    call deduce_global_divmax(current_state)    
    call complete_nonblocking_halo_swap(current_state, halo_swap_state, perform_local_data_copy_for_p, copy_halo_buffer_to_p)

    psource=current_state%p%data
    if (first_run) then
      ! If first timestep then initial guess is zero
      current_state%p%data=0.0_DEFAULT_PRECISION
      first_run=.false.
    else
      ! Initial guess is set to previous timesteps p
      current_state%p%data=prev_p
    end if
    
    if (symm_prob) then
      call cg_solver(current_state, A, current_state%p%data, psource, i_strt, i_end, j_strt, j_end, k_end)
    else
      call bicgstab(current_state, A, current_state%p%data, psource, i_strt, i_end, j_strt, j_end, k_end)
    end if

    prev_p=current_state%p%data
  end subroutine timestep_callback

  !> Called as MONC is shutting down and frees the halo swap state and deallocates local data
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call finalise_halo_communication(halo_swap_state)
    deallocate(psource, prev_p, A%u, A%d, A%p, A%lu_u, A%lu_d, A%vol)
  end subroutine finalisation_callback  

  !> Performs the BiCGStab KS method
  !! @param current_state The current model state
  !! @param A The matrix
  !! @param x The solution
  !! @param b The RHS
  subroutine bicgstab(current_state, A, x, b, i_strt, i_end, j_strt, j_end, k_end)
    type(model_state_type), target, intent(inout) :: current_state
    type(matrix_type), intent(inout) :: A    
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: x
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: b
    integer, intent(in) :: i_strt, i_end, j_strt, j_end, k_end

    integer :: it, i, j, k
    real(kind=DEFAULT_PRECISION) :: sc_err, alf, omg, nrm, my_rho, bet, tt, ts, ss, err, init_err, inner_prod_results(3)
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX) + &
         current_state%local_grid%halo_size(Z_INDEX) * 2, &
         current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2, &
         current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2) :: Ax, r, cr, pp, v, t, s, cs

    ! Calculate scale factor for error
    sc_err = sqrt(inner_prod(current_state, b, b, i_strt, i_end, j_strt, j_end, k_end))
    sc_err = max(sc_err, 0.0001_DEFAULT_PRECISION)

    ! Calculate initial residual
    call calc_Ax(current_state, A, x, Ax)

    do i = i_strt, i_end
      do j = j_strt, j_end
        do k = 2, k_end
          r(k,j,i) = b(k,j,i) - Ax(k,j,i)
          cr(k,j,i) = r(k,j,i)
        end do
      end do
    end do

    my_rho = inner_prod(current_state, r, r, i_strt, i_end, j_strt, j_end, k_end)
    err = sqrt(my_rho)/sc_err
    init_err = err

    alf = 1.0_DEFAULT_PRECISION
    omg = 1.0_DEFAULT_PRECISION
    nrm = 1.0_DEFAULT_PRECISION

    if (err .ge. tolerance) then
      do it=1, max_iterations
        if (it > 1) my_rho = inner_prod(current_state, r, cr, i_strt, i_end, j_strt, j_end, k_end)
        bet = (my_rho/nrm) * (alf/omg)
        if (it == 1) then
          call precond(current_state, A, pp, r, preconditioner_iterations)
        else
          do i = i_strt, i_end
            do j = j_strt, j_end
              do k = 2, k_end
                t(k,j,i) = r(k,j,i) - bet*omg*v(k,j,i)
              end do
            end do
          end do
          call precond(current_state, A, s, t, preconditioner_iterations)
          do i = i_strt, i_end
            do j = j_strt, j_end
              do k = 2, k_end
                pp(k,j,i) = s(k,j,i) + bet*pp(k,j,i)
              end do
            end do
          end do
        end if       
        call calc_Ax(current_state, A, pp, v)
        nrm = inner_prod(current_state, cr, v, i_strt, i_end, j_strt, j_end, k_end)
        alf = my_rho / nrm

        do i = i_strt, i_end
          do j = j_strt, j_end
            do k = 2, k_end
              s(k,j,i) = r(k,j,i) - alf*v(k,j,i)
            end do
          end do
        end do

        call precond(current_state, A, cs, s, preconditioner_iterations)
        call calc_Ax(current_state, A, cs, t)

        inner_prod_results=inner_prod_three_way(current_state, t, s, i_strt, i_end, j_strt, j_end, k_end)
        tt = inner_prod_results(1)
        ts = inner_prod_results(2)
        ss = inner_prod_results(3)
        omg = ts/tt
        x = x + alf*pp + omg*cs
        do i = i_strt, i_end
          do j = j_strt, j_end
            do k = 2, k_end
              r(k,j,i) = s(k,j,i) - omg*t(k,j,i)
            end do
          end do
        end do
        nrm = my_rho

        if (abs(omg) < TINY) then
          call log_log(LOG_WARN, "Convergence problem, omega="//conv_to_string(omg))
        endif

        err = sqrt(ss - 2*omg*ts + omg**2 *tt)/sc_err        
        if (err < tolerance) exit        
      end do
    end if
    
    if (err > tolerance) then
      call log_log(LOG_WARN, "Convergence failed, RNorm="//conv_to_string(err, exponent=.true.))
    else if (current_state%parallel%my_rank==0 .and. log_get_logging_level() .eq. LOG_DEBUG) then
      call log_log(LOG_DEBUG, "Converged in "//trim(conv_to_string(it))//" iterations with RNorm="//&
           trim(conv_to_string(err, 5, .true.))//" initial norm="//trim(conv_to_string(init_err, 5, .true.)))
    end if
  end subroutine bicgstab

    !> Performs the preconditioned conjugate gradient method
  !! @param current_state The current model state
  !! @param A The matrix
  !! @param x The solution
  !! @param b The RHS
  subroutine cg_solver(current_state, A, x, b, i_strt, i_end, j_strt, j_end, k_end)
    type(model_state_type), target, intent(inout) :: current_state
    type(matrix_type), intent(inout) :: A
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: x
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: b
    integer, intent(in) :: i_strt, i_end, j_strt, j_end, k_end

    integer :: it, k, i, j
    real(kind=DEFAULT_PRECISION) :: sc_err, alf, bet, err, init_err, rho
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX) + &
         current_state%local_grid%halo_size(Z_INDEX) * 2, &
         current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2, &
         current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2) :: Ax, r, z, p

    ! first rescale RHS for symmetry (this could be done when p_source is calculated
    do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
      do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)

        r(1,j,i) = 0.0_DEFAULT_PRECISION
        do k=2,current_state%local_grid%size(Z_INDEX)
           r(k,j,i) = b(k,j,i) * A%vol(k)
        end do
      end do
    end do
 
    ! Calculate scale factor for error

    call calc_Ax(current_state, A, x, Ax)

    sc_err = sqrt(inner_prod(current_state, r, r, i_strt, i_end, j_strt, j_end, k_end))
    sc_err = max(sc_err, 0.0001_DEFAULT_PRECISION)
    r = r - Ax
    init_err = sqrt(inner_prod(current_state, r, r, i_strt, i_end, j_strt, j_end, k_end))/sc_err

    do it=1, max_iterations
       if( it == 1 ) then
          call precond(current_state, A, p, r, preconditioner_iterations)
          rho = inner_prod(current_state, p, r, i_strt, i_end, j_strt, j_end, k_end)
          alf = rho
       else
          call precond(current_state, A, z, r, preconditioner_iterations)
          alf = inner_prod(current_state, z, r, i_strt, i_end, j_strt, j_end, k_end)
          bet = alf/rho
          rho = alf
          p   = z + bet*p
       end if

       call calc_Ax(current_state, A, p, Ax)
       alf = alf/inner_prod(current_state, p, Ax, i_strt, i_end, j_strt, j_end, k_end)
       x   = x + alf*p
       r   = r - alf*Ax

       err = sqrt(inner_prod(current_state, r, r, i_strt, i_end, j_strt, j_end, k_end))/sc_err
       if (err < tolerance) exit
    end do

    if( current_state%parallel%my_rank == 0 ) print*,it, err, init_err

    if (err > tolerance) then
      call log_log(LOG_WARN, "Convergence failed, RNorm="//conv_to_string(err, exponent=.true.))
    else if (current_state%parallel%my_rank==0 .and. log_get_logging_level() .eq. LOG_DEBUG) then
      call log_log(LOG_DEBUG, "Converged in "//trim(conv_to_string(it))//" iterations with RNorm="//&
           trim(conv_to_string(err, 5, .true.))//" initial norm="//trim(conv_to_string(init_err, 5, .true.)))
    end if
   end subroutine cg_solver

  !> Jacobi preconditioner
  !! @param current_state The current model state
  !! @param A The matrix
  !! @param s Written into as result of preconditioning
  !! @param r Input values to preconditioner
  !! @param preits Number of iterations of the preconditioner to perform per call
  subroutine precond(current_state, A, s, r, preits)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: r
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: s
    integer, intent(in) :: preits
    type(matrix_type), intent(inout) :: A

    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX) + &
         current_state%local_grid%halo_size(Z_INDEX) * 2, current_state%local_grid%size(Y_INDEX) + &
         current_state%local_grid%halo_size(Y_INDEX) * 2, current_state%local_grid%size(X_INDEX) + &
         current_state%local_grid%halo_size(X_INDEX) * 2) :: t
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: s_k
    integer :: it, i, j, k

    if (preits .lt. 0) then
      s=r
      return
    end if    

    do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
      do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
        s(1,j,i) = 0.0_DEFAULT_PRECISION
        k=2 
        s(k,j,i)=r(k,j,i)*A%LU_d(k)
        do k=3,current_state%local_grid%size(Z_INDEX)
          s(k,j,i)=(r(k,j,i) - A%d(k)*s(k-1,j,i))*A%lu_d(k)
        end do
        do k=current_state%local_grid%size(Z_INDEX)-1, 2, -1
          s(k,j,i)=s(k,j,i) - A%lu_u(k)*s(k+1,j,i)
        end do        
      end do
    end do

    do it=1, preits
      call calc_Ax(current_state, A, s, t)
      do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
      do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
          k=2 
          s_k(k)=(r(k,j,i) - t(k,j,i))*A%lu_d(k) 
          do k=3,current_state%local_grid%size(Z_INDEX) 
            s_k(k)=(r(k,j,i) - t(k,j,i) - A%d(k)*s_k(k-1))*A%lu_d(k) 
          end do
          k=current_state%local_grid%size(Z_INDEX)
          s(k,j,i)=s(k,j,i)+s_k(k) 
          do k=current_state%local_grid%size(Z_INDEX)-1, 2, -1 
            s_k(k)=s_k(k) - A%lu_u(k)*s_k(k+1) 
            s(k,j,i)=s(k,j,i) + relaxation*s_k(k) 
          end do
        end do
      end do
    end do
  end subroutine precond

  !> Calculates A * x
  !! @param current_state The current model state
  !! @param A The matrix
  !! @param x Vector to multiply with
  !! @param Ax Result of A*x
  subroutine calc_Ax(current_state, A, x, Ax)
    type(model_state_type), target, intent(inout) :: current_state
    type(matrix_type), intent(in) :: A
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), target, intent(inout) :: x, Ax

    integer :: i, k, j, n, istart, iend, jstart, jend
    type(field_data_wrapper_type) :: source_data   

    source_data%data=>x

    call initiate_nonblocking_halo_swap(current_state, halo_swap_state, &
         copy_calc_Ax_to_halo_buffer, source_data=(/source_data/))

    Ax(1,:,:) = 0.0_DEFAULT_PRECISION
    if (symm_prob) then
      do n=1, 5
        if (n==1) then
          istart=current_state%local_grid%local_domain_start_index(X_INDEX)+1
          iend=current_state%local_grid%local_domain_end_index(X_INDEX)-1
          jstart=current_state%local_grid%local_domain_start_index(Y_INDEX)+1
          jend=current_state%local_grid%local_domain_end_index(Y_INDEX)-1
        else if (n==2) then
          istart=current_state%local_grid%local_domain_start_index(X_INDEX)
          iend=current_state%local_grid%local_domain_start_index(X_INDEX)
        else if (n==3) then
          istart=current_state%local_grid%local_domain_end_index(X_INDEX)
          iend=current_state%local_grid%local_domain_end_index(X_INDEX)
        else if (n==4) then
          jstart=current_state%local_grid%local_domain_start_index(Y_INDEX)
          jend=current_state%local_grid%local_domain_start_index(Y_INDEX)
          istart=current_state%local_grid%local_domain_start_index(X_INDEX)
          iend=current_state%local_grid%local_domain_end_index(X_INDEX)
        else if (n==5) then
          jstart=current_state%local_grid%local_domain_end_index(Y_INDEX)
          jend=current_state%local_grid%local_domain_end_index(Y_INDEX)
        end if
        do i=istart, iend
          do j=jstart, jend
            k=2
            Ax(k,j,i)=A%vol(k)*(A%n*(x(k,j,i+1)+x(k,j,i-1))+A%e*(x(k,j+1,i)+x(k,j-1,i)))+ A%u(k)*x(k+1,j,i)+A%p(k)*x(k,j,i)
            do k=3,current_state%local_grid%size(Z_INDEX)-1
              Ax(k,j,i)=A%vol(k)*(A%n*(x(k,j,i+1)+x(k,j,i-1))+A%e*(x(k,j+1,i)+x(k,j-1,i)))+&
                   A%u(k)*x(k+1,j,i)+A%d(k)*x(k-1,j,i)+A%p(k)*x(k,j,i)
            end do
            k=current_state%local_grid%size(Z_INDEX)
            Ax(k,j,i) = A%vol(k)*(A%n*(x(k,j,i+1)+x(k,j,i-1))+A%e*(x(k,j+1,i)+x(k,j-1,i)))+ A%d(k)*x(k-1,j,i)+A%p(k)*x(k,j,i)
          end do
        end do
        if (n==1) then
          call complete_nonblocking_halo_swap(current_state, halo_swap_state, perform_local_data_copy_for_calc_Ax, &
               copy_halo_buffer_to_calc_Ax, source_data=(/source_data/))
        end if
      end do
    else
      do n=1, 5
        if (n==1) then
          istart=current_state%local_grid%local_domain_start_index(X_INDEX)+1
          iend=current_state%local_grid%local_domain_end_index(X_INDEX)-1
          jstart=current_state%local_grid%local_domain_start_index(Y_INDEX)+1
          jend=current_state%local_grid%local_domain_end_index(Y_INDEX)-1
        else if (n==2) then
          istart=current_state%local_grid%local_domain_start_index(X_INDEX)
          iend=current_state%local_grid%local_domain_start_index(X_INDEX)
        else if (n==3) then
          istart=current_state%local_grid%local_domain_end_index(X_INDEX)
          iend=current_state%local_grid%local_domain_end_index(X_INDEX)
        else if (n==4) then
          jstart=current_state%local_grid%local_domain_start_index(Y_INDEX)
          jend=current_state%local_grid%local_domain_start_index(Y_INDEX)
          istart=current_state%local_grid%local_domain_start_index(X_INDEX)
          iend=current_state%local_grid%local_domain_end_index(X_INDEX)
        else if (n==5) then
          jstart=current_state%local_grid%local_domain_end_index(Y_INDEX)
          jend=current_state%local_grid%local_domain_end_index(Y_INDEX)
        end if
        do i=istart, iend
          do j=jstart, jend
            k=2
            Ax(k,j,i)=A%n*(x(k,j,i+1)+x(k,j,i-1))+A%e*(x(k,j+1,i)+x(k,j-1,i))+ A%u(k)*x(k+1,j,i)+A%p(k)*x(k,j,i)
            do k=3,current_state%local_grid%size(Z_INDEX)-1
              Ax(k,j,i)=A%n*(x(k,j,i+1)+x(k,j,i-1))+A%e*(x(k,j+1,i)+x(k,j-1,i))+&
                   A%u(k)*x(k+1,j,i)+A%d(k)*x(k-1,j,i)+A%p(k)*x(k,j,i)
            end do
            k=current_state%local_grid%size(Z_INDEX)
            Ax(k,j,i) = A%n*(x(k,j,i+1)+x(k,j,i-1))+A%e*(x(k,j+1,i)+x(k,j-1,i))+ A%d(k)*x(k-1,j,i)+A%p(k)*x(k,j,i)
          end do
        end do
        if (n==1) then
          call complete_nonblocking_halo_swap(current_state, halo_swap_state, perform_local_data_copy_for_calc_Ax, &
               copy_halo_buffer_to_calc_Ax, source_data=(/source_data/))
        end if
      end do
    endif
  end subroutine calc_Ax

  !> Returns the global inner product of two vectors, ignoring the halo cells
  !! @param current_state The current model state
  !! @param x First vector
  !! @praam y Second vector
  !! @returns Global inner product of the two input vectors
  real(kind=DEFAULT_PRECISION) function inner_prod(current_state, x, y, i_strt, i_end, j_strt, j_end, k_end)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: x, y
    integer, intent(in) :: i_strt, i_end, j_strt, j_end, k_end

    real(kind=DEFAULT_PRECISION) :: local_sum, global_sum
    integer :: ierr, i, j, k

    local_sum=0.0_DEFAULT_PRECISION

     do i=i_strt, i_end
      do j=j_strt, j_end
        do k=2, k_end
          local_sum=local_sum+x(k,j,i)*y(k,j,i)
        end do        
      end do      
    end do

    call mpi_allreduce(local_sum, global_sum, 1, PRECISION_TYPE, MPI_SUM, current_state%parallel%monc_communicator, ierr)
    inner_prod=global_sum
  end function inner_prod

  !> Returns the global inner product of a pair of vectors, ignoring the halo cells for three separate pairs. This call
  !! is for optimisation to bunch up the comms in a BiCGStab solver per iteration
  !! @param current_state The current model state
  !! @param x First vector
  !! @praam y Second vector
  !! @returns Global inner product of the two input vectors
  function inner_prod_three_way(current_state, t, s, i_strt, i_end, j_strt, j_end, k_end)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: i_strt, i_end, j_strt, j_end, k_end
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: t, s
    real(kind=DEFAULT_PRECISION), dimension(3) :: inner_prod_three_way

    real(kind=DEFAULT_PRECISION), dimension(3) :: local_sum, global_sum
    integer :: ierr, i, j, k

    local_sum(1)=0.0_DEFAULT_PRECISION
    local_sum(2)=0.0_DEFAULT_PRECISION
    local_sum(3)=0.0_DEFAULT_PRECISION

    do i=i_strt, i_end
      do j=j_strt, j_end
        do k=2, k_end
          local_sum(1)=local_sum(1)+t(k,j,i)*t(k,j,i)
          local_sum(2)=local_sum(2)+t(k,j,i)*s(k,j,i)
          local_sum(3)=local_sum(3)+s(k,j,i)*s(k,j,i)
        end do
      end do
    end do

    call mpi_allreduce(local_sum, global_sum, 3, PRECISION_TYPE, MPI_SUM, current_state%parallel%monc_communicator, ierr)
    inner_prod_three_way=global_sum
  end function inner_prod_three_way 

  !> Sets the values of the provided matrix to solve the poisson equation
  !! @param grid_configuration Configuration of the vertical and horizontal grids
  !! @param A The matrix that the values are written into
  !! @param z_size Number of elements in a column
  subroutine set_matrix_for_poisson(grid_configuration, A, z_size)
    type(grid_configuration_type), intent(inout) :: grid_configuration
    type(matrix_type), intent(inout) :: A
    integer, intent(in) :: z_size    

    integer :: k
    real(kind=DEFAULT_PRECISION) :: d_sc, concat_scalars    

    A%n=grid_configuration%horizontal%cx*grid_configuration%horizontal%cx
    A%s=A%n
    A%e=grid_configuration%horizontal%cy*grid_configuration%horizontal%cy
    A%w=A%e
    concat_scalars=A%n+A%s+A%e+A%w
    do k=2, z_size      
      if (symm_prob) then
         A%vol(k)=grid_configuration%vertical%dz(k)
         d_sc=1.0/grid_configuration%vertical%rhon(k)
      else
         d_sc=grid_configuration%vertical%rdz(k) / grid_configuration%vertical%rhon(k)
         A%vol(k)=1.0
      endif

      if (k==z_size) then
        A%u(k)=0.0_DEFAULT_PRECISION
      else
        A%u(k)=grid_configuration%vertical%rho(k)*grid_configuration%vertical%rdzn(k+1)
      end if
      if (k==2) then
        A%d(k)=0.0_DEFAULT_PRECISION
      else
        A%d(k)=grid_configuration%vertical%rho(k-1)*grid_configuration%vertical%rdzn(k)
      end if
      A%p(k) = d_sc * (-(A%u(k) + A%d(k))) - concat_scalars * A%vol(k)
      A%u(k)=d_sc * A%u(k)
      A%d(k)=d_sc * A%d(k)
    end do
     k=2 
     A%lu_d(k)=1.0_DEFAULT_PRECISION/A%p(k) 
     A%lu_u(k)=A%lu_d(k)*A%u(k) 
     do k=3, z_size 
       A%lu_d(k)=1.0_DEFAULT_PRECISION/(A%p(k) - A%d(k)*A%lu_u(k-1))
       A%lu_u(k)=A%u(k)*A%lu_d(k) 
     end do
  end subroutine set_matrix_for_poisson

  !> Determines the global divmax which is written into the current state
  !! @param current_state The current model state
  subroutine deduce_global_divmax(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    
    integer :: ierr
    
    call mpi_allreduce(current_state%local_divmax, current_state%global_divmax, 1, PRECISION_TYPE, MPI_MAX, &
         current_state%parallel%monc_communicator, ierr)
  end subroutine deduce_global_divmax   

  !> Copies the p field data to halo buffers for a specific process in a dimension and halo cell
  !! @param current_state The current model state
  !! @param neighbour_descriptions Description of the neighbour halo swapping status
  !! @param dim Dimension to copy from
  !! @param source_index The source index of the dimension we are reading from in the prognostic field
  !! @param pid_location Location of the neighbouring process in the local stored data structures
  !! @param current_page The current (next) buffer page to copy into
  !! @param source_data Optional source data which is read from
  subroutine copy_p_to_halo_buffer(current_state, neighbour_description, dim, source_index, &
       pid_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, pid_location, source_index
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, current_state%p%data, &
         dim, source_index, current_page(pid_location))

    current_page(pid_location)=current_page(pid_location)+1
  end subroutine copy_p_to_halo_buffer

  !> Copies the source field data to halo buffers for a specific process in a dimension and halo cell - for the calc_Ax halo swaps
  !! @param current_state The current model state
  !! @param neighbour_descriptions Description of the neighbour halo swapping status
  !! @param dim Dimension to copy from
  !! @param source_index The source index of the dimension we are reading from in the prognostic field
  !! @param pid_location Location of the neighbouring process in the local stored data structures
  !! @param current_page The current (next) buffer page to copy into
  !! @param source_data Optional source data which is read from
  subroutine copy_calc_Ax_to_halo_buffer(current_state, neighbour_description, dim, source_index, &
       pid_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, pid_location, source_index
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    type(field_data_wrapper_type) :: selected_source

    selected_source=source_data(1)

    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, selected_source%data, &
         dim, source_index, current_page(pid_location))

    current_page(pid_location)=current_page(pid_location)+1
  end subroutine copy_calc_Ax_to_halo_buffer

  !> Copies the halo buffer to halo location for the p field
  !! @param current_state The current model state
  !! @param neighbour_description The halo swapping description of the neighbour we are accessing the buffer of
  !! @param dim The dimension we receive for
  !! @param target_index The target index for the dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been read and copied already)
  !! @param source_data Optional source data which is written into
  subroutine copy_halo_buffer_to_p(current_state, neighbour_description, dim, target_index, &
       neighbour_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, current_state%p%data, &
         dim, target_index, current_page(neighbour_location))

    current_page(neighbour_location)=current_page(neighbour_location)+1
  end subroutine copy_halo_buffer_to_p

  !> Copies the halo buffer to halo location for the source field as required in the calc_Ax procedure
  !! @param current_state The current model state
  !! @param neighbour_description The halo swapping description of the neighbour we are accessing the buffer of
  !! @param dim The dimension we receive for
  !! @param target_index The target index for the dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been read and copied already)
  !! @param source_data Optional source data which is written into
  subroutine copy_halo_buffer_to_calc_Ax(current_state, neighbour_description, dim, target_index, &
       neighbour_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    type(field_data_wrapper_type) :: selected_source

    selected_source=source_data(1)

    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, selected_source%data, &
         dim, target_index, current_page(neighbour_location))

    current_page(neighbour_location)=current_page(neighbour_location)+1
  end subroutine copy_halo_buffer_to_calc_Ax

  !> Does local data copying for P variable halo swap
  !! @param current_state The current model state_mod
  !! @param source_data Optional source data which is written into
  subroutine perform_local_data_copy_for_p(current_state, halo_depth, involve_corners, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: halo_depth
    logical, intent(in) :: involve_corners
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call perform_local_data_copy_for_field(current_state%p%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
  end subroutine perform_local_data_copy_for_p

  !> Does a local data copy for halo swapping cells with wrap around (to maintain periodic boundary condition)
  !! @param current_state The current model state
  !! @param source_data Optional source data which is written into
  subroutine perform_local_data_copy_for_calc_Ax(current_state, halo_depth, involve_corners, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: halo_depth
    logical, intent(in) :: involve_corners
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    type(field_data_wrapper_type) :: selected_source

    selected_source=source_data(1)

    call perform_local_data_copy_for_field(selected_source%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
  end subroutine perform_local_data_copy_for_calc_Ax

  !> Creates a problem matrix, allocating the required data based upon the column size
  !! @param z_size Number of elements in the vertical column
  !! @returns The allocated matrix ready to be used
  function create_problem_matrix(z_size)
    integer, intent(in) :: z_size
    type(matrix_type) :: create_problem_matrix

    allocate(create_problem_matrix%u(z_size), create_problem_matrix%d(z_size), create_problem_matrix%p(z_size), &
         create_problem_matrix%lu_u(z_size), create_problem_matrix%lu_d(z_size), create_problem_matrix%vol(z_size))
  end function create_problem_matrix

  !> Completes the psrce calculation by waiting on all outstanding psrce communications to complete and then combine the
  !! received values with the P field for U and V
  !! @param current_state The current model state
  !! @param y_halo_size The halo size in the Y dimension
  !! @param x_halo_size The halo size in the X dimension
  subroutine complete_psrce_calculation(current_state, y_halo_size, x_halo_size)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: y_halo_size, x_halo_size

    integer :: ierr, combined_handles(2), i, j, k

    combined_handles(1)=current_state%psrce_x_hs_recv_request
    combined_handles(2)=current_state%psrce_y_hs_recv_request
    call mpi_waitall(2, combined_handles, MPI_STATUSES_IGNORE, ierr)

    do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
      do k=2,current_state%local_grid%size(Z_INDEX)
#ifdef U_ACTIVE
        current_state%p%data(k,j,x_halo_size+1)=current_state%p%data(k,j,x_halo_size+1)-&
               current_state%psrce_recv_buffer_x(k-1,j-x_halo_size)
#endif
#ifdef V_ACTIVE
        if (j .gt. y_halo_size+1) current_state%p%data(k, j, x_halo_size+1)=current_state%p%data(k, j, x_halo_size+1)-&
             current_state%global_grid%configuration%horizontal%cy * current_state%sv%data(k, j-1, x_halo_size+1)
#endif
      end do
    end do
      
#ifdef V_ACTIVE
    do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
      do k=2,current_state%local_grid%size(Z_INDEX)
        current_state%p%data(k,y_halo_size+1,i)=current_state%p%data(k,y_halo_size+1,i)-&
             current_state%psrce_recv_buffer_y(k-1,i-y_halo_size)
      end do
    end do
#endif

    combined_handles(1)=current_state%psrce_x_hs_send_request
    combined_handles(2)=current_state%psrce_y_hs_send_request
    call mpi_waitall(2, combined_handles, MPI_STATUSES_IGNORE, ierr)
  end subroutine complete_psrce_calculation 
end module iterativesolver_mod
