!> PETSc solver component to call out to PETSc for solving the Poisson equation for pressure
#include "petsc/finclude/petscvecdef.h"
#include "petsc/finclude/petscmatdef.h"
#include "petsc/finclude/petsckspdef.h"
#include "petsc/finclude/petscdmdadef.h"
module petsc_solver_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use state_mod, only : model_state_type
  use monc_component_mod, only : component_descriptor_type
  use grids_mod, only : local_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use logging_mod, only : LOG_DEBUG, LOG_INFO, LOG_ERROR, log_log, log_master_log, log_get_logging_level, log_is_master
  use optionsdatabase_mod, only : options_get_integer, options_get_real, options_get_string
  use conversions_mod, only : conv_to_string
  use communication_types_mod, only : halo_communication_type, neighbour_description_type, field_data_wrapper_type
  use halo_communication_mod, only : copy_buffer_to_field, copy_field_to_buffer, perform_local_data_copy_for_field, &
       init_halo_communication, finalise_halo_communication, blocking_halo_swap, get_single_field_per_halo_cell
  use petscvec
  use petscmat
  use petscksp
  use petscsys
  use petscdmda
  implicit none

#ifndef TEST_MODE
  private
#endif

  KSP :: ksp
  DM :: da
  integer :: z_start, z_end, y_start, y_end, x_start, x_end
  type(halo_communication_type), save :: halo_swap_state

  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: p_source, prev_p
  real(kind=DEFAULT_PRECISION) :: cx, cy
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: rdzn, rdz, rho, rhon, dz, dzn

  public petsc_solver_get_descriptor

  contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The PETSc solvre component descriptor
  type(component_descriptor_type) function petsc_solver_get_descriptor()
    petsc_solver_get_descriptor%name="petsc_solver"
    petsc_solver_get_descriptor%version=0.1
    petsc_solver_get_descriptor%initialisation=>init_callback
    petsc_solver_get_descriptor%timestep=>timestep_callback
    petsc_solver_get_descriptor%finalisation=>finalisation_callback
  end function petsc_solver_get_descriptor

 subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    PetscErrorCode :: ierr

    call PetscFinalize(ierr)
  end subroutine finalisation_callback  

  !> Called upon model initialisation. Will basically read from the options database and set options in
  !! the database that are appropriate
  !! @param current_state The current model state_mod
  subroutine init_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    PetscErrorCode :: ierr
    KSPType :: ksp_type
    PC :: pc_instance
    PCType :: pc_type
    integer :: y_lengths(current_state%parallel%dim_sizes(Y_INDEX)), x_lengths(current_state%parallel%dim_sizes(X_INDEX))
    integer new_communicator, my_transposed_rank, x_rank_val, y_rank_val

    call apply_dimension_bounds(current_state%global_grid%size(Y_INDEX), current_state%parallel%dim_sizes(Y_INDEX), y_lengths)
    call apply_dimension_bounds(current_state%global_grid%size(X_INDEX), current_state%parallel%dim_sizes(X_INDEX), x_lengths)

    allocate(p_source(current_state%local_grid%size(Z_INDEX)-1, current_state%local_grid%size(Y_INDEX), &
         current_state%local_grid%size(X_INDEX)), prev_p(current_state%local_grid%size(Z_INDEX)-1, &
         current_state%local_grid%size(Y_INDEX), current_state%local_grid%size(X_INDEX)))

    cx=current_state%global_grid%configuration%horizontal%cx
    cy=current_state%global_grid%configuration%horizontal%cy
    allocate(rdzn(current_state%local_grid%size(Z_INDEX)-1), rdz(current_state%local_grid%size(Z_INDEX)-1), &
         rho(current_state%local_grid%size(Z_INDEX)-1), rhon(current_state%local_grid%size(Z_INDEX)-1), &
         dz(current_state%local_grid%size(Z_INDEX)-1), dzn(current_state%local_grid%size(Z_INDEX)-1))
    rdzn=current_state%global_grid%configuration%vertical%rdzn(2:)
    rdz=current_state%global_grid%configuration%vertical%rdz(2:)
    rho=current_state%global_grid%configuration%vertical%rho(2:)
    rhon=current_state%global_grid%configuration%vertical%rhon(2:)
    dz=current_state%global_grid%configuration%vertical%dz(2:)
    dzn=current_state%global_grid%configuration%vertical%dzn(2:)

    z_start=current_state%local_grid%halo_size(Z_INDEX)+2
    z_end=current_state%local_grid%halo_size(Z_INDEX)+current_state%local_grid%size(Z_INDEX)
    y_start=current_state%local_grid%halo_size(Y_INDEX)+1
    y_end=current_state%local_grid%halo_size(Y_INDEX)+current_state%local_grid%size(Y_INDEX)
    x_start=current_state%local_grid%halo_size(X_INDEX)+1
    x_end=current_state%local_grid%halo_size(X_INDEX)+current_state%local_grid%size(X_INDEX)

    y_rank_val = current_state%parallel%my_rank / current_state%parallel%dim_sizes(X_INDEX)
    x_rank_val = mod(current_state%parallel%my_rank, current_state%parallel%dim_sizes(X_INDEX))    
    my_transposed_rank = x_rank_val*current_state%parallel%dim_sizes(Y_INDEX) + y_rank_val

    call MPI_Comm_split(current_state%parallel%monc_communicator, 1, my_transposed_rank, new_communicator, ierr)
    PETSC_COMM_WORLD = new_communicator

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)    
    call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
    call DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_STAR, &
         current_state%global_grid%size(Z_INDEX)-1, current_state%global_grid%size(Y_INDEX), &
         current_state%global_grid%size(X_INDEX), 1, current_state%parallel%dim_sizes(Y_INDEX), &
         current_state%parallel%dim_sizes(X_INDEX), 1, 1, (/ current_state%global_grid%size(Z_INDEX)-1 /), &
         y_lengths, x_lengths, da, ierr)
    call KSPSetDM(ksp, da, ierr)
    call KSPSetComputeRHS(ksp, compute_RHS, PETSC_NULL_OBJECT, ierr)
    call KSPSetComputeOperators(ksp, compute_matrix, PETSC_NULL_OBJECT, ierr)
    call KSPSetComputeInitialGuess(ksp, compute_initial_guess, PETSC_NULL_OBJECT, ierr)
    call PetscOptionsSetValue("-options_left", "no", ierr)
    if (trim(options_get_string(current_state%options_database, "solver_type")) .ne. "auto") then
      call PetscOptionsSetValue("-ksp_type", trim(options_get_string(current_state%options_database, "solver_type")), ierr)
    end if
    if (trim(options_get_string(current_state%options_database, "preconditioner_type")) .ne. "auto") then
      call PetscOptionsSetValue("-pc_type", trim(options_get_string(current_state%options_database, "preconditioner_type")), ierr)
    end if
    if (trim(options_get_string(current_state%options_database, "norm_type")) .ne. "auto") then
      if (trim(options_get_string(current_state%options_database, "norm_type")) .eq. "preconditioned" .or. &
           trim(options_get_string(current_state%options_database, "norm_type")) .eq. "unpreconditioned" .or. &
           trim(options_get_string(current_state%options_database, "norm_type")) .eq. "natural" .or. &
           trim(options_get_string(current_state%options_database, "norm_type")) .eq. "none") then
        call PetscOptionsSetValue("-ksp_norm_type", trim(options_get_string(current_state%options_database, "norm_type")), ierr)
      else
        call log_master_log(LOG_ERROR, "Configured PETSc norm type of '"//&
             trim(options_get_string(current_state%options_database, "norm_type"))//"' not recognised")
      end if      
    end if    
    call PetscOptionsSetValue("-ksp_rtol", conv_to_string(options_get_real(current_state%options_database, "tolerance")), ierr)
    if (options_get_integer(current_state%options_database, "max_iterations") .gt. 0) then
      call PetscOptionsSetValue("-ksp_max_it", &
           conv_to_string(options_get_integer(current_state%options_database, "max_iterations")), ierr)
    end if
    call KSPSetFromOptions(ksp, ierr)
    call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)

    if (log_is_master() .and. log_get_logging_level() == LOG_DEBUG) then
      call KSPMonitorSet(ksp,KSP_monitor,PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION,ierr)
    end if    
    prev_p=0.0_DEFAULT_PRECISION
    call init_halo_communication(current_state, get_single_field_per_halo_cell, halo_swap_state, 1, .false.)
    call KSPGetType(ksp, ksp_type, ierr)
    call KSPGetPC(ksp, pc_instance, ierr)
    call PCGetType(pc_instance, pc_type, ierr)
    call log_master_log(LOG_INFO, "PETSc iterative solver initialised, using "//trim(pc_type)//" preconditioner with "//&
         trim(ksp_type)//" solver")
  end subroutine init_callback

  !> Monitor procedure called on each iteration, this is provided only at debugging level
  !! @param ksp The PETSc KSP data structure which contains the context of the solver
  !! @param n The iteration
  !! @param rnorm The relative norm
  !! @param dummy Dummy argument to maintain consistency with the PETSc C interface
  !! @param ierr PETSc error code
  subroutine KSP_monitor(ksp, n, rnorm, dummy, ierr)
    KSP              ksp
    PetscErrorCode ierr
    PetscInt n,dummy
    PetscReal rnorm

    call log_log(LOG_DEBUG, "Iteration="//trim(conv_to_string(n))//" Rnorm="//trim(conv_to_string(rnorm, &
         exponent_small_numbers=.true.)))
  end subroutine KSP_monitor  

  !> Timestep hook which is called at each timestep to solve the pressure equation via PETSc. This will call to set and 
  !! precondition the matrix on the first call (first timestep) only, and will set the RHS and X initial guess on each
  !! call as one would expect. At the end a halo swap is performed on p, which is needed for pstep 
  !! @param current_state The current model state_mod
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    
    PetscErrorCode :: ierr
    Vec x
    PetscScalar, pointer :: xx(:)
    integer :: its, i, j, k
    PetscReal :: rnorm
    KSPConvergedReason :: reason

    call complete_psrce_calculation(current_state, current_state%local_grid%halo_size(Y_INDEX), &
         current_state%local_grid%halo_size(X_INDEX))    

    p_source=current_state%p%data(z_start:z_end, y_start:y_end, x_start:x_end)
    call KSPSolve(ksp, PETSC_NULL_OBJECT, PETSC_NULL_OBJECT, ierr)
    call KSPGetSolution(ksp, x, ierr)
    call VecGetArrayReadF90(x, xx, ierr)
    
    call copy_petsc_pointer_to_data(xx, z_start, z_end, y_start, y_end, x_start, x_end, current_state%p%data)
    call VecRestoreArrayReadF90(x, xx, ierr)
    
    prev_p=current_state%p%data(z_start:z_end, y_start:y_end, x_start:x_end)

    if (log_is_master() .and. log_get_logging_level() == LOG_DEBUG) then
      call KSPGetIterationNumber(ksp, its, ierr)
      call KSPGetResidualNorm(ksp, rnorm, ierr)
      call KSPGetConvergedReason(ksp, reason, ierr)
      call log_log(LOG_DEBUG, "Converged in "//trim(conv_to_string(its))//" rnorm="//trim(conv_to_string(&
           rnorm, exponent_small_numbers=.true.))// " ("//trim(determine_convegence_reason(reason))//")")
    end if
    call blocking_halo_swap(current_state, halo_swap_state, copy_p_to_halo_buffer, &
         perform_local_data_copy_for_p, copy_halo_buffer_to_p)
  end subroutine timestep_callback

  character(len=STRING_LENGTH) function determine_convegence_reason(r_code)
    integer, intent(in) :: r_code

    if (r_code == 1 .or. r_code == 2) then
      determine_convegence_reason="relative tolerance of residual norm"
    else if (r_code == 9 .or. r_code == 3) then
      determine_convegence_reason="residual norm less than absolute tolerance"
    else if (r_code == 4) then
      determine_convegence_reason="number of iterations exceeded maximum"
    else if (r_code == -3) then
      determine_convegence_reason="diverged as number of iterations exceeded maximum before any convergence criteria"
    else if (r_code == -4) then
      determine_convegence_reason="diverged as divergence tolerance exceeded"
    else if (r_code .gt. 0) then
      determine_convegence_reason="convergence code of "//conv_to_string(r_code)
    else if (r_code .lt. 0) then
      determine_convegence_reason="divergence code of "//conv_to_string(r_code)
    end if    
  end function determine_convegence_reason  

  !> Determines the initial guess, this is called from within PETSc and sets it to be the last value of P or zero if
  !! it is the first timestep
  !! @param ksp The PETSc KSP data structure which contains the context of the solver
  !! @param X The X vector to create an initial guess for
  !! @param dummy Dummy argument to maintain consistency with the PETSc C interface
  !! @param ierr PETSc error code
  subroutine compute_initial_guess(ksp, x, dummy, ierr)
    PetscErrorCode :: ierr
    KSP :: ksp
    Vec :: x
    integer :: dummy(*)

    DM dm
    PetscScalar, pointer :: xx(:)
    PetscInt :: zs, ys, xs, zm, ym, xm

    call KSPGetDM(ksp, dm, ierr)
    call DMDAGetCorners(dm, zs, ys, xs, zm, ym, xm, ierr)
    call VecGetArrayF90(x, xx, ierr)
    call copy_data_to_petsc_pointer(xx, zs, ys, xs, zm, ym, xm, prev_p)
    call VecRestoreArrayF90(x, xx, ierr)
  end subroutine compute_initial_guess  

  !> Callback issued from the PETSc library to compute the RHS, this is called every timestep (as we have a different RHS)
  !! @param ksp The KSP data structure which contains the solver context
  !! @param b The RHS
  !! @param dummy Dummy argument for compatability with the C PETSc interface
  !! @param ierr Error code
  subroutine compute_RHS(ksp, b, dummy, ierr)
    PetscErrorCode :: ierr
    KSP :: ksp
    Vec :: b
    integer :: dummy(*)
    
    DM dm
    PetscScalar, pointer :: xx(:)
    PetscInt :: zs, ys, xs, zm, ym, xm

    call KSPGetDM(ksp, dm, ierr)
    call DMDAGetCorners(dm, zs, ys, xs, zm, ym, xm, ierr)
    call VecGetArrayF90(b, xx, ierr)
    call copy_data_to_petsc_pointer(xx, zs, ys, xs, zm, ym, xm, p_source)
    call VecRestoreArrayF90(b, xx, ierr)
  end subroutine compute_RHS

  !> Copies data into a data pointer which is provided by PETSc, doing it this ways allows us to give a shape to the pointer
  !! data and hence the data copy is simpler code wise
  !! @param x The target pointer to copy into
  !! @param zs Local starting point in Z
  !! @param ys Local starting point in Y
  !! @param xs Local starting point in X
  !! @param zm Number of points held locally in Z
  !! @param ym Number of points held locally in Y
  !! @param xm Number of points held locally in X
  !! @param src_values The data values to copy from into this target data array
  subroutine copy_data_to_petsc_pointer(x, zs, ys, xs, zm, ym, xm, src_values)
    integer, intent(in) :: zs, ys, xs, zm, ym, xm
    PetscScalar, intent(inout) :: x(zs:zs+zm-1, ys:ys+ym-1, xs:xs+xm-1)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: src_values

    integer :: i, j, k

    do i=xs, xs+xm-1
      do j=ys, ys+ym-1
        do k=zs, zs+zm-1
          x(k, j, i)=src_values((k-zs)+1, (j-ys)+1, (i-xs)+1)
        end do
      end do
    end do
  end subroutine copy_data_to_petsc_pointer

  !> Copies the data in a pointer that was provided by PETSc into some target data, doing it this way means we can give a shape
  !! to the pointer
  !! @param x The pointer provided by PETSc that we copy from
  !! @param z_start Starting point for the target array in Z
  !! @param z_end Ending point for the target array in Z
  !! @param y_start Starting point for the target array in Y
  !! @param y_end Ending point for the target array in Y
  !! @param x_start Starting point for the target array in X
  !! @param x_end Ending point for the target array in X
  !! @param target_values Target array to copy into
  subroutine copy_petsc_pointer_to_data(x, z_start_idx, z_end_idx, y_start_idx, y_end_idx, x_start_idx, x_end_idx, target_values)
     integer, intent(in) :: z_start_idx, z_end_idx, y_start_idx, y_end_idx, x_start_idx, x_end_idx
     PetscScalar :: x((z_end_idx-z_start_idx)+1, (y_end_idx-y_start_idx)+1, (x_end_idx-x_start_idx)+1)
     real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: target_values

     target_values(z_start_idx:z_end_idx, y_start_idx:y_end_idx, x_start_idx:x_end_idx)=x
  end subroutine copy_petsc_pointer_to_data

  !> Call back issued from within PETSc to create the matrix, this is only called once (the first PETSc run)
  !! @param ksp The solver data structure which provides context for this run
  !! @param A Matrix A this is the linear operator
  !! @param B Matrix B this is the preconditioning matrix
  !! @param dummy Dummy argument to ensure consistency with C PETSc interface
  !! @param ierr The PETSc error code
  subroutine compute_matrix(ksp, A, B, dummy, ierr)
    PetscErrorCode :: ierr
    KSP :: ksp
    Mat :: A, B
    integer :: dummy(*)

    DM dm
    PetscInt :: zs, ys,xs, zm, ym, xm, mz, my, mx
    real(kind=DEFAULT_PRECISION) ::  v(7), Hx, Hy
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::  HzU, HzD, d_sc  
    MatStencil :: row(4),col(4,7)
    integer :: i, j, k, num

    call KSPGetDM(ksp, dm, ierr)
    call DMDAGetInfo(dm,PETSC_NULL_INTEGER, mz, my, mx, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
     PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
    call DMDAGetCorners(dm, zs, ys, xs, zm, ym, xm, ierr)

    allocate(HzU(0:mz-1), HzD(0:mz-1), d_sc(0:mz-1))

    do k=0, mz-1
      d_sc(k)=rdz(k+1) / rhon(k+1)
      if (k .gt. 0) then
        HzD(k)=rho(k) * rdzn(k+1)
      else
        HzD(k)=0.0_DEFAULT_PRECISION
      end if
      if (k .lt. mz-1) then
        HzU(k)=rho(k+1) * rdzn(k+2)
      else
        HzU(k)=0.0_DEFAULT_PRECISION
      end if
    end do   

    Hx = cx * cx
    Hy = cy * cy  
    do i=xs, xs+xm-1
      do j=ys,ys+ym-1
        do k=zs, zs+zm-1
          row(MatStencil_i) = k
          row(MatStencil_j) = j
          row(MatStencil_k) = i

          v(1)=Hy
          col(MatStencil_i, 1)=k
          col(MatStencil_j, 1)=j-1
          col(MatStencil_k, 1)=i
          v(2)=Hx
          col(MatStencil_i, 2)=k
          col(MatStencil_j, 2)=j
          col(MatStencil_k, 2)=i-1
          v(3)=d_sc(k) * (-(HzU(k) + HzD(k))) - (2.0*(Hx + Hy))
          col(MatStencil_i, 3)=k
          col(MatStencil_j, 3)=j
          col(MatStencil_k, 3)=i
          v(4)=Hx
          col(MatStencil_i, 4)=k
          col(MatStencil_j, 4)=j 
          col(MatStencil_k, 4)=i+1
          v(5)=Hy
          col(MatStencil_i, 5)=k
          col(MatStencil_j, 5)=j+1
          col(MatStencil_k, 5)=i
          num=6
          if (k .gt. zs) then
            v(num)=(d_sc(k)*HzD(k))
            col(MatStencil_i, num)=k-1
            col(MatStencil_j, num)=j 
            col(MatStencil_k, num)=i
            num=num+1
          end if
          if (k .lt. mz-1) then
            v(num)=(d_sc(k)*HzU(k))
            col(MatStencil_i, num)=k+1
            col(MatStencil_j, num)=j
            col(MatStencil_k, num)=i
            num=num+1
          end if
          call MatSetValuesStencil(B, 1, row, num-1, col, v, INSERT_VALUES, ierr)
        end do
      end do
    end do
    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY, ierr)
    if ( A .ne. B) then
      call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
    endif
    deallocate(HzU, HzD, d_sc)
  end subroutine compute_matrix

  !> Applies dimension bounds to determine the number of elements held locally for each process in that dimension
  !! @param dim_size The global size of the dimension
  !! @param dim_processes The number of processes in that dimension
  !! @param length_array Output array of size dim_processes which contains the number of elements held locally for each process
  subroutine apply_dimension_bounds(dim_size, dim_processes, length_array)
    integer, intent(in) :: dim_size, dim_processes
    integer, intent(out) :: length_array(:)

    integer :: i, dimension_division, dimension_extra

    dimension_division = dim_size / dim_processes
    length_array= dimension_division    
    
    dimension_extra = dim_size - (dimension_division * dim_processes)
    if (dimension_extra .gt. 0) then
      do i=1, dimension_extra
        length_array(i)=length_array(i)+1
      end do
    end if
  end subroutine apply_dimension_bounds

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
end module petsc_solver_mod
