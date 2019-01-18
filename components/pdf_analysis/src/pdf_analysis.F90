!> Calculates fields related to distributions of data on full-domain horizontal 2d slices 
module pdf_analysis_mod
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type
  use state_mod, only : model_state_type
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use optionsdatabase_mod, only : options_has_key, options_get_logical, options_get_integer, options_get_string, options_get_real
  use mpi, only : MPI_SUM, MPI_IN_PLACE, MPI_INT, MPI_REAL, MPI_DOUBLE, MPI_Comm
  use logging_mod, only : LOG_INFO, LOG_DEBUG, LOG_ERROR, log_master_log, log_is_master
  use conversions_mod, only : conv_to_string
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: start_x, end_x, start_y, end_y, xsize, ysize

  real(kind=DEFAULT_PRECISION) :: uppercrit, dwnpercrit
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: tmp_all
  integer, dimension(:), allocatable :: gpts_on_proc, & ! number of horizontal grid points on each process
                                        displacements   ! displacement for mpi_gatherv
                                                        ! these are available on all processes

  integer :: tpts  ! total number of horizontal grid points on full domain
  integer :: lpts  ! local number of horizontal grid points on 

  logical :: show_critical_w  ! stdout diagnostic logical

  integer :: diagnostic_generation_frequency

  public pdf_analysis_get_descriptor

contains


  !> Returns the component descriptor of pdf analysis  module
  !! @returns Component descriptor of pdf analysis
  type(component_descriptor_type) function pdf_analysis_get_descriptor()
    pdf_analysis_get_descriptor%name="pdf_analysis"
    pdf_analysis_get_descriptor%version=0.1
    pdf_analysis_get_descriptor%initialisation=>init_callback
    pdf_analysis_get_descriptor%timestep=>timestep_callback
    pdf_analysis_get_descriptor%finalisation=>finalisation_callback

    pdf_analysis_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    pdf_analysis_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(pdf_analysis_get_descriptor%published_fields(2))
    pdf_analysis_get_descriptor%published_fields(1)="critical_updraft_local"
    pdf_analysis_get_descriptor%published_fields(2)="critical_downdraft_local"

  end function pdf_analysis_get_descriptor


  !> Called on MONC initialisation, will allocate appropriate data structures
  !! @param current_state The current model state
  subroutine init_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: ierr, inc

    tpts = current_state%global_grid%size(X_INDEX)*current_state%global_grid%size(Y_INDEX)

    start_x = current_state%local_grid%local_domain_start_index(X_INDEX)
    end_x   = current_state%local_grid%local_domain_end_index(X_INDEX)
    start_y = current_state%local_grid%local_domain_start_index(Y_INDEX)
    end_y   = current_state%local_grid%local_domain_end_index(Y_INDEX)

    xsize = end_x - start_x + 1
    ysize = end_y - start_y + 1

    lpts = xsize*ysize

    uppercrit  = options_get_real(current_state%options_database, "uppercrit")
    dwnpercrit = options_get_real(current_state%options_database, "dwnpercrit")

    show_critical_w = options_get_logical(current_state%options_database, "show_critical_w")

    !> Allocate space for the global 2d field only on a single process
!    if (current_state%parallel%my_rank == 0) 
     allocate(tmp_all(tpts))
!    else
!      allocate(tmp_all(1))
!    end if

    !> Allocate and collect horizontal local sizes, send to all proceses
    allocate(gpts_on_proc(current_state%parallel%processes))
    call mpi_allgather(lpts, 1, MPI_INT, gpts_on_proc, 1, MPI_INT, current_state%parallel%monc_communicator, ierr)

    !> Allocate and initialize displacement values
    allocate(displacements(current_state%parallel%processes)) 
    displacements(1) = 0
    do inc = 2, current_state%parallel%processes
      displacements(inc) = displacements(inc-1) + gpts_on_proc(inc-1)
    end do ! loop over processes

    !> Allocate critial fields in current_state if a cold start
    if (.not. current_state%continuation_run) then
      allocate(current_state%global_grid%configuration%vertical%w_dwn(current_state%local_grid%size(Z_INDEX)),&
               current_state%global_grid%configuration%vertical%w_up(current_state%local_grid%size(Z_INDEX)))
    end if

    ! Save the sampling_frequency to force diagnostic calculation on select time steps
    diagnostic_generation_frequency=options_get_integer(current_state%options_database, "sampling_frequency")

 end subroutine init_callback


  !> Will sort the values across the whole domain and calculate the value corresponding to the percentile threshold
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    !> Current forumulation only handles vertical velocity percentiles.
    !! Future enhancements may employ this component to perform additional 
    !! operations that require access to full horizontal fields, such as
    !! pdf calculations.

    if (mod(current_state%timestep, diagnostic_generation_frequency) == 0) call calculate_w_percentiles(current_state)

  end subroutine timestep_callback


  !> Frees up the temporary data for the tm_allp
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(tmp_all)) deallocate(tmp_all)    
  end subroutine finalisation_callback


  !> Calculates the w percentiles over the whole domain and stores these in the w up/dwn percentile arrays of current_state
  !! @param current_state The current model state
  subroutine calculate_w_percentiles(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(lpts) :: tmp_var

    integer :: i, j, k, num_neg, num_pos, dd_thresh_pos, ud_thresh_pos
    integer :: max_up_k, min_dwn_k 
    real(kind=DEFAULT_PRECISION), dimension((lpts+1)/2) :: T
    real(kind=DEFAULT_PRECISION), dimension((tpts+1)/2) :: Tall
    real(kind=DEFAULT_PRECISION)                        :: max_up, min_dwn, &
                                                           max_up_th, min_dwn_th
    integer :: ierr
    real(kind=DEFAULT_PRECISION), dimension(ysize,xsize) :: l2d ! local 2d data

    !> initialize diagnostic thresholds
    max_up_th  = 0.0_DEFAULT_PRECISION
    min_dwn_th = 0.0_DEFAULT_PRECISION
    max_up     = max_up_th
    min_dwn    = min_dwn_th
    max_up_k   = 0
    min_dwn_k  = 0

    !> reset thresholds
    current_state%global_grid%configuration%vertical%w_dwn(:) = 0.0_DEFAULT_PRECISION 
    current_state%global_grid%configuration%vertical%w_up(:)  = 0.0_DEFAULT_PRECISION

    !> Loop over levels
    do k = 2, current_state%local_grid%size(Z_INDEX)

       !> specify local data area
       l2d(:,:)=0.5_DEFAULT_PRECISION*(   current_state%w%data(k  , start_y:end_y, start_x:end_x)   &
                                        + current_state%w%data(k-1, start_y:end_y, start_x:end_x)    )

       !> Reshape to 1-D array
       tmp_var=pack(l2d,.true.)

       !> Perform sort of data on local process
       call MergeSort(tmp_var,lpts,T)

       !> Gather 2d field to single process
       call mpi_gatherv(tmp_var, lpts, PRECISION_TYPE, tmp_all, gpts_on_proc, displacements, PRECISION_TYPE, &
                        0, current_state%parallel%monc_communicator, ierr )

       !> Perform global operations
       if (current_state%parallel%my_rank == 0) then

         !> Sort the global data on single process
         call MergeSort(tmp_all,tpts,Tall)

         !> Determine threshold updraft and downdraft values
         num_neg = count(tmp_all < 0.0_DEFAULT_PRECISION)
         num_pos = count(tmp_all > 0.0_DEFAULT_PRECISION)

         dd_thresh_pos = int(num_neg * dwnpercrit) 
         ud_thresh_pos = tpts - int(num_pos * uppercrit) + 1

         if ( dd_thresh_pos == 0 ) dd_thresh_pos = 1
         if ( ud_thresh_pos == 0 .or. num_pos == 0 ) ud_thresh_pos = tpts

         current_state%global_grid%configuration%vertical%w_dwn(k) = tmp_all(dd_thresh_pos)
         current_state%global_grid%configuration%vertical%w_up(k)  = tmp_all(ud_thresh_pos)

         !> do some stdout diagnostic work
         if (show_critical_w) then
           if ( tmp_all(dd_thresh_pos) < min_dwn_th ) then
             min_dwn_th = tmp_all(dd_thresh_pos)
             min_dwn = tmp_all(1)  ! sorted array
             min_dwn_k = k
           end if 
           if ( tmp_all(ud_thresh_pos) > max_up_th ) then
             max_up_th = tmp_all(ud_thresh_pos)
             max_up = tmp_all(tpts)  ! sorted array
             max_up_k = k
           end if
         end if ! show_critical_w

       end if ! global operations section

    end do ! loop over k

    !> Inform all processes of calculated thresholds
    call mpi_bcast(current_state%global_grid%configuration%vertical%w_dwn(:), current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%global_grid%configuration%vertical%w_up(:),  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)


    !> Display some diagnostics, if requested
    if (show_critical_w) then
      call log_master_log(LOG_INFO, 'Time:  '//trim(conv_to_string(current_state%time))//' s')
      call log_master_log(LOG_INFO, 'Maximum updraft threshold:   '&
                          //trim(conv_to_string(max_up_th))//' found at level '//trim(conv_to_string(max_up_k)) )
      call log_master_log(LOG_INFO, 'Maximum updraft:   '&
                          //trim(conv_to_string(max_up))//' at level '//trim(conv_to_string(max_up_k)) )
      call log_master_log(LOG_INFO, 'Minimum downdraft threshold:   '&
                          //trim(conv_to_string(min_dwn_th))//' found at level '//trim(conv_to_string(min_dwn_k)) )   
      call log_master_log(LOG_INFO, 'Minimum downdraft:   '&
                          //trim(conv_to_string(min_dwn))//' at level '//trim(conv_to_string(min_dwn_k)) )
    end if ! show_critical_w

  end subroutine calculate_w_percentiles   


  !> Combines with MergeSort sorting algorithm taken from:
  !  https://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
  !  and modified to match local type and renamed to avoid confusion with intrinsic merge
  !  All parameters based on MergeSort.  No need to modify anything.
  subroutine MergeSortMerge(A,NA,B,NB,C,NC)
 
    integer, intent(in) :: NA,NB,NC                              ! Normal usage: NA+NB = NC
    real(kind=DEFAULT_PRECISION), intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
    real(kind=DEFAULT_PRECISION), intent(in)     :: B(NB)
    real(kind=DEFAULT_PRECISION), intent(in out) :: C(NC)
 
    integer :: I,J,K
 
    I = 1; J = 1; K = 1;
    do while(I <= NA .and. J <= NB)
      if (A(I) <= B(J)) then
         C(K) = A(I)
         I = I+1
      else
         C(K) = B(J)
         J = J+1
      endif
      K = K + 1
    enddo
    do while (I <= NA)
      C(K) = A(I)
      I = I + 1
      K = K + 1
    enddo
    return
 
  end subroutine mergesortmerge
 
  !> Combines with MergeSortMerge sorting algorithm taken from:
  !  https://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
  !  and modified to match local type
  !! @A array of values to be sorted, returned sorted
  !! @N size of A
  !! @T I don't really understand T
  recursive subroutine MergeSort(A,N,T)
 
    integer, intent(in) :: N
    real(kind=DEFAULT_PRECISION), dimension(N), intent(in out) :: A
    real(kind=DEFAULT_PRECISION), dimension((N+1)/2), intent (out) :: T
 
    integer :: NA,NB
    real(kind=DEFAULT_PRECISION) :: V
 
    if (N < 2) return
    if (N == 2) then
      if (A(1) > A(2)) then
        V = A(1)
        A(1) = A(2)
        A(2) = V
      endif
      return
    endif      
    NA=(N+1)/2
    NB=N-NA
 
    call MergeSort(A,NA,T)
    call MergeSort(A(NA+1),NB,T)
 
    if (A(NA) > A(NA+1)) then
      T(1:NA)=A(1:NA)
      call MergeSortMerge(T,NA,A(NA+1),NB,A,N)
    endif
    return
 
  end subroutine MergeSort


  !> Field information retrieval callback, this returns information for a specific components published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve information for
  !! @param field_information Populated with information about the field
  subroutine field_information_retrieval_callback(current_state, name, field_information)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_information

    field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
    field_information%number_dimensions=1
    field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
    field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
    field_information%enabled=.true.

  end subroutine field_information_retrieval_callback



  !> Field value retrieval callback, this returns the value of a specific published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve the value for
  !! @param field_value Populated with the value of the field
  subroutine field_value_retrieval_callback(current_state, name, field_value)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_value_type), intent(out) :: field_value

    if      (name .eq. "critical_updraft_local") then
      allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)), &
               source=current_state%global_grid%configuration%vertical%w_up(:))
    else if (name .eq. "critical_downdraft_local") then
      allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)), &
               source=current_state%global_grid%configuration%vertical%w_dwn(:))
    end if 
  end subroutine field_value_retrieval_callback

end module pdf_analysis_mod
