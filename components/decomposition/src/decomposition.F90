!> Parallel decomposition to determine the grid points and data columns that
!! are located on this process.
module decomposition_mod
  use datadefn_mod, only : STRING_LENGTH
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type
  use grids_mod, only : local_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use collections_mod, only : map_type
  use conversions_mod, only : conv_to_string
  use logging_mod, only : LOG_INFO, LOG_DEBUG, LOG_WARN, log_get_logging_level,&
       log_master_log, log_log
  use optionsdatabase_mod, only : options_get_string, options_get_integer
  use mpi
  implicit none

#ifndef TEST_MODE
  private
#endif

  public decomposition_get_descriptor

contains

  !> Provides the component descriptor for the core to register
  !! @returns The descriptor describing this component
  type(component_descriptor_type) function decomposition_get_descriptor()
    decomposition_get_descriptor%name="decomposition"
    decomposition_get_descriptor%version=0.1
    decomposition_get_descriptor%initialisation=>init_callback
  end function decomposition_get_descriptor

  !> The initialisation hook. Will set up the appropriate decomposition
  !! @param current_state The current model state_mod
  subroutine init_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    character(len=STRING_LENGTH) :: method

    method=options_get_string(current_state%options_database, "decomposition_method")

    if (method .eq. "onedim") then
       current_state%parallel%decomposition_procedure => one_dim_decomposition
    else if (method .eq. "twodim") then
       current_state%parallel%decomposition_procedure => two_dim_decomposition
    else
       current_state%parallel%decomposition_procedure => serial_decomposition
       if (method .ne. "serial") call log_log(LOG_WARN, "Decomposition method "//trim(method)//&
            " not recognised so defaulting to serial")
    end if
  end subroutine init_callback

  !> Decomposition into two dimensions
  !!
  !! Will currently balance between the x and y dimension as evenly as possible with any extra (odd numbered)
  !! processes being added to the X dimension
  !! @param current_state The current model state
  subroutine two_dim_decomposition(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: ierr
    integer, dimension(2) :: coords, distributed_dims

    current_state%local_grid%active = current_state%global_grid%active
    current_state%local_grid%dimensions = current_state%global_grid%dimensions

    if (current_state%global_grid%dimensions .ne. 3) then
       call log_master_log(LOG_WARN, "Two dimension decomposition selected with only "// &
            trim(conv_to_string(current_state%global_grid%dimensions)) //" so defaulting to one dimension instead")
       call one_dim_decomposition(current_state)
       return
    end if

    distributed_dims= (/0,0/)
    call mpi_dims_create(current_state%parallel%processes, 2, distributed_dims, ierr)

    current_state%parallel%dim_sizes(Z_INDEX)=1
    current_state%parallel%dim_sizes(Y_INDEX) = distributed_dims(1)
    current_state%parallel%dim_sizes(X_INDEX) = distributed_dims(2)

    if (.not. does_two_dim_work_with_domain_size(current_state, current_state%parallel%dim_sizes(Y_INDEX), &
         current_state%parallel%dim_sizes(X_INDEX))) then
       call log_master_log(LOG_WARN, "Defaulting to one dimension decomposition due to solution size too small")
       call one_dim_decomposition(current_state)
       return
    end if
    
    call mpi_cart_create(current_state%parallel%monc_communicator, 2, (/distributed_dims(1), distributed_dims(2)/),&
         (/.false., .false./), .false., current_state%parallel%neighbour_comm, ierr)
    call mpi_cart_coords(current_state%parallel%neighbour_comm, current_state%parallel%my_rank, 2, coords, ierr)

    call apply_halo_information_and_allocate_neighbours(current_state)
    call apply_z_dimension_information(current_state)

    ! sets the limits for each rank in X and Y
    call apply_dimension_bounds(current_state, Y_INDEX, merge(current_state%global_grid%size(Y_INDEX), 1, &
         current_state%global_grid%active(Y_INDEX)), distributed_dims(1), coords(1))
    call apply_dimension_bounds(current_state, X_INDEX, merge(current_state%global_grid%size(X_INDEX), 1, &
         current_state%global_grid%active(X_INDEX)), distributed_dims(2), coords(2))

    ! get all information about the neighbours and corner neighbours
    call apply_two_dim_neighbour_information(current_state, distributed_dims(1), distributed_dims(2), coords)

    call apply_data_start_end_bounds(current_state%local_grid)

    if (log_get_logging_level() .ge. LOG_DEBUG) then
       call log_log(LOG_DEBUG, "PID "//trim(conv_to_string(current_state%parallel%my_rank))//": y="//&
            trim(conv_to_string(current_state%parallel%my_coords(Y_INDEX)))//&
            " x="//trim(conv_to_string(current_state%parallel%my_coords(X_INDEX)))//" ny-1="//&
            trim(conv_to_string(current_state%local_grid%neighbours(Y_INDEX,1)))//&
            " ny+1="//trim(conv_to_string(current_state%local_grid%neighbours(Y_INDEX,3)))//" nx-1="//&
            trim(conv_to_string(current_state%local_grid%neighbours(X_INDEX,1)))//" nx+1="//&
            trim(conv_to_string(current_state%local_grid%neighbours(X_INDEX,3))))
    end if
    call display_decomposition_information(current_state, "TwoDim")
  end subroutine two_dim_decomposition

  !> Will apply the two dimensional neighbour information. This implements the wrap around
  !! aspect if there is not an immediate neighbour in a dimensional direction.
  !! Currently assumes that the halo is on the entire neighbour process
  !! @param current_state The current model state_mod
  !! @param yProcs Number of processes in the y dimension
  !! @param xProcs Number of processes in the x dimension
  !! @param coords Y-X coordinates of the current process
  subroutine apply_two_dim_neighbour_information(current_state, y_procs, x_procs, coords)

    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: y_procs, x_procs, coords(2)

    integer :: ierr, i, halo_depth

    current_state%parallel%wrapped_around=.false.

    call mpi_cart_shift(current_state%parallel%neighbour_comm, 0, 1, current_state%local_grid%neighbours(Y_INDEX,1), &
         current_state%local_grid%neighbours(Y_INDEX,3), ierr)
    if (current_state%local_grid%neighbours(Y_INDEX,1) .lt. 0) then
      call mpi_cart_rank(current_state%parallel%neighbour_comm,&
           (/y_procs-1, coords(2)/),  current_state%local_grid%neighbours(Y_INDEX,1), ierr)
      current_state%parallel%wrapped_around(Y_INDEX, 1)=.true.
    end if
    
    if (current_state%local_grid%neighbours(Y_INDEX,3) .lt. 0) then
      call mpi_cart_rank(current_state%parallel%neighbour_comm, &
         (/0, coords(2)/),  current_state%local_grid%neighbours(Y_INDEX,3), ierr)
      current_state%parallel%wrapped_around(Y_INDEX, 2)=.true.
    end if    

    ! find out who is on the left and on the right
    call mpi_cart_shift(current_state%parallel%neighbour_comm, 1, 1, &
         current_state%local_grid%neighbours(X_INDEX,1), &
         current_state%local_grid%neighbours(X_INDEX,3), ierr)
    if (current_state%local_grid%neighbours(X_INDEX,1) .lt. 0) then
      call mpi_cart_rank(current_state%parallel%neighbour_comm, &
         (/coords(1), x_procs-1/),  current_state%local_grid%neighbours(X_INDEX,1), ierr)
      current_state%parallel%wrapped_around(X_INDEX, 1)=.true.
    end if
    if (current_state%local_grid%neighbours(X_INDEX,3) .lt. 0) then
      call mpi_cart_rank(current_state%parallel%neighbour_comm, &
         (/coords(1), 0/),  current_state%local_grid%neighbours(X_INDEX,3), ierr)
      current_state%parallel%wrapped_around(X_INDEX, 2)=.true.
    end if

    ! the second column is on the same rank than the first one
    halo_depth = options_get_integer(current_state%options_database, "halo_depth")
    
    do i = 1, halo_depth -1
       current_state%local_grid%neighbours(Y_INDEX,i+1) = &
            current_state%local_grid%neighbours(Y_INDEX,i)
       current_state%local_grid%neighbours(Y_INDEX,i+halo_depth + 1) = &
            current_state%local_grid%neighbours(Y_INDEX,i+halo_depth)
       current_state%local_grid%neighbours(X_INDEX,i+1) = &
            current_state%local_grid%neighbours(X_INDEX,i)
       current_state%local_grid%neighbours(X_INDEX,i+halo_depth + 1) = &
            current_state%local_grid%neighbours(X_INDEX,i+halo_depth)
    end do
    ! obtain the rank of the corner neighbours
    call mpi_cart_rank(current_state%parallel%neighbour_comm, &
         (/merge(coords(1)-1, y_procs-1, coords(1) .ge. 1),   &
         merge(coords(2)-1, x_procs-1, coords(2) .ge. 1)/),   &
         current_state%local_grid%corner_neighbours(1,1), ierr)
    call mpi_cart_rank(current_state%parallel%neighbour_comm, &
         (/merge(coords(1)-1, y_procs-1, coords(1) .ge. 1),   &
         merge(coords(2)+1, 0, coords(2) .lt. x_procs-1)/),   &
         current_state%local_grid%corner_neighbours(2,1), ierr)
    call mpi_cart_rank(current_state%parallel%neighbour_comm, &
         (/merge(coords(1)+1, 0, coords(1) .lt. y_procs-1),   &
         merge(coords(2)-1, x_procs-1, coords(2) .ge. 1)/),   &
         current_state%local_grid%corner_neighbours(3,1), ierr)
    call mpi_cart_rank(current_state%parallel%neighbour_comm, &
         (/merge(coords(1)+1, 0, coords(1) .lt. y_procs-1),   &
         merge(coords(2)+1, 0, coords(2) .lt. x_procs-1)/),   &
         current_state%local_grid%corner_neighbours(4,1), ierr)
    
    !! TODO: hardcoded for halo depth two?
    current_state%local_grid%corner_neighbours(:,2)=current_state%local_grid%corner_neighbours(:,1)
  end subroutine apply_two_dim_neighbour_information

  !> Determines whether or not the planned two dimensional decomposition will fit into the
  !! X and Y global solutions size
  !! @param current_state The current model state
  !! @param y_dims Parallel process dimensions in Y
  !! @param x_dims Parallel process dimensions in X
  logical function does_two_dim_work_with_domain_size(current_state, y_dims, x_dims)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: y_dims, x_dims

    if (current_state%global_grid%active(Y_INDEX)) then
       if (floor(real(current_state%global_grid%size(Y_INDEX)) / y_dims) .lt. 2) then
          does_two_dim_work_with_domain_size=.false.
          return
       end if
    end if
    if (current_state%global_grid%active(X_INDEX)) then
       if (floor(real(current_state%global_grid%size(X_INDEX)) / x_dims) .lt. 2) then
          does_two_dim_work_with_domain_size=.false.
          return
       end if
    end if
    does_two_dim_work_with_domain_size=.true.
  end function does_two_dim_work_with_domain_size

  !> Applys the local bounds (start, end and size) for a specific dimension
  !! @param current_state The current model state_mod
  !! @param dim The dimension that we are applying
  !! @param dimSize The global size of the dimension
  !! @param dimProcesses Number of processes in the dimension
  !! @param dimMyRank My rank in the dimension
  subroutine apply_dimension_bounds(current_state, dim, dim_size, dim_processes, dim_my_rank)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, dim_size, dim_my_rank, dim_processes

    integer :: dimension_division, dimension_extra

    dimension_division = dim_size / dim_processes
    dimension_extra = dim_size - (dimension_division * dim_processes)

    current_state%local_grid%start(dim) = dimension_division*dim_my_rank + merge(dimension_extra, dim_my_rank, &
         dimension_extra .lt. dim_my_rank) + 1
    current_state%local_grid%end(dim) = (current_state%local_grid%start(dim)-1) + dimension_division + &
         merge(1, 0, dim_my_rank .lt. dimension_extra)
    current_state%local_grid%size(dim)=(current_state%local_grid%end(dim) - current_state%local_grid%start(dim)) + 1
    current_state%parallel%my_coords(dim) = dim_my_rank
  end subroutine apply_dimension_bounds

  !> Applys Z dimension information. As we always decompose into columns then Z is never
  !! split and as such is just the entire Z dimension regardless
  !! @param current_state The current model state_mod
  subroutine apply_z_dimension_information(current_state)
    type(model_state_type), intent(inout) :: current_state

    current_state%local_grid%start(Z_INDEX) = 1
    current_state%local_grid%end(Z_INDEX) = current_state%global_grid%size(Z_INDEX)
    current_state%local_grid%size(Z_INDEX) = current_state%global_grid%size(Z_INDEX)
    current_state%parallel%my_coords(Z_INDEX) = 0
    current_state%parallel%dim_sizes(Z_INDEX) = 1
  end subroutine apply_z_dimension_information

  !> One dimension decomposition
  !!
  !! Will decompose across the X or Y axis depending upon which is greater - it selects the greater
  !! to maximise the decomposition
  !! @param current_state The current model state_mod
  subroutine one_dim_decomposition(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: x_size, y_size, split_size, dimension_division, dimension_extra

    call apply_halo_information_and_allocate_neighbours(current_state)

    current_state%parallel%wrapped_around=.false.

    x_size = merge(current_state%global_grid%size(X_INDEX), 1, current_state%global_grid%active(X_INDEX))
    y_size = merge(current_state%global_grid%size(Y_INDEX), 1, current_state%global_grid%active(Y_INDEX))

    split_size = merge(x_size, y_size, x_size .gt. y_size)
    dimension_division = split_size / current_state%parallel%processes
    dimension_extra = split_size - (dimension_division * current_state%parallel%processes)

    current_state%local_grid%active = current_state%global_grid%active
    current_state%local_grid%dimensions = current_state%global_grid%dimensions
    current_state%parallel%neighbour_comm = current_state%parallel%monc_communicator
    call apply_z_dimension_information(current_state)         

    if (x_size .gt. y_size) then
      ! Decompose in X
      current_state%local_grid%start(Y_INDEX)=1
      current_state%local_grid%end(Y_INDEX)=y_size
      current_state%local_grid%size(Y_INDEX)=y_size
      current_state%parallel%my_coords(Y_INDEX)=0
      current_state%parallel%dim_sizes(Y_INDEX)=1
      current_state%local_grid%start(X_INDEX)=dimension_division*current_state%parallel%my_rank+merge(&
           dimension_extra, current_state%parallel%my_rank, dimension_extra .lt. current_state%parallel%my_rank) + 1
      current_state%local_grid%end(X_INDEX)=(current_state%local_grid%start(X_INDEX)-1) + dimension_division + merge(&
           1, 0, current_state%parallel%my_rank .lt. dimension_extra)
      current_state%local_grid%size(X_INDEX)=(current_state%local_grid%end(X_INDEX) - current_state%local_grid%start(X_INDEX)) + 1
      current_state%parallel%my_coords(X_INDEX) = current_state%parallel%my_rank
      current_state%parallel%dim_sizes(X_INDEX) = current_state%parallel%processes
      current_state%local_grid%neighbours(Y_INDEX,:) = current_state%parallel%my_rank
      ! Currently assume same PID in single direction (TODO: relax this restraint)
      current_state%local_grid%neighbours(X_INDEX,1:2) = merge(current_state%parallel%my_rank, current_state%parallel%processes, &
           current_state%parallel%my_rank .gt. 0) - 1
      current_state%parallel%wrapped_around(X_INDEX, 1)=&
           current_state%local_grid%neighbours(X_INDEX,1)==current_state%parallel%processes-1
      current_state%local_grid%neighbours(X_INDEX,3:4) = merge(current_state%parallel%my_rank+1, 0, &
           current_state%parallel%my_rank .lt. current_state%parallel%processes-1)
      current_state%parallel%wrapped_around(X_INDEX, 2)=current_state%local_grid%neighbours(X_INDEX,3)==0
      current_state%local_grid%corner_neighbours(1,:)=current_state%local_grid%neighbours(X_INDEX,1)
      current_state%local_grid%corner_neighbours(3,:)=current_state%local_grid%neighbours(X_INDEX,1)
      current_state%local_grid%corner_neighbours(2,:)=current_state%local_grid%neighbours(X_INDEX,3)
      current_state%local_grid%corner_neighbours(4,:)=current_state%local_grid%neighbours(X_INDEX,3)
      current_state%parallel%wrapped_around(Y_INDEX, :)=.true.
    else 
      ! Decompose in Y
      current_state%local_grid%start(X_INDEX)=1
      current_state%local_grid%end(X_INDEX)=x_size
      current_state%local_grid%size(X_INDEX)=x_size
      current_state%parallel%my_coords(X_INDEX)=0
      current_state%parallel%dim_sizes(X_INDEX)=1
      current_state%local_grid%start(Y_INDEX)=dimension_division*current_state%parallel%my_rank+merge(&
           dimension_extra, current_state%parallel%my_rank, dimension_extra .lt. current_state%parallel%my_rank) + 1
      current_state%local_grid%end(Y_INDEX)=(current_state%local_grid%start(Y_INDEX)-1) + dimension_division + merge(&
           1, 0, current_state%parallel%my_rank .lt. dimension_extra)
      current_state%local_grid%size(Y_INDEX)=(current_state%local_grid%end(Y_INDEX) - current_state%local_grid%start(Y_INDEX)) + 1
      current_state%parallel%my_coords(Y_INDEX)=current_state%parallel%my_rank
      current_state%parallel%dim_sizes(Y_INDEX) = current_state%parallel%processes
      current_state%local_grid%neighbours(X_INDEX,:) = current_state%parallel%my_rank
      current_state%local_grid%neighbours(Y_INDEX,1:2) = merge(current_state%parallel%my_rank, current_state%parallel%processes, &
           current_state%parallel%my_rank .gt. 0) - 1
      current_state%parallel%wrapped_around(Y_INDEX, 1)=&
           current_state%local_grid%neighbours(Y_INDEX,1)==current_state%parallel%processes-1
      current_state%local_grid%neighbours(Y_INDEX,3:4) = merge(current_state%parallel%my_rank+1, 0, &
           current_state%parallel%my_rank .lt. current_state%parallel%processes-1)
      current_state%parallel%wrapped_around(Y_INDEX, 2)=current_state%local_grid%neighbours(Y_INDEX,3)==0
      current_state%local_grid%corner_neighbours(1:2,:)=current_state%local_grid%neighbours(Y_INDEX,1)
      current_state%local_grid%corner_neighbours(3:4,:)=current_state%local_grid%neighbours(Y_INDEX,3)
      current_state%parallel%wrapped_around(X_INDEX, :)=.true.
    end if

    call apply_data_start_end_bounds(current_state%local_grid)
    call display_decomposition_information(current_state, "OneDim")
  end subroutine one_dim_decomposition

  !> Serial decomposition
  !!
  !! Simply places the entire domain onto the one process
  !! @param current_state The current model state_mod
  subroutine serial_decomposition(current_state)
    type(model_state_type), intent(inout) :: current_state

    call apply_halo_information_and_allocate_neighbours(current_state)

    current_state%local_grid%active = current_state%global_grid%active
    current_state%local_grid%dimensions = current_state%global_grid%dimensions
    current_state%parallel%neighbour_comm = current_state%parallel%monc_communicator
    call apply_z_dimension_information(current_state)

    current_state%local_grid%start(Y_INDEX)=1
    current_state%local_grid%start(X_INDEX)=1
    current_state%parallel%my_coords(Y_INDEX)=0
    current_state%parallel%my_coords(X_INDEX)=0
    current_state%parallel%dim_sizes(Y_INDEX)=0
    current_state%parallel%dim_sizes(X_INDEX)=0

    current_state%local_grid%end(X_INDEX)= &
         merge(current_state%global_grid%size(X_INDEX), 1, current_state%global_grid%active(X_INDEX))
    current_state%local_grid%size(X_INDEX)=current_state%local_grid%end(X_INDEX)

    current_state%local_grid%end(Y_INDEX)=&
         merge(current_state%global_grid%size(Y_INDEX), 1, current_state%global_grid%active(Y_INDEX))
    current_state%local_grid%size(Y_INDEX)=current_state%local_grid%end(Y_INDEX)

    current_state%local_grid%neighbours(X_INDEX,:) = current_state%parallel%my_rank
    current_state%local_grid%neighbours(Y_INDEX,:) = current_state%parallel%my_rank
    current_state%local_grid%corner_neighbours=current_state%parallel%my_rank

    current_state%parallel%wrapped_around=.true.

    call apply_data_start_end_bounds(current_state%local_grid)
    call display_decomposition_information(current_state, "Serial")
  end subroutine serial_decomposition

  !> Dumps out information about the decomposition at DEBUG level
  !! @param current_state The current model state_mod
  !! @param decompName The name of the decomposition which has been selected
  subroutine display_decomposition_information(current_state, decomp_name)
    type(model_state_type), intent(inout) :: current_state
    character(len=*), intent(in) :: decomp_name

    if (log_get_logging_level() .le. LOG_DEBUG) then
       call log_log(LOG_DEBUG, decomp_name//" rank="//trim(conv_to_string(current_state%parallel%my_rank))//&
            ": W=("//trim(conv_to_string(current_state%local_grid%start(Z_INDEX)))//"->"&
            //trim(conv_to_string(current_state%local_grid%end(Z_INDEX)))//") V=("//trim(conv_to_string(&
            current_state%local_grid%start(Y_INDEX)))//"->"//trim(conv_to_string(current_state%local_grid%end(&
            Y_INDEX)))//") U=("//trim(conv_to_string(current_state%local_grid%start(X_INDEX)))//&
            "->"//trim(conv_to_string(current_state%local_grid%end(X_INDEX)))//")")
       call log_log(LOG_DEBUG, "Neighbours of "//trim(conv_to_string(current_state%parallel%my_rank))//": y=(["//&
            trim(conv_to_string(current_state%local_grid%neighbours(Y_INDEX,1)))//"],["//&
            trim(conv_to_string(current_state%local_grid%neighbours(Y_INDEX,2)))//"]) x=(["//&
            trim(conv_to_string(current_state%local_grid%neighbours(X_INDEX,1)))//","//&
            trim(conv_to_string(current_state%local_grid%neighbours(X_INDEX,2)))//"],["//&
            trim(conv_to_string(current_state%local_grid%neighbours(X_INDEX,3)))//","//&
            trim(conv_to_string(current_state%local_grid%neighbours(X_INDEX,4)))//"])")
       call log_master_log(LOG_DEBUG, "Halo size z="//&
            trim(conv_to_string(current_state%local_grid%halo_size(Z_INDEX)))//&
            " y="//trim(conv_to_string(current_state%local_grid%halo_size(Y_INDEX)))//" x="//&
            trim(conv_to_string(current_state%local_grid%halo_size(X_INDEX))))
    end if
    call log_master_log(LOG_INFO, "Decomposed "//trim(conv_to_string(current_state%parallel%processes))//&
         " processes via '"//decomp_name// "' into z="//trim(conv_to_string(current_state%parallel%dim_sizes(&
         Z_INDEX)))//" y="//trim(conv_to_string(current_state%parallel%dim_sizes(Y_INDEX)))//" x="//&
         trim(conv_to_string(current_state%parallel%dim_sizes(X_INDEX))))
  end subroutine display_decomposition_information

  !> Calculates and applys the start and end bounds of the local data. Held locally is the halo and
  !! local data, this signifies for each dimension where the local data (without halo) starts and ends
  !! @param localGrid The populated local grid
  subroutine apply_data_start_end_bounds(local_grid)
    type(local_grid_type), intent(inout) :: local_grid

    integer :: i,n_dim

    !Total number of dimensions
    n_dim = 3

    do i = 1,n_dim
       local_grid%local_domain_start_index(i) = local_grid%halo_size(i) + 1
       local_grid%local_domain_end_index(i) = local_grid%halo_size(i) + local_grid%size(i)
    end do
    ! call log_master_log(LOG_DEBUG,"start_index(1)= "//&
    !      trim(conv_to_string(local_grid%local_domain_start_index(1)))//" end_index(1)="//&
    !      trim(conv_to_string(local_grid%local_domain_end_index(1)))//" start_index(2)="//&
    !      trim(conv_to_string(local_grid%local_domain_start_index(2)))//" end_index(2)="//&
    !      trim(conv_to_string(local_grid%local_domain_end_index(2)))//" start_index(3)="//&
    !      trim(conv_to_string(local_grid%local_domain_start_index(3)))//" end_index(3)="//&
    !      trim(conv_to_string(local_grid%local_domain_end_index(3))))

  end subroutine apply_data_start_end_bounds

  !> Will apply halo size information and allocated the neighbour array data.
  !! @param localGrid - The local grid
  subroutine apply_halo_information_and_allocate_neighbours(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: n_dim, n_corners,total_halo_size_XY_dim
    integer :: halo_depth

    ! Read halo_depth value
    halo_depth = options_get_integer(current_state%options_database, "halo_depth")

    ! There is no halo_depth in Z direction
    current_state%local_grid%halo_size(Z_INDEX) = 0
    current_state%local_grid%halo_size(Y_INDEX) = halo_depth
    current_state%local_grid%halo_size(X_INDEX) = halo_depth

    ! Left and right or Up and down
    total_halo_size_XY_dim  = current_state%local_grid%halo_size(X_INDEX)*2
    ! Total number of dimensions
    n_dim = 3
    ! Total number of corners
    n_corners = 4
    allocate(current_state%local_grid%neighbours(n_dim,total_halo_size_XY_dim), &
         current_state%local_grid%corner_neighbours(n_corners,halo_depth))
  end subroutine apply_halo_information_and_allocate_neighbours
end module decomposition_mod
