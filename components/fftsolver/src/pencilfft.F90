!> This Pencil FFT performs 3D forward and backwards FFTs using pencil decomposition. It uses FFTW for the actual FFT kernel
!! and this module contains all the data decomposition around this. There is no FFT required in Z, so this performs FFTs in
!! Y and X (in that order forward and reversed backwards.) The data decomposition is the complex aspect, there is the concept of
!! forward and backwards transformations. Forward transformations will go from pencil Z to Y to X and the backwards transformations
!! undo these, so go from X to Y to Z.
!! Note that we use quite a lot of buffer space here, this could be cut down if Y=X dimensions so some optimisation on memory
!! could be done there in that case
module pencil_fft_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX, global_grid_type
  use state_mod, only : model_state_type
  use fftw_mod, only : C_DOUBLE_COMPLEX, C_PTR, FFTW_BACKWARD, FFTW_FORWARD, FFTW_ESTIMATE, fftw_plan_many_dft_r2c, &
       fftw_plan_many_dft_c2r, fftw_execute_dft_c2r, fftw_execute_dft_r2c, fftw_destroy_plan
  use mpi, only : MPI_DOUBLE_COMPLEX, MPI_INT, MPI_COMM_SELF
  implicit none

#ifndef TEST_MODE
  private
#endif

  !> Describes a specific pencil transposition, from one pencil decomposition to another
  type pencil_transposition
     integer :: my_pencil_size(3), process_decomposition_layout(3), my_process_location(3), dim
     integer, dimension(:), allocatable :: send_sizes, send_offsets, recv_sizes, recv_offsets
     integer, dimension(:,:), allocatable :: recv_dims, send_dims
  end type pencil_transposition

  integer, parameter :: FORWARD=1, BACKWARD=2   !< Transposition directions
  integer :: dim_y_comm, dim_x_comm     !< Communicators for each dimension
  ! Transpositions from one pencil to another
  type(pencil_transposition) :: y_from_z_transposition, x_from_y_transposition, y_from_x_transposition, z_from_y_transposition, &
       y_from_z_2_transposition, x_from_y_2_transposition, y_from_x_2_transposition, z_from_y_2_transposition

  ! Temporary buffers used in transposition
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), contiguous, pointer :: real_buffer1, real_buffer2, real_buffer3, &
       fft_in_y_buffer , fft_in_x_buffer       
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:), contiguous, pointer :: buffer1, buffer2

  ! Pointers to FFTW plans and whether these have been initialised (only initialised once)
  type(C_PTR) :: fftw_plan(4)
  logical :: fftw_plan_initialised(4)=.false.
  
  public initialise_pencil_fft, finalise_pencil_fft, perform_forward_3dfft, perform_backwards_3dfft
contains

  !> Initialises the pencil FFT functionality, this will create the transposition structures needed
  !! @param current_state The current model state
  !! @param my_y_start My global start in fourier space for Y
  !! @param my_x_start My global start in fourier space for X
  !! @returns Size of local dimensions in fourier space for this process
  function initialise_pencil_fft(current_state, my_y_start, my_x_start)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(out) :: my_y_start, my_x_start
    integer :: initialise_pencil_fft(3) 

    integer :: ierr, y_distinct_sizes(current_state%parallel%dim_sizes(Y_INDEX)), &
         x_distinct_sizes(current_state%parallel%dim_sizes(X_INDEX))

    my_y_start=deduce_my_global_start(current_state, Y_INDEX)
    my_x_start=deduce_my_global_start(current_state, X_INDEX)

    if (current_state%parallel%dim_sizes(Y_INDEX) .gt. 1 .and. current_state%parallel%dim_sizes(X_INDEX) .gt. 1) then
      call mpi_cart_sub(current_state%parallel%neighbour_comm, (/1,0/), dim_y_comm, ierr)
      call mpi_cart_sub(current_state%parallel%neighbour_comm, (/0,1/), dim_x_comm, ierr)

      call mpi_allgather(current_state%local_grid%size(Y_INDEX), 1, MPI_INT, y_distinct_sizes, 1, MPI_INT, dim_y_comm, ierr)
      call mpi_allgather(current_state%local_grid%size(X_INDEX), 1, MPI_INT, x_distinct_sizes, 1, MPI_INT, dim_x_comm, ierr)
    else if (current_state%parallel%dim_sizes(Y_INDEX) .gt. 1) then
      dim_y_comm=current_state%parallel%monc_communicator
      dim_x_comm=MPI_COMM_SELF
      call mpi_allgather(current_state%local_grid%size(Y_INDEX), 1, MPI_INT, y_distinct_sizes, 1, MPI_INT, dim_y_comm, ierr)
      x_distinct_sizes=current_state%local_grid%size(X_INDEX)
    else if (current_state%parallel%dim_sizes(X_INDEX) .gt. 1) then      
      dim_y_comm=MPI_COMM_SELF
      dim_x_comm=current_state%parallel%monc_communicator
      y_distinct_sizes=current_state%local_grid%size(Y_INDEX)
      call mpi_allgather(current_state%local_grid%size(X_INDEX), 1, MPI_INT, x_distinct_sizes, 1, MPI_INT, dim_x_comm, ierr)
    else
      dim_y_comm=MPI_COMM_SELF
      dim_x_comm=MPI_COMM_SELF
      y_distinct_sizes=current_state%local_grid%size(Y_INDEX)
      x_distinct_sizes=current_state%local_grid%size(X_INDEX)
    end if

    call initialise_transpositions(current_state, y_distinct_sizes, x_distinct_sizes)
    call initialise_buffers()

    initialise_pencil_fft=z_from_y_transposition%my_pencil_size    
  end function initialise_pencil_fft  

  !> Cleans up allocated buffer memory
  subroutine finalise_pencil_fft(monc_communicator)
    integer, intent(in) :: monc_communicator
    integer :: ierr, i

    do i=1,size(fftw_plan_initialised)
      if (fftw_plan_initialised(i)) then
        call fftw_destroy_plan(fftw_plan(i))
      end if      
    end do    

    if (dim_y_comm .ne. MPI_COMM_SELF .and. dim_y_comm .ne. monc_communicator) call mpi_comm_free(dim_y_comm, ierr)
    if (dim_x_comm .ne. MPI_COMM_SELF .and. dim_x_comm .ne. monc_communicator) call mpi_comm_free(dim_x_comm, ierr)
    deallocate(buffer1, buffer2, real_buffer1, real_buffer2, real_buffer3, fft_in_y_buffer , fft_in_x_buffer)
  end subroutine finalise_pencil_fft  

  !> Performs a forward 3D FFT and currently results in target data which is the X, Z, Y oriented pencil
  !! Note that the source_data here takes no account for the halo, it is up to caller to exclude this.
  !! This does no FFT in Z, but transposes to Y, does FFT in Y, then transposes to X and 
  !! performs an FFT in that dimension. Pencil decomposition is used which has already been set up.
  !! @param current_state The current model state
  !! @param source_data The source real data to in the time domain
  !! @param target_data Frequency domain real representation of the time domain source which is allocated here
  subroutine perform_forward_3dfft(current_state, source_data, target_data)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(out) :: target_data

    call transpose_and_forward_fft_in_y(current_state, source_data, buffer1, real_buffer1)
    real_buffer1=real_buffer1/current_state%global_grid%size(Y_INDEX)
    call transpose_and_forward_fft_in_x(current_state, real_buffer1, buffer2, real_buffer2)
    real_buffer2=real_buffer2/current_state%global_grid%size(X_INDEX)

    call transpose_to_pencil(y_from_x_transposition, (/X_INDEX, Z_INDEX, Y_INDEX/), dim_x_comm, BACKWARD, &
         real_buffer2, real_buffer3)
    call transpose_to_pencil(z_from_y_transposition, (/Y_INDEX, X_INDEX, Z_INDEX/), dim_y_comm, BACKWARD, &
       real_buffer3, target_data)     
  end subroutine perform_forward_3dfft

  !> Performs a backwards 3D FFT and currently results in target data which is the X, Z, Y oriented pencil
  !! Note that the source_data here takes no account for the halo, it is up to caller to exclude this.
  !! This does no FFT in Z, but transposes to Y, does FFT in Y, then transposes to X and 
  !! performs an FFT in that dimension. Pencil decomposition is used which has already been set up.
  !! @param current_state The current model state
  !! @param source_data The source real data to in the frequency domain
  !! @param target_data Time domain complex representation of the frequency domain source
  subroutine perform_backwards_3dfft(current_state, source_data, target_data)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(out) :: target_data

    call transpose_to_pencil(y_from_z_2_transposition, (/Z_INDEX, Y_INDEX, X_INDEX/), dim_y_comm, FORWARD, &
       source_data, real_buffer3)
    call transpose_to_pencil(x_from_y_2_transposition, (/Y_INDEX, X_INDEX, Z_INDEX/), dim_x_comm, FORWARD, &
       real_buffer3, real_buffer2)

    call transpose_and_backward_fft_in_x(current_state, real_buffer2, buffer2, real_buffer1)
    call transpose_and_backward_fft_in_y(current_state, real_buffer1, buffer1, target_data)
  end subroutine perform_backwards_3dfft

  !> Initialises memory for the buffers used in the FFT
  subroutine initialise_buffers()
    allocate(buffer1(y_from_z_transposition%my_pencil_size(Y_INDEX)/2+1, y_from_z_transposition%my_pencil_size(X_INDEX), &
         y_from_z_transposition%my_pencil_size(Z_INDEX)), &
         real_buffer1((y_from_z_transposition%my_pencil_size(Y_INDEX)/2+1)*2, y_from_z_transposition%my_pencil_size(X_INDEX), &
         y_from_z_transposition%my_pencil_size(Z_INDEX)), &
         buffer2(x_from_y_transposition%my_pencil_size(X_INDEX)/2+1, x_from_y_transposition%my_pencil_size(Z_INDEX), &
         x_from_y_transposition%my_pencil_size(Y_INDEX)), &
         real_buffer2((x_from_y_transposition%my_pencil_size(X_INDEX)/2+1)*2, x_from_y_transposition%my_pencil_size(Z_INDEX), &
         x_from_y_transposition%my_pencil_size(Y_INDEX)), &
         fft_in_y_buffer(y_from_z_transposition%my_pencil_size(Y_INDEX), y_from_z_transposition%my_pencil_size(X_INDEX), &
         y_from_z_transposition%my_pencil_size(Z_INDEX)), &
         fft_in_x_buffer(x_from_y_transposition%my_pencil_size(X_INDEX), x_from_y_transposition%my_pencil_size(Z_INDEX), &
         x_from_y_transposition%my_pencil_size(Y_INDEX)), &
         real_buffer3(y_from_x_transposition%my_pencil_size(Y_INDEX), y_from_x_transposition%my_pencil_size(X_INDEX), &
         y_from_x_transposition%my_pencil_size(Z_INDEX)))
  end subroutine initialise_buffers

  !> Initialises the pencil transpositions, from a pencil in one dimension to that in another
  !! @param current_state The current model state
  !! @param y_distinct_sizes Y sizes per process 
  !! @param x_distinct_sizes X sizes per process
  subroutine initialise_transpositions(current_state, y_distinct_sizes, x_distinct_sizes)
    type(model_state_type), intent(inout) :: current_state
    integer, dimension(:), intent(in) :: y_distinct_sizes, x_distinct_sizes

    type(pencil_transposition) :: z_pencil

    z_pencil=create_initial_transposition_description(current_state)

    ! Transpositions
    y_from_z_transposition=create_transposition(current_state%global_grid, z_pencil, Y_INDEX, y_distinct_sizes, &
         FORWARD, (/ -1 /))   
    x_from_y_transposition=create_transposition(current_state%global_grid, y_from_z_transposition, X_INDEX, &
         x_distinct_sizes, FORWARD, (/ Y_INDEX /))
    y_from_x_transposition=create_transposition(current_state%global_grid, x_from_y_transposition, Y_INDEX, &
         normal_to_extended_process_dim_sizes(x_distinct_sizes), BACKWARD, (/ Y_INDEX, X_INDEX /))    
    z_from_y_transposition=create_transposition(current_state%global_grid, y_from_x_transposition, Z_INDEX, &
         normal_to_extended_process_dim_sizes(y_distinct_sizes), BACKWARD, (/ Y_INDEX, X_INDEX /))  

    y_from_z_2_transposition=create_transposition(current_state%global_grid, z_from_y_transposition, Y_INDEX, &
          normal_to_extended_process_dim_sizes(y_distinct_sizes), FORWARD, (/ Y_INDEX, X_INDEX /))   
    x_from_y_2_transposition=create_transposition(current_state%global_grid, y_from_z_2_transposition, X_INDEX, &
          normal_to_extended_process_dim_sizes(x_distinct_sizes), FORWARD, (/ Y_INDEX, X_INDEX /))
    y_from_x_2_transposition=create_transposition(current_state%global_grid, x_from_y_2_transposition, Y_INDEX, &
         x_distinct_sizes, BACKWARD, (/ Y_INDEX /))
    z_from_y_2_transposition=create_transposition(current_state%global_grid, y_from_x_2_transposition, Z_INDEX, &
          y_distinct_sizes, BACKWARD, (/ -1 /))   
  end subroutine initialise_transpositions

  !> Creates a specific pencil transposition description. It is maybe more a decomposition description, but the main
  !! complexity comes from the transposition from existing decomposition to new decomposition so therefore it is
  !! called transposition. The new pencil decomposition depends not only on the dimension to split on, but also the existing
  !! pencil decomposition. The new decomposed dimension (i.e. the existing pencil dimension) * other local dimensions is used
  !! as the sending size, receiving though requires knowledge about the data size on the source process so others will send
  !! their pencil dimension size to this process.
  !! @param new_pencil_dim The dimension to use as the new pencil decomposition
  !! @param existing_pencil_dim The dimension used in the current decomposition
  !! @param existing_pencil_process_layout Number of processes per dimension for the current decomposition
  !! @param existing_my_location The current processes block location per dimension for the current decomposition
  !! @param existing_pencil_size Pencil size per dimension for the current decomposition
  !! @param process_dim_sizes Sizes of the pencil dimension from other processes that is used to calculate receive count
  !! @param direction Whether we are transposing forwards or backwards, backwards is just an inverse
  !! @param extended_dimensions The dimensions that this process extends from n to (n/2+1)*2 (i.e. result of fft complex->real)
  type(pencil_transposition) function create_transposition(global_grid, existing_transposition, new_pencil_dim,&
       process_dim_sizes, direction, extended_dimensions)
    type(global_grid_type), intent(inout) :: global_grid
    type(pencil_transposition), intent(in) :: existing_transposition
    integer, dimension(:), intent(in) :: process_dim_sizes
    integer, intent(in) :: new_pencil_dim, direction, extended_dimensions(:)

    create_transposition%process_decomposition_layout=determine_pencil_process_dimensions(&
         new_pencil_dim, existing_transposition%dim, existing_transposition%process_decomposition_layout)

    create_transposition%my_process_location=determine_my_pencil_location(new_pencil_dim, &
         existing_transposition%dim, existing_transposition%my_process_location)

    create_transposition%my_pencil_size=determine_pencil_size(new_pencil_dim, create_transposition%process_decomposition_layout,&
         create_transposition%my_process_location, existing_transposition, global_grid, extended_dimensions)

    allocate(create_transposition%send_dims(3, create_transposition%process_decomposition_layout(existing_transposition%dim)), &
           create_transposition%recv_dims(3, create_transposition%process_decomposition_layout(existing_transposition%dim)))
    if (direction == FORWARD) then            
      call determine_my_process_sizes_per_dim(existing_transposition%dim, &
           existing_transposition%my_pencil_size, create_transposition%process_decomposition_layout, &
           global_grid, extended_dimensions, create_transposition%send_dims)
      call determine_matching_process_dimensions(new_pencil_dim, existing_transposition%dim, process_dim_sizes, &
           create_transposition%my_pencil_size, create_transposition%process_decomposition_layout, create_transposition%recv_dims)
    else
      call determine_my_process_sizes_per_dim(new_pencil_dim, create_transposition%my_pencil_size, &
           existing_transposition%process_decomposition_layout, global_grid, extended_dimensions, create_transposition%recv_dims)
      call determine_matching_process_dimensions(existing_transposition%dim, new_pencil_dim, process_dim_sizes, &
           existing_transposition%my_pencil_size, existing_transposition%process_decomposition_layout, &
           create_transposition%send_dims)
    end if

    allocate(create_transposition%send_sizes(size(create_transposition%send_dims, 2)), &
         create_transposition%send_offsets(size(create_transposition%send_sizes)), &
         create_transposition%recv_sizes(size(create_transposition%recv_dims, 2)), &
         create_transposition%recv_offsets(size(create_transposition%recv_sizes)))

    call concatenate_dimension_sizes(create_transposition%send_dims, create_transposition%send_sizes)
    call determine_offsets_from_size(create_transposition%send_sizes, create_transposition%send_offsets)

    call concatenate_dimension_sizes(create_transposition%recv_dims, create_transposition%recv_sizes)
    call determine_offsets_from_size(create_transposition%recv_sizes, create_transposition%recv_offsets)
    create_transposition%dim=new_pencil_dim
  end function create_transposition

  !> Performs the transposition and forward FFT in the y dimension then converts back to real numbers. The Y size is
  !! (n/2+1)*2 due to the complex to real transformation after the FFT.
  !! @param current_state The current model state
  !! @param source_data Input buffer, Z pencil oriented z,y,x
  !! @param buffer Complex buffer which the FFT writes into
  !! @param real_buffer Output buffer, Y pencil, oriented y,x,z
  subroutine transpose_and_forward_fft_in_y(current_state, source_data, buffer, real_buffer)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:),  intent(out) :: real_buffer
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:),  contiguous, pointer, intent(out) :: buffer

    ! Transpose globally from Z pencil to Y pencil
    call transpose_to_pencil(y_from_z_transposition, (/Z_INDEX, Y_INDEX, X_INDEX/), dim_y_comm, FORWARD, &
       source_data, fft_in_y_buffer)
    
    call perform_r2c_fft(fft_in_y_buffer, buffer, y_from_z_transposition%my_pencil_size(Y_INDEX), &
         y_from_z_transposition%my_pencil_size(X_INDEX) * y_from_z_transposition%my_pencil_size(Z_INDEX), 1)
    call convert_complex_to_real(buffer, real_buffer)
  end subroutine transpose_and_forward_fft_in_y

  !> Performs the backwards FFT in X and then transposes to Y pencil. The FFT requires complex numbers which are converted to real, 
  !! so the this real to complex operation is performed first. If n is the logical size of the FFT row, then the input
  !! size is n+2, complex number size is n/2+1 and we get n reals out.
  !! @param current_state The current model state
  !! @param source_data Input buffer, X pencil oriented x,z,y
  !! @param buffer Complex buffer which is fed into the FFT
  !! @param real_buffer Output buffer, Y pencil, oriented y,x,z
  subroutine transpose_and_backward_fft_in_x(current_state, source_data, buffer, real_buffer)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:),  intent(out) :: real_buffer
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), contiguous, pointer, intent(out) :: buffer

    call convert_real_to_complex(source_data, buffer)
    call perform_c2r_fft(buffer, fft_in_x_buffer, x_from_y_2_transposition%my_pencil_size(X_INDEX)-2, &
         x_from_y_2_transposition%my_pencil_size(Y_INDEX) * x_from_y_2_transposition%my_pencil_size(Z_INDEX), 2) 

    ! Transpose globally from X pencil to Y pencil
    call transpose_to_pencil(y_from_x_2_transposition, (/X_INDEX, Z_INDEX, Y_INDEX/), dim_x_comm, BACKWARD, &
       fft_in_x_buffer, real_buffer)
  end subroutine transpose_and_backward_fft_in_x

  !> Performs the transposition and forward FFT in the x dimension. After the FFT the complex space is converted back into
  !! real numbers. The X size is (n/2+1)*2 due to this transformation.
  !! @param current_state The current model state
  !! @param buffer1 Input buffer, Y pencil after the Y dimension FFT oriented y,x,z
  !! @param buffer Complex buffer which results from the FFT
  !! @param buffer2 Output buffer, X pencil after this X FFT, oriented x,z,y
  subroutine transpose_and_forward_fft_in_x(current_state, buffer1, buffer2, buffer3)
    type(model_state_type), target, intent(inout) :: current_state
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:),  contiguous, pointer, intent(out) :: buffer2
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: buffer1, buffer3

    ! Go from global Y pencil to global X pencil
    call transpose_to_pencil(x_from_y_transposition, (/Y_INDEX, X_INDEX, Z_INDEX/), dim_x_comm, FORWARD, &
       buffer1, fft_in_x_buffer)   

    call perform_r2c_fft(fft_in_x_buffer, buffer2, x_from_y_transposition%my_pencil_size(X_INDEX), &
         x_from_y_transposition%my_pencil_size(Y_INDEX) * x_from_y_transposition%my_pencil_size(Z_INDEX), 3)

    call convert_complex_to_real(buffer2, buffer3)
  end subroutine transpose_and_forward_fft_in_x

  !> Performs the backwards FFT in Y and then transposes to Z pencil. The FFT requires complex numbers which are converted to real, 
  !! so the this real to complex operation is performed first. If n is the logical size of the FFT row, then the input
  !! size is n+2, complex number size is n/2+1 and we get n reals out.
  !! @param current_state The current model state
  !! @param source_data Input buffer, Y pencil oriented y,x,z
  !! @param buffer Complex buffer which is fed into the FFT
  !! @param real_buffer Output buffer, Z pencil, oriented z,y,x
  subroutine transpose_and_backward_fft_in_y(current_state, source_data, buffer, real_buffer)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:),  intent(out) :: real_buffer
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), contiguous, pointer, intent(out) :: buffer

    call convert_real_to_complex(source_data, buffer)
   
    call perform_c2r_fft(buffer, fft_in_y_buffer,  y_from_x_2_transposition%my_pencil_size(Y_INDEX)-2, &
         y_from_x_2_transposition%my_pencil_size(X_INDEX) * y_from_x_2_transposition%my_pencil_size(Z_INDEX), 4)

    ! Go from global Y pencil to global Z pencil
    call transpose_to_pencil(z_from_y_2_transposition, (/Y_INDEX, X_INDEX, Z_INDEX/), dim_y_comm, BACKWARD, &
       fft_in_y_buffer, real_buffer)
  end subroutine transpose_and_backward_fft_in_y

  !> Transposes globally to a new pencil decomposition. This goes from the source dimensions a,b,c to b,c,a (forwards) or c,a,b
  !! (backwards.) It requires multiple steps, first the local data is transposed to c,b,a regardless of direction. 
  !! then it is communicated via alltoall, each process then assembles its own b,c,a or c,a,b data via contiguising 
  !! across blocks as the data layout is nonlinear.
  !! @param transposition_description Description of the transposition
  !! @param source_dims Dimensions of the current pencil that we wish to transpose from, will go from abc to bca
  !! @param communicator The MPI communicator associated with the group of processes who will swap data
  !! @param direction Whether this is going forwards or backwards, it makes a difference to the data arrangement
  !! @param source_data Source data (abc)
  !! @param target_data Target data (bca)
  subroutine transpose_to_pencil(transposition_description, source_dims, communicator, direction, source_data, target_data)
    type(pencil_transposition), intent(in) :: transposition_description
    integer, intent(in) :: source_dims(3), communicator, direction
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(out) :: target_data

    integer :: ierr
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: real_temp
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: real_temp2
    
    
    allocate(real_temp(size(source_data,3), size(source_data,2), size(source_data,1)), &
         real_temp2(product(transposition_description%my_pencil_size)+1))

    call rearrange_data_for_sending(real_source=source_data, real_target=real_temp)    

    call mpi_alltoallv(real_temp, transposition_description%send_sizes, transposition_description%send_offsets, &
         PRECISION_TYPE, real_temp2, transposition_description%recv_sizes, transposition_description%recv_offsets, &
         PRECISION_TYPE, communicator, ierr)
    call contiguise_data(transposition_description, (/source_dims(3), source_dims(2), source_dims(1)/), direction, &
         source_real_buffer=real_temp2, target_real_buffer=target_data)
    deallocate(real_temp, real_temp2)   
  end subroutine transpose_to_pencil  

  !> Contiguises from c,b,a to b,c,a (forwards) or c,a,b (backwards) where these are defined by the source_dims argument. 
  !! It is not as simple as just swapping the required dimensions, as this is after the mpi alltoall and each block lies 
  !! after the previous block running sequentially in a.
  !! @param transposition_description Transposition descriptor
  !! @param source_dims Representation a,b,c of source data, will contiguise to b,a,c
  !! @param direction Whether we wish to contiguise forwards or backwards
  !! @param source_real_buffer Source real data to transform
  !! @param target_real_buffer Target real data which is the result of the operation
  subroutine contiguise_data(transposition_description, source_dims, direction, source_real_buffer, target_real_buffer)
    integer, intent(in) :: source_dims(3), direction
    type(pencil_transposition), intent(in) :: transposition_description
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: source_real_buffer
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(out) :: target_real_buffer
    
    integer :: number_blocks, i, j, k, n, index_prefix, index_prefix_dim, block_offset, source_index

    number_blocks=size(transposition_description%recv_sizes)
    index_prefix=0
    block_offset=0
    index_prefix_dim=merge(2,1, direction == FORWARD)
    do i=1,number_blocks
      if (i .ge. 2) then
        index_prefix=index_prefix+transposition_description%recv_dims(source_dims(index_prefix_dim), i-1)
        block_offset=block_offset+transposition_description%recv_sizes(i-1)
      end if
      !Transformation is either cba -> bca (forward) or cab (backwards)       
      do j=1, transposition_description%recv_dims(source_dims(3), i) ! a
        do k=1, transposition_description%recv_dims(source_dims(1), i) ! c
          do n=1, transposition_description%recv_dims(source_dims(2), i) ! b
            source_index=block_offset+(j-1)* transposition_description%recv_dims(source_dims(1), i)* &
                 transposition_description%recv_dims(source_dims(2), i)+ (n-1)* &
                 transposition_description%recv_dims(source_dims(1), i)+k            
            if (direction == FORWARD) then
              target_real_buffer(index_prefix+n, k, j)=source_real_buffer(source_index) ! bca
            else
              target_real_buffer(index_prefix+k, j, n)=source_real_buffer(source_index) ! cab
            end if
          end do
        end do
      end do      
    end do    
  end subroutine contiguise_data

  !> Actually performs a forward real to complex FFT
  !! @param source_data Source (real) data in the time domain
  !! @param transformed_data Resulting complex data in the frequency domain
  !! @param row_size Number of elements for each FFT
  !! @param num_rows The number of FFTs to perform on the next data elements in the source_data
  !! @param plan_id Id number of the plan that tracks whether we need to create it or can reuse the existing one
  subroutine perform_r2c_fft(source_data, transformed_data, row_size, num_rows, plan_id)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), contiguous, pointer, intent(inout) :: source_data
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), contiguous, pointer, intent(inout) :: transformed_data
    integer, intent(in) :: row_size, num_rows, plan_id
    
    if (.not. fftw_plan_initialised(plan_id)) then
      fftw_plan(plan_id) = fftw_plan_many_dft_r2c(1, (/row_size/), num_rows, source_data, (/row_size/), 1, row_size, &
           transformed_data, (/row_size/), 1, row_size/2+1, FFTW_ESTIMATE)
      fftw_plan_initialised(plan_id)=.true.
    end if
    call fftw_execute_dft_r2c(fftw_plan(plan_id), source_data, transformed_data)
  end subroutine perform_r2c_fft

  !> Performs the complex to real (backwards) FFT
  !! @param source_data Source (complex) data in the frequency domain
  !! @param transformed_data Resulting real data in the time domain
  !! @param row_size Number of elements for each FFT
  !! @param num_rows The number of FFTs to perform on the next data elements in the source_data
  !! @param plan_id Id number of the plan that tracks whether we need to create it or can reuse the existing one
  subroutine perform_c2r_fft(source_data, transformed_data, row_size, num_rows, plan_id)
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), contiguous, pointer, intent(inout) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), contiguous, pointer, intent(inout) :: transformed_data
    integer, intent(in) :: row_size, num_rows, plan_id   

    if (.not. fftw_plan_initialised(plan_id)) then
      ! n is the size of the FFT (in real, not complex->real coords.) There are row_size/2+1 between entries for the input
      ! (complex) data and row_size between entries for the output data
      fftw_plan(plan_id) = fftw_plan_many_dft_c2r(1, (/row_size/), num_rows, source_data, (/row_size/2+1/), 1, row_size/2+1, &
           transformed_data, (/row_size/), 1, row_size, FFTW_ESTIMATE)
      fftw_plan_initialised(plan_id)=.true.
    end if
    call fftw_execute_dft_c2r(fftw_plan(plan_id), source_data, transformed_data)
  end subroutine perform_c2r_fft

  !> Rearranges data for sending, transposing a,b,c into c,b,a . This is done as alltoall splits on dimension c
  !! so to go from one pencil to another we assume here that a is the existing pencil as it is contiguous
  !! @param real_source Source data to transpose from
  !! @param real_target Target data to transpose to
  subroutine rearrange_data_for_sending(real_source, real_target)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: real_source
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(out) :: real_target

    integer :: i

    do i=1, size(real_source,2)
      real_target(:,i,:)=transpose(real_source(:,i,:))
    end do
  end subroutine rearrange_data_for_sending

  !> Determines the number of elements to on my process per dimension which either need to be sent to (forwards transformation) or
  !! received from (backwards) each target process (in the row or column)
  !! This depends on the existing pencil decomposition, as effectively we are breaking that contigulity and
  !! decomposing it into n blocks in that dimension now (provided by new_pencil_procs_per_dim)
  !! @param existing_pencil_dim The pencil dimension that we are transforming from
  !! @param existing_pencil_size Existing pencil decomposition sizes per dimension
  !! @param new_pencil_procs_per_dim For the target decomposition the number of processes per dimension
  !! @param global_grid Description of the global grid which we use for sizing information
  !! @param extended_dimensions List of dimensions where we extend from n to n+2 (i.e. result of FFT complex-> real transformation)
  subroutine determine_my_process_sizes_per_dim(existing_pencil_dim, existing_pencil_size, new_pencil_procs_per_dim, &
       global_grid, extended_dimensions, specific_sizes_per_dim)
    integer, intent(in) :: existing_pencil_dim, existing_pencil_size(:), new_pencil_procs_per_dim(:), extended_dimensions(:)
    type(global_grid_type), intent(inout) :: global_grid
    integer, dimension(:,:), intent(inout) :: specific_sizes_per_dim

    integer :: i, split_size, split_remainder, j, s

    do i=1,3
      if (i == existing_pencil_dim) then
        s=global_grid%size(i)
        if (is_extended_dimension(i, extended_dimensions)) s=s+2
        split_size = s / new_pencil_procs_per_dim(i)
        split_remainder = s - split_size * new_pencil_procs_per_dim(i)
        do j=1,new_pencil_procs_per_dim(existing_pencil_dim)
          specific_sizes_per_dim(i,j)=merge(split_size+1, split_size, j .le. split_remainder)
        end do        
      else
        specific_sizes_per_dim(i,:) = existing_pencil_size(i)
      end if
    end do    
  end subroutine determine_my_process_sizes_per_dim

  !> Simple helper function to deduce send or receive offsets from the sizes
  !! @param source_sizes Sizes that we are using to build the offsets
  subroutine determine_offsets_from_size(source_sizes, determined_offsets)
    integer, intent(in) :: source_sizes(:)
    integer, dimension(:), intent(inout) :: determined_offsets

    integer :: i

    determined_offsets(1)=0
    do i=2,size(source_sizes)
      determined_offsets(i)=determined_offsets(i-1)+source_sizes(i-1)
    end do    
  end subroutine determine_offsets_from_size  

  !> Determines the number of processes in each dimension for the target decomposition. This depends heavily
  !! on the existing decomposition, as we basically contiguise our pencil dimension and decompose the existing
  !! pencil dimension. The third dimension remains unchanged
  !! @param new_pencil_dim New pencil dimension
  !! @param existing_pencil_dim Current decomposition pencil dimension
  !! @param existing_pencil_procs Current decomposition process layout
  function determine_pencil_process_dimensions(new_pencil_dim, existing_pencil_dim, existing_pencil_procs)
    integer, intent(in) :: new_pencil_dim, existing_pencil_dim, existing_pencil_procs(3)
    integer :: determine_pencil_process_dimensions(3)

    integer :: i

    do i=1,3
      if (i == new_pencil_dim) then
        determine_pencil_process_dimensions(i)=1
      else if (i == existing_pencil_dim) then
        determine_pencil_process_dimensions(i)=existing_pencil_procs(new_pencil_dim)
      else
        determine_pencil_process_dimensions(i)=existing_pencil_procs(i)
      end if
    end do    
  end function determine_pencil_process_dimensions

  !> Determines my location for each dimension in the new pencil decomposition. I.e. which block I am
  !! operating on
  !! @param new_pencil_dim New pencil decomposition dimension
  !! @param existing_pencil_dim Current pencil dimension
  !! @param existing_locations Location for the current decomposition
  function determine_my_pencil_location(new_pencil_dim, existing_pencil_dim, existing_locations)
    integer, intent(in) :: new_pencil_dim, existing_pencil_dim, existing_locations(3)
    integer :: determine_my_pencil_location(3)

    integer :: i

    do i=1,3
      if (i == new_pencil_dim) then
        determine_my_pencil_location(i)=1
      else if (i == existing_pencil_dim) then
        determine_my_pencil_location(i)=existing_locations(new_pencil_dim)
      else
        determine_my_pencil_location(i)=existing_locations(i)
      end if
    end do    
  end function determine_my_pencil_location  

  !> Concatenates sizes in multiple dimensions for each target process (in a row or column) into a product of
  !! that. This represents all the dimension sizes per process
  !! @param dims The sizes, per dimension and per process that we will fold into target process
  subroutine concatenate_dimension_sizes(dims, concatenated_dim_sizes)
    integer, dimension(:,:), intent(in) :: dims
    integer, dimension(:), intent(inout) :: concatenated_dim_sizes

    integer :: i

    do i=1,size(dims, 2)
      concatenated_dim_sizes(i)=product(dims(:,i))
    end do    
  end subroutine concatenate_dimension_sizes  

  !> Determines the sizes per dimension on the matching process either to receive from (forward transposition) or send to 
  !! (backwards transposition) each source process. Not only does this depend on the
  !! my pencil sizes, but it also depends on the amount of data that the source process has to send over
  !! @param new_pencil_dim The dimension for the new pencil decomposition
  !! @param existing_pencil_dim Dimension for the existing pencil decomposition
  !! @param proc_sizes Size of dimension on the source processes (index in array corresponds to source PID)
  !! @param my_pencil_size My (new) pencil size per dimension
  !! @param pencil_processes_per_dim The process layout per dimension
  subroutine determine_matching_process_dimensions(new_pencil_dim, existing_pencil_dim, proc_sizes, &
       my_pencil_size, pencil_processes_per_dim, specific_sizes_per_dim)
    integer, intent(in) :: new_pencil_dim, existing_pencil_dim, proc_sizes(:), my_pencil_size(:), pencil_processes_per_dim(:)
    integer, dimension(:,:), intent(inout) :: specific_sizes_per_dim

    integer :: i, j

    do i=1,pencil_processes_per_dim(existing_pencil_dim)
      do j=1,3
        if (j==new_pencil_dim) then
          specific_sizes_per_dim(j, i)=proc_sizes(i)
        else
          specific_sizes_per_dim(j, i)=my_pencil_size(j)
        end if
      end do      
    end do    
  end subroutine determine_matching_process_dimensions

  !> Creates an initial transposition representation of the Z pencil that MONC is normally decomposed in. This is then
  !! fed into the create transposition procedure which will generate transpositions to other pencils
  !! @param current_state The current model state
  type(pencil_transposition) function create_initial_transposition_description(current_state)
    type(model_state_type), intent(inout) :: current_state

    create_initial_transposition_description%dim=Z_INDEX
    create_initial_transposition_description%process_decomposition_layout=current_state%parallel%dim_sizes
    create_initial_transposition_description%my_process_location=current_state%parallel%my_coords
    create_initial_transposition_description%my_pencil_size=current_state%local_grid%size
  end function create_initial_transposition_description

  !> Deduces the size of my (local) pencil based upon the new decomposition. This depends heavily on the current
  !! pencil decomposition, the new pencil dimension is the global size, the existing pencil dimension becomes 
  !! decomposed based on the number of processes in that dimension. The third dimension remains unchanged.
  !! @param new_pencil_dim Dimension for the new pencil decomposition
  !! @param pencil_process_layout The processes per dimension layout for the new decomposition
  !! @param my_pencil_location My location in the block layout
  !! @param existing_pencil_dim Current decomposition dimension
  !! @param existing_pencil_size Current decomposition sizes
  !! @param global_grid Description of the global grid which we use for sizing information
  !! @param extended_dimensions List of dimensions where we extend from n to n+2 (i.e. result of FFT complex-> real transformation)
  function determine_pencil_size(new_pencil_dim, pencil_process_layout, my_pencil_location, existing_transposition,&
       global_grid, extended_dimensions)

    type(pencil_transposition), intent(in) :: existing_transposition
    integer, intent(in) :: new_pencil_dim, pencil_process_layout(3), my_pencil_location(3), extended_dimensions(:)
    type(global_grid_type), intent(inout) :: global_grid
    integer :: determine_pencil_size(3)

    integer :: i, split_size, split_remainder, s

    do i=1,3
      if (i == new_pencil_dim) then
        if (is_extended_dimension(i, extended_dimensions)) then
          ! If complex and Y dim then /2+1 for the global size
          determine_pencil_size(i)=(global_grid%size(new_pencil_dim)/2+1)*2
        else
          determine_pencil_size(i)=global_grid%size(new_pencil_dim)
        end if
      else if (i == existing_transposition%dim) then
        s=global_grid%size(i)
        ! If complex and Y dim then use s/2+1 for the size to split
        if (is_extended_dimension(i, extended_dimensions)) s=(s/2+1)*2
        split_size=s/pencil_process_layout(i)
        split_remainder=s - split_size * pencil_process_layout(i)
        determine_pencil_size(i)=merge(split_size+1, split_size, my_pencil_location(i)+1 .le. split_remainder)        
      else        
        determine_pencil_size(i)=existing_transposition%my_pencil_size(i)       
      end if      
    end do    
  end function determine_pencil_size

  !> Determines whether or not the specific dimension is in the list of extended dimensions
  !! @param dimension The dimension to test for
  !! @param extended_dimensions Array of dimensions that will be searched
  !! @returns Whether the dimension is found in the array
  logical function is_extended_dimension(dimension, extended_dimensions)
    integer, intent(in) :: dimension, extended_dimensions(:)

    integer :: i
    do i=1,size(extended_dimensions)
      if (extended_dimensions(i) == dimension) then
        is_extended_dimension=.true.
        return
      end if      
    end do
    is_extended_dimension=.false.
  end function is_extended_dimension
  
  !> Transforms real process dimension sizes into their real after FFT complex->real transformation. The way this works is that
  !! it goes from n to (n/2+1)*2 numbers which is distributed amongst the processes deterministically
  !! @param process_dim_sizes Real process dimension sizes
  !! @returns The extended process dimension sizes
  function normal_to_extended_process_dim_sizes(process_dim_sizes)
    integer, dimension(:), intent(in) :: process_dim_sizes
    integer, dimension(size(process_dim_sizes)) :: normal_to_extended_process_dim_sizes

    integer :: temp_total, split_size, remainder

    temp_total=(sum(process_dim_sizes) /2 + 1) * 2
    split_size=temp_total/size(process_dim_sizes)
    remainder=temp_total - split_size*size(process_dim_sizes)

    normal_to_extended_process_dim_sizes=split_size
    normal_to_extended_process_dim_sizes(1:remainder)=split_size+1    
  end function normal_to_extended_process_dim_sizes

  !> Converts complex representation to its real data counterpart and is called after each forward FFT. 
  !! After a r2c FFT, there are n/2+1 complex numbers - which means that there will be more real numbers in Fourier space
  !! than are provided into the forward FFT call (due to the extra +1). Note that the real size n will always be complex size * 2
  !! This always unpacks the complex dimension in the first dimension
  !! @param complex_data Complex data in Z,Y,X orientation to be unpacked into its real representation
  !! @param real_data The real representation is written into here
  subroutine convert_complex_to_real(complex_data, real_data)
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), intent(in) :: complex_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(out) :: real_data

    integer :: i, j, k

    do i=1,size(real_data,3)
      do j=1,size(real_data,2)
        do k=1,size(real_data,1),2
          real_data(k,j,i)=real(real(complex_data((k+1)/2,j,i)), kind=DEFAULT_PRECISION)
          real_data(k+1,j,i)=real(aimag(complex_data((k+1)/2,j,i)), kind=DEFAULT_PRECISION)
        end do
      end do
    end do
  end subroutine convert_complex_to_real

  !> Converts reals into their complex representation, this is called for backwards FFTs as we need to feed in complex numbers
  !! to force FFTW to do a backwards. It is a relatively simple transformation, as n goes into n/2 complex numbers and as this
  !! is the result of the `convert_complex_to_real` procedure, n always divides evenly.
  !! This is always applied to the first dimension of the real data
  !! @param real_data The source real data to pack into the complex data, it is oriented Z,Y,X
  !! @param complex_data Target complex data which the real data is packaged into
  subroutine convert_real_to_complex(real_data, complex_data)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: real_data
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), contiguous, pointer, intent(out) :: complex_data
    
    integer :: i, j, k

    complex_data=cmplx(0.0d0, 0.0d0, kind=C_DOUBLE_COMPLEX)

    do i=1,size(real_data,3)
      do j=1,size(real_data,2)
        do k=1,size(real_data,1),2
          complex_data((k+1)/2,j,i)=cmplx(real_data(k,j,i), real_data(k+1,j,i), kind=C_DOUBLE_COMPLEX)
        end do
      end do
    end do
  end subroutine convert_real_to_complex

  !> Determines my global start coordinate in Fourier space. 
  !! This is required for cos y and cos x calculation which is fed into the tridiagonal solver. After the forward FFTs, 
  !! each process has ((n/2+1)/p+r) * 2 elements, where p is the number of processes and r
  !! is the uneven process remainder (1 or 0 depending on p). Therefore some processes will have t elements, and some t-2 elements
  !! to feed into the solver
  !! @param current_state The current model state
  !! @param dimension The dimension that we are calculating this for (Y or X)
  !! @returns My global start in Fourier space
  integer function deduce_my_global_start(current_state, dimension)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dimension

    integer complex_size, distributed_size, remainder, larger_nums, smaller_nums

    complex_size=(current_state%global_grid%size(dimension)/2+1)*2
    distributed_size=complex_size / current_state%parallel%dim_sizes(dimension)
    remainder=complex_size - distributed_size * current_state%parallel%dim_sizes(dimension) 
    larger_nums=min(remainder, current_state%parallel%my_coords(dimension))
    smaller_nums=current_state%parallel%my_coords(dimension)-remainder
    deduce_my_global_start=((distributed_size+1)*larger_nums + merge(distributed_size*smaller_nums, 0, smaller_nums .gt. 0)) + 1
  end function deduce_my_global_start
end module pencil_fft_mod    
    
