!> This is a thread pool and the single management "main" thread will spawn out free threads in the pool
!! to perform specific work. If there are no free threads then it will block until one becomes available. It uses 
!! ForThreads, which is a wrapper around pthreads. The thread pool works by creating a number of threads and then passing the
!! work to these threads, rather than creating a new thread for each piece of work. 
module threadpool_mod
  use forthread_mod, only : forthread_init, forthread_create, forthread_kill, forthread_destroy, forthread_mutex_init, &
       forthread_mutex_lock, forthread_mutex_unlock, forthread_mutex_destroy, forthread_cond_init, forthread_cond_wait, &
       forthread_cond_broadcast, forthread_cond_signal, forthread_cond_destroy, forthread_join
  use conversions_mod, only : conv_to_string
  use logging_mod, only : LOG_ERROR, LOG_WARN, log_log
  use configuration_parser_mod, only : io_configuration_type
  use mpi, only : MPI_THREAD_MULTIPLE, MPI_THREAD_SERIALIZED
  implicit none

#ifndef TEST_MODE
  private
#endif

  !> The thread call procedure interface
  interface
     subroutine thread_procedure(arguments, data_buffer)
       integer, dimension(:), intent(in) :: arguments
       character, dimension(:), allocatable, intent(inout), optional :: data_buffer
     end subroutine thread_procedure
  end interface

  !> Wraps the thread procedure with the call itself and the data to pass to it
  type threaded_procedure_container_type
     procedure(thread_procedure), pointer, nopass :: proc
     integer, dimension(:), allocatable :: arguments
     character, dimension(:), allocatable :: data_buffer
  end type threaded_procedure_container_type  

  integer, parameter :: DEFAULT_THREAD_POOL_SIZE=10 !< Number of threads in the pool
  logical, volatile, dimension(:), allocatable :: thread_busy, thread_start
  integer, volatile, dimension(:), allocatable :: thread_ids, thread_pass_data
  integer, volatile, dimension(:), allocatable :: activate_thread_condition_variables, activate_thread_mutex
  type(threaded_procedure_container_type), volatile, dimension(:), allocatable :: thread_entry_containers
  integer, volatile :: netcdfmutex !< Mutex used for controling NetCDF access
  integer, volatile :: next_suggested_idle_thread
  logical, volatile :: threadpool_active

  integer, volatile :: active_threads, total_number_of_threads, active_scalar_mutex

  public threadpool_init, threadpool_finalise, threadpool_start_thread, check_thread_status, threadpool_lock_netcdf_access, &
       threadpool_unlock_netcdf_access, threadpool_deactivate, threadpool_is_idle
contains

  !> Initialises the thread pool and marks each thread as idle
  subroutine threadpool_init(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    integer :: n

    call check_thread_status(forthread_init())
    call check_thread_status(forthread_mutex_init(netcdfmutex, -1))
    if (io_configuration%number_of_threads .ge. 1) then
      total_number_of_threads=io_configuration%number_of_threads
    else
      if (io_configuration%my_io_rank==0) then
        call log_log(LOG_WARN, "No setting for IO server thread pool size which must be 1 or more so using default size")
      end if
      total_number_of_threads=DEFAULT_THREAD_POOL_SIZE
    end if
    allocate(thread_busy(total_number_of_threads), thread_start(total_number_of_threads), &
         thread_ids(total_number_of_threads), thread_pass_data(total_number_of_threads), &
         activate_thread_condition_variables(total_number_of_threads), activate_thread_mutex(total_number_of_threads), &
         thread_entry_containers(total_number_of_threads))
    threadpool_active=.true.
    active_threads=total_number_of_threads
    next_suggested_idle_thread=1
    call check_thread_status(forthread_mutex_init(active_scalar_mutex, -1))
    do n=1, total_number_of_threads
      call check_thread_status(forthread_cond_init(activate_thread_condition_variables(n), -1))
      call check_thread_status(forthread_mutex_init(activate_thread_mutex(n), -1))
      thread_busy(n)=.false.
      thread_start(n)=.false.        
      thread_pass_data(n)=n
      call check_thread_status(forthread_create(thread_ids(n), -1, threadpool_thread_entry_procedure, thread_pass_data(n)))
    end do
  end subroutine threadpool_init

  !> Aquires the NetCDF thread lock, NetCDF is not thread safe so we need to manage thread calls to it
  subroutine threadpool_lock_netcdf_access()
#ifdef ENFORCE_THREAD_SAFETY
    call check_thread_status(forthread_mutex_lock(netcdfmutex))
#endif
  end subroutine threadpool_lock_netcdf_access

  !> Releases the NetCDF thread lock, NetCDF is not thread safe so we need to manage thread calls to it
  subroutine threadpool_unlock_netcdf_access()
#ifdef ENFORCE_THREAD_SAFETY
    call check_thread_status(forthread_mutex_unlock(netcdfmutex))
#endif
  end subroutine threadpool_unlock_netcdf_access

  !> Starts an idle thread from the pool to execute a specific procedure with some data. If there is no thread available
  !! then this will block until one becomes idle
  !! @param proc The procedure for the thread to execute
  !! @param data Data to pass into the thread
  subroutine threadpool_start_thread(proc, arguments, data_buffer)
    procedure(thread_procedure) :: proc
    integer, dimension(:), intent(in) :: arguments
    character, dimension(:), allocatable, intent(in), optional :: data_buffer

    integer :: idle_thread_id

    if (.not. threadpool_active) call log_log(LOG_ERROR, "Attemping to start IO thread on deactivated thread pool")

    idle_thread_id=find_idle_thread()
    if (idle_thread_id .ne. -1) then
      thread_busy(idle_thread_id)=.true.              
      thread_entry_containers(idle_thread_id)%proc=>proc
      allocate(thread_entry_containers(idle_thread_id)%arguments(size(arguments)))
      thread_entry_containers(idle_thread_id)%arguments=arguments
      if (present(data_buffer)) allocate(thread_entry_containers(idle_thread_id)%data_buffer(size(data_buffer)), &
           source=data_buffer)
      ! Send the signal to the thread to wake up and start
      call check_thread_status(forthread_mutex_lock(activate_thread_mutex(idle_thread_id)))
      thread_start(idle_thread_id)=.true.
      call check_thread_status(forthread_cond_signal(activate_thread_condition_variables(idle_thread_id)))
      call check_thread_status(forthread_mutex_unlock(activate_thread_mutex(idle_thread_id)))      
    end if
  end subroutine threadpool_start_thread

  !> Entry point called by each thread creation in the pool, this calls out to the actual procedure to execute and doing
  !! it this way allows us to perform some actions before or after which can help with the management of the pool
  !! @param thread_id The thread pool id (index) of this thread
  subroutine threadpool_thread_entry_procedure(thread_id)
    integer :: thread_id

    do while (threadpool_active)
      call check_thread_status(forthread_mutex_lock(activate_thread_mutex(thread_id)))
      do while (.not. thread_start(thread_id) .and. threadpool_active)
        call check_thread_status(forthread_cond_wait(activate_thread_condition_variables(thread_id), &
             activate_thread_mutex(thread_id)))
      end do
      call check_thread_status(forthread_mutex_unlock(activate_thread_mutex(thread_id)))      
      if (.not. threadpool_active) return
      thread_busy(thread_id)=.true.
      thread_start(thread_id)=.false.
      
      call check_thread_status(forthread_mutex_lock(active_scalar_mutex))
      active_threads=active_threads-1
      call check_thread_status(forthread_mutex_unlock(active_scalar_mutex))
      if (allocated(thread_entry_containers(thread_id)%data_buffer)) then
        call thread_entry_containers(thread_id)%proc(thread_entry_containers(thread_id)%arguments, &
             data_buffer=thread_entry_containers(thread_id)%data_buffer)
        deallocate(thread_entry_containers(thread_id)%data_buffer)
      else
        call thread_entry_containers(thread_id)%proc(thread_entry_containers(thread_id)%arguments)
      end if
      deallocate(thread_entry_containers(thread_id)%arguments)
      call check_thread_status(forthread_mutex_lock(active_scalar_mutex))
      active_threads=active_threads+1
      call check_thread_status(forthread_mutex_unlock(active_scalar_mutex))
      thread_busy(thread_id)=.false.
    end do
  end subroutine threadpool_thread_entry_procedure

  !> Determines whether the thread pool is idle or not (i.e. all threads are idle and waiting for work)
  !! @returns Whether the thread pool is idle
  logical function threadpool_is_idle()

    call check_thread_status(forthread_mutex_lock(active_scalar_mutex))
    threadpool_is_idle = active_threads==total_number_of_threads
    call check_thread_status(forthread_mutex_unlock(active_scalar_mutex))
  end function threadpool_is_idle  

  !> This waits for all busy threads to complete and then shuts all the pthreads down. The deactivation and finalisation
  !! procedures are split out as we want to deactivate the pool (to ensure no threads are running actions), finalise these
  !! actions which might involve destroying mutexes, and then destroying the threading environment in finalisation
  subroutine threadpool_deactivate()
    integer :: i
    integer, pointer :: retval

    allocate(retval)

    threadpool_active=.false.
    do i=1, total_number_of_threads
      call check_thread_status(forthread_mutex_lock(activate_thread_mutex(i)))
      call check_thread_status(forthread_cond_signal(activate_thread_condition_variables(i)))
      call check_thread_status(forthread_mutex_unlock(activate_thread_mutex(i)))
      call check_thread_status(forthread_join(thread_ids(i),retval))
      call check_thread_status(forthread_mutex_destroy(activate_thread_mutex(i)))
      call check_thread_status(forthread_cond_destroy(activate_thread_condition_variables(i)))
    end do
  end subroutine threadpool_deactivate  

  !> Finalises the thread pool
  subroutine threadpool_finalise()
      call check_thread_status(forthread_mutex_destroy(netcdfmutex))
      call check_thread_status(forthread_mutex_destroy(active_scalar_mutex))
      deallocate(thread_busy, thread_start, thread_ids, thread_pass_data, activate_thread_condition_variables, &
           activate_thread_mutex, thread_entry_containers)
    call check_thread_status(forthread_destroy())
  end subroutine threadpool_finalise

  !> Finds an idle thread, if one is not available then will block until one becomes free
  !! @returns The id of the idle thread which can be used
  integer function find_idle_thread()
    find_idle_thread=get_index_of_idle_thread()
    do while (find_idle_thread == -1)
      find_idle_thread=get_index_of_idle_thread()
    end do
  end function find_idle_thread
  
  !> Specifically gets the index of the next idle thread or -1 if they are all busy. This starts from a next suggested idle
  !! thread and will wrap around, as often the next thread will be idle rather than searching from the beginning again
  !! @returns The index of the next idle thread or -1 if there is none
  integer function get_index_of_idle_thread()
    integer :: i

    do i=next_suggested_idle_thread, total_number_of_threads
      if (.not. thread_busy(i)) then
        get_index_of_idle_thread=i
        next_suggested_idle_thread=i+1
        if (next_suggested_idle_thread .gt. total_number_of_threads) next_suggested_idle_thread=1
        return
      end if          
    end do
    next_suggested_idle_thread=1
    get_index_of_idle_thread=-1
  end function get_index_of_idle_thread

  !> Checks the error status of any thread operation and reports an error if it failed
  !! @param ierr The error/success flag returned from the ForThreads library, which itself is returned from pthreads
  subroutine check_thread_status(ierr)
    integer, intent(in) :: ierr

    if (ierr .ne. 0) then
      call log_log(LOG_ERROR, "Pthreads error in IO server, error code="//conv_to_string(ierr))
    end if    
  end subroutine check_thread_status  
end module threadpool_mod
