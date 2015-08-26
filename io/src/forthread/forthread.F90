module forthread_mod
  use iso_c_binding
  use forthread_data
  use forthread_types
  use forthread_ciface_mod
  implicit none

contains

  integer function forthread_init()
    integer :: info

    allocate(routine_table(init_size))
    routine_table_size = init_size

    call thread_init(info)

    call thread_mutex_init(routine_table_mutex,-1,info)
    forthread_init=info
  end function forthread_init

  integer function forthread_destroy()
    integer :: info

    deallocate(routine_table)
    routine_table_size = 0
    call thread_mutex_destroy(routine_table_mutex,info)
    call thread_destroy(info)
    forthread_destroy=info
  end function forthread_destroy

  integer function forthread_create(thread_id,attr_id,run,arg)
    integer, intent(out) :: thread_id
    integer, intent(in) :: attr_id
    procedure(i_run) :: run !type i_run
    integer, target  :: arg

    integer :: i, info
    procedure(i_start_routine), bind(c), pointer :: start_routinep
    type(ptr_t_run), dimension(:), pointer       :: tmp
    type(t_run), pointer :: runp
    call thread_mutex_lock(routine_table_mutex,info)

    call thread_alloc(thread_id,info)
    if (thread_id.gt.routine_table_size) then
      nullify(tmp)
      allocate(tmp(routine_table_size*2))
      do i=1,routine_table_size
        tmp(i) = routine_table(i)
      enddo
      deallocate(routine_table)
      routine_table => tmp
      routine_table_size = routine_table_size*2
    endif
    allocate(routine_table(thread_id)%t)
    routine_table(thread_id)%t%run => run
    routine_table(thread_id)%t%arg => arg
    start_routinep => start_routine

    call thread_create(thread_id,attr_id,c_funloc(start_routinep),&
         c_loc(routine_table(thread_id)%t),info)

    call thread_mutex_unlock(routine_table_mutex,info)
    forthread_create=info
  end function forthread_create

  integer function forthread_detach(thread_id)
    integer, intent(in) :: thread_id

    integer :: info

    call thread_detach(thread_id,info)
    forthread_detach=info
  end function forthread_detach

  integer function forthread_equal(t1,t2)
    integer, intent(in) :: t1, t2

    integer :: info

    call thread_equal(t1,t2,info)
    forthread_equal=info
  end function forthread_equal

  ! Exits the current thread
  subroutine forthread_exit(val)
    integer, pointer :: val

    call thread_exit(c_loc(val))
  end subroutine forthread_exit

  integer function forthread_join(thread_id,val)
    integer, intent(in) :: thread_id
    integer, pointer:: val

    integer :: info

    type(c_ptr)                 :: value_ptr
    call thread_join(thread_id,value_ptr,info)
    call c_f_pointer(value_ptr,val)
    forthread_join=info
  end function forthread_join

  integer function forthread_cancel(thread_id)
    integer, intent(in) :: thread_id

    integer :: info

    call thread_cancel(thread_id,info)
    forthread_cancel=info
  end function forthread_cancel

  integer function forthread_kill(thread_id,sig)
    integer, intent(in) :: thread_id, sig

    integer :: info

    call thread_kill(thread_id,sig,info)
    forthread_kill=info
  end function forthread_kill

  integer function forthread_once_init(once_ctrl_id)
    integer, intent(out) :: once_ctrl_id

    integer :: info

    call thread_once_init(once_ctrl_id,info)
    forthread_once_init=info
  end function forthread_once_init

  integer function forthread_once(once_ctrl_id,init_routine)
    integer, intent(in) :: once_ctrl_id
    procedure(i_once) :: init_routine
    ! dangerous but works! (gfortran)
    ! TODO test in other compilers

    integer :: info

    call thread_once(once_ctrl_id,c_funloc(init_routine),info)
    forthread_once=info
  end function forthread_once

  ! TODO implement thread_atfork

  integer function forthread_getconcurrency(currlevel)
    integer, intent(out) :: currlevel

    integer :: info

    call thread_getconcurrency(currlevel,info)
    forthread_getconcurrency=info
  end function forthread_getconcurrency

  integer function forthread_setconcurrency(newlevel)
    integer, intent(in) :: newlevel

    integer :: info

    call thread_setconcurrency(newlevel,info)
    forthread_setconcurrency=info
  end function forthread_setconcurrency

#ifndef __DARWIN
  integer function forthread_getcpuclockid(thread,clock_id)
    integer, intent(in) :: thread
    integer, intent(out) :: clock_id

    integer :: info

    call thread_getcpuclockid(thread,clock_id,info)
    forthread_getcpuclockid=info
  end function forthread_getcpuclockid
#endif

  integer function forthread_getschedparam(thread,policy,param)
    integer, intent(in) :: thread
    integer, intent(out) :: policy
    type(sched_param), intent(out) :: param

    integer :: info

    call thread_getschedparam(thread,policy,param,info)
    forthread_getschedparam=info
  end function forthread_getschedparam

  integer function forthread_setschedparam(thread,policy,param)
    integer, intent(in) :: thread, policy
    type(sched_param), intent(in) :: param

    integer :: info

    call thread_setschedparam(thread,policy,param,info)
    forthread_setschedparam=info
  end function forthread_setschedparam

#ifndef __DARWIN
  integer function forthread_setschedprio(thread,prio)
    integer, intent(in) :: thread, prio

    integer :: info

    call thread_setschedprio(thread,prio,info)
    forthread_setschedprio=info
  end function forthread_setschedprio
#endif

  integer function forthread_setcancelstate(state,oldstate)
    integer, intent(in) :: state
    integer, intent(out) :: oldstate

    integer :: info

    call thread_setcancelstate(state,oldstate,info)
    forthread_setcancelstate=info
  end function forthread_setcancelstate

  integer function forthread_setcanceltype(ctype,oldctype)
    integer, intent(in) :: ctype
    integer, intent(out) :: oldctype

    integer :: info

    call thread_setcanceltype(ctype,oldctype,info)
    forthread_setcanceltype=info
  end function forthread_setcanceltype

  !*****************************************!
  !*   sharing private data in threads     *!
  !*****************************************!

  integer function forthread_key_delete(key_id)
    integer, intent(in) :: key_id

    integer :: info

    call thread_key_delete(key_id,info)
    forthread_key_delete=info
  end function forthread_key_delete

  integer function forthread_key_create(key_id,destructor)
    integer, intent(out) :: key_id
    procedure(i_destructor) :: destructor
    ! dangerous but works! (gfortran)
    ! TODO test in other compilers

    integer :: info

    call thread_key_create(key_id,c_funloc(destructor),info)
    forthread_key_create=info
  end function forthread_key_create

  ! no wrappers provided for the following two routines
  !void thread_getspecific(int *key, void **value, int *info);

  !void thread_setspecific(int *key, void **value, int *info);



  !*****************************************!
  !*             mutex routines            *!
  !*****************************************!


  integer function forthread_mutex_destroy(mutex_id)
    integer, intent(in) :: mutex_id

    integer :: info

    call thread_mutex_destroy(mutex_id,info)
    forthread_mutex_destroy=info
  end function forthread_mutex_destroy

  integer function forthread_mutex_init(mutex_id,attr_id)
    integer, intent(out) :: mutex_id
    integer, intent(in) :: attr_id

    integer :: info

    call thread_mutex_init(mutex_id,attr_id,info)
    forthread_mutex_init=info
  end function forthread_mutex_init

  integer function forthread_mutex_lock(mutex_id)
    integer, intent(in) :: mutex_id

    integer :: info

    call thread_mutex_lock(mutex_id,info)
    forthread_mutex_lock=info
  end function forthread_mutex_lock

  integer function forthread_mutex_trylock(mutex_id)
    integer, intent(in) :: mutex_id

    integer :: info

    call thread_mutex_trylock(mutex_id,info)
    forthread_mutex_trylock=info
  end function forthread_mutex_trylock

  integer function forthread_mutex_unlock(mutex_id)
    integer, intent(in) :: mutex_id

    integer :: info

    call thread_mutex_unlock(mutex_id,info)
    forthread_mutex_unlock=info
  end function forthread_mutex_unlock

  integer function forthread_mutex_getprioceiling(mutex,prioceiling)
    integer, intent(in) :: mutex
    integer, intent(out) :: prioceiling

    integer :: info

    call thread_mutex_getprioceiling(mutex,prioceiling,info)
    forthread_mutex_getprioceiling=info
  end function forthread_mutex_getprioceiling

  integer function forthread_mutex_setprioceiling(mutex,prioceiling,old_ceiling)
    integer, intent(in) :: mutex, prioceiling
    integer, intent(out) :: old_ceiling

    integer :: info

    call thread_mutex_setprioceiling(mutex,prioceiling,old_ceiling,info)
    forthread_mutex_setprioceiling=info
  end function forthread_mutex_setprioceiling

#ifndef __DARWIN
  integer function forthread_mutex_timedlock(mutex,abs_timeout)
    integer,        intent(in) :: mutex
    type(timespec), intent(in) :: abs_timeout

    integer :: info

    call thread_mutex_timedlock(mutex,abs_timeout,info)
    forthread_mutex_timedlock=info
  end function forthread_mutex_timedlock
#endif

  !*****************************************!
  !*    condition variable routines        *!
  !*****************************************!

  integer function forthread_cond_destroy(cond_id)
    integer, intent(in) :: cond_id

    integer :: info

    call thread_cond_destroy(cond_id,info)
    forthread_cond_destroy=info
  end function forthread_cond_destroy

  integer function forthread_cond_init(cond_id,attr_id)
    integer, intent(out) :: cond_id
    integer, intent(in) :: attr_id

    integer :: info

    call thread_cond_init(cond_id,attr_id,info)
    forthread_cond_init=info
  end function forthread_cond_init

  integer function forthread_cond_timedwait(mutex,abstime)
    integer, intent(in) :: mutex
    type(timespec), intent(in) :: abstime

    integer :: info

    call thread_cond_timedwait(mutex,abstime,info)
    forthread_cond_timedwait=info
  end function forthread_cond_timedwait

  integer function forthread_cond_wait(cond_id,mutex_id)
    integer, intent(in) :: cond_id, mutex_id

    integer :: info

    call thread_cond_wait(cond_id,mutex_id,info)
    forthread_cond_wait=info
  end function forthread_cond_wait

  integer function forthread_cond_broadcast(cond_id)
    integer, intent(in) :: cond_id

    integer :: info

    call thread_cond_broadcast(cond_id,info)
    forthread_cond_broadcast=info
  end function forthread_cond_broadcast

  integer function forthread_cond_signal(cond_id)
    integer, intent(in) :: cond_id

    integer :: info

    call thread_cond_signal(cond_id,info)
    forthread_cond_signal=info
  end function forthread_cond_signal

#ifndef __DARWIN
  !****************************************!
  !*    barrier variable routines         *!
  !****************************************!

  integer function forthread_barrier_destroy(barrier_id)
    integer, intent(in) :: barrier_id

    integer :: info

    call thread_barrier_destroy(barrier_id,info)
    forthread_barrier_destroy=info
  end function forthread_barrier_destroy

  integer function forthread_barrier_init(barrier_id,attr_id,tcount)
    integer, intent(out) :: barrier_id
    integer, intent(in) :: attr_id, tcount

    integer :: info

    call thread_barrier_init(barrier_id,attr_id,tcount,info)
    forthread_barrier_init=info
  end function forthread_barrier_init

  integer function forthread_barrier_wait(barrier_id)
    integer, intent(in) :: barrier_id

    integer :: info

    call thread_barrier_wait(barrier_id,info)
    forthread_barrier_wait=info
  end function forthread_barrier_wait


  !*************************************!
  !*    spin variable routines         *!
  !*************************************!


  integer function forthread_spin_destroy(spinlock_id)
    integer, intent(in) :: spinlock_id

    integer :: info

    call thread_spin_destroy(spinlock_id,info)
    forthread_spin_destroy=info
  end function forthread_spin_destroy

  integer function forthread_spin_init(spinlock_id,pshared)
    integer, intent(out) :: spinlock_id
    integer, intent(in) :: pshared

    integer :: info

    call thread_spin_init(spinlock_id,pshared,info)
    forthread_spin_init=info
  end function forthread_spin_init

  integer function forthread_spin_lock(lock_id)
    integer, intent(in) :: lock_id

    integer :: info

    call thread_spin_lock(lock_id,info)
    forthread_spin_lock=info
  end function forthread_spin_lock

  integer function forthread_spin_trylock(lock_id)
    integer, intent(in) :: lock_id

    integer :: info

    call thread_spin_trylock(lock_id,info)
    forthread_spin_trylock=info
  end function forthread_spin_trylock

  integer function forthread_spin_unlock(lock_id)
    integer, intent(in) :: lock_id

    integer :: info

    call thread_spin_unlock(lock_id,info)
    forthread_spin_unlock=info
  end function forthread_spin_unlock

#endif

  !*************************************!
  !*    rwlock variable routines       *!
  !*************************************!


  integer function forthread_rwlock_destroy(rwlock_id)
    integer, intent(in) :: rwlock_id

    integer :: info

    call thread_rwlock_destroy(rwlock_id,info)
    forthread_rwlock_destroy=info
  end function forthread_rwlock_destroy

  integer function forthread_rwlock_init(rwlock_id,attr_id)
    integer, intent(out) :: rwlock_id
    integer, intent(in) :: attr_id

    integer :: info

    call thread_rwlock_init(rwlock_id,attr_id,info)
    forthread_rwlock_init=info
  end function forthread_rwlock_init

  integer function forthread_rwlock_rdlock(lock_id)
    integer, intent(in) :: lock_id

    integer :: info

    call thread_rwlock_rdlock(lock_id,info)
    forthread_rwlock_rdlock=info
  end function forthread_rwlock_rdlock

  integer function forthread_rwlock_tryrdlock(lock_id)
    integer, intent(in) :: lock_id

    integer :: info

    call thread_rwlock_tryrdlock(lock_id,info)
    forthread_rwlock_tryrdlock=info
  end function forthread_rwlock_tryrdlock

  integer function forthread_rwlock_wrlock(lock_id)
    integer, intent(in) :: lock_id

    integer :: info

    call thread_rwlock_wrlock(lock_id,info)
    forthread_rwlock_wrlock=info
  end function forthread_rwlock_wrlock

  integer function forthread_rwlock_trywrlock(lock_id)
    integer, intent(in) :: lock_id

    integer :: info

    call thread_rwlock_trywrlock(lock_id,info)
    forthread_rwlock_trywrlock=info
  end function forthread_rwlock_trywrlock

  integer function forthread_rwlock_unlock(lock_id)
    integer, intent(in) :: lock_id

    integer  :: info

    call thread_rwlock_unlock(lock_id,info)
    forthread_rwlock_unlock=info
  end function forthread_rwlock_unlock

#ifndef __DARWIN
  integer function forthread_rwlock_timedrdlock(lock_id,abs_timeout)
    integer, intent(in) :: lock_id
    type(timespec), intent(in) :: abs_timeout

    integer :: info

    call thread_rwlock_timedrdlock(lock_id,abs_timeout,info)
    forthread_rwlock_timedrdlock=info
  end function forthread_rwlock_timedrdlock

  integer function forthread_rwlock_timedwrlock(lock_id,abs_timeout)
    integer, intent(in) :: lock_id
    type(timespec), intent(in) :: abs_timeout

    integer :: info

    call thread_rwlock_timedwrlock(lock_id,abs_timeout,info)
    forthread_rwlock_timedwrlock=info
  end function forthread_rwlock_timedwrlock
#endif

  !*****************************************!
  !*      attribute object routines        *!
  !*****************************************!

  integer function forthread_attr_destroy(attr)
    integer, intent(in) :: attr

    integer :: info

    call thread_attr_destroy(attr,info)
    forthread_attr_destroy=info
  end function forthread_attr_destroy

  integer function forthread_attr_init(attr)
    integer, intent(in) :: attr

    integer :: info

    call thread_attr_init(attr,info)
    forthread_attr_init=info
  end function forthread_attr_init

  integer function forthread_attr_getdetachstate(attr,detachstate)
    integer, intent(in) :: attr
    integer, intent(out) :: detachstate

    integer  :: info

    call thread_attr_getdetachstate(attr,detachstate,info)
    forthread_attr_getdetachstate=info
  end function forthread_attr_getdetachstate

  integer function forthread_attr_setdetachstate(attr,detachstate)
    integer, intent(in)  :: attr, detachstate
    integer :: info

    call thread_attr_setdetachstate(attr,detachstate,info)
    forthread_attr_setdetachstate=info
  end function forthread_attr_setdetachstate

  integer function forthread_attr_getguardsize(attr,guardsize)
    integer, intent(in)  :: attr
    integer(size_t), intent(out) :: guardsize

    integer :: info

    call thread_attr_getguardsize(attr,guardsize,info)
    forthread_attr_getguardsize=info
  end function forthread_attr_getguardsize

  integer function forthread_attr_setguardsize(attr,guardsize)
    integer, intent(in) :: attr
    integer(size_t), intent(in) :: guardsize

    integer :: info

    call thread_attr_setguardsize(attr,guardsize,info)
    forthread_attr_setguardsize=info
  end function forthread_attr_setguardsize

  integer function forthread_attr_getinheritsched(attr,inheritsched)
    integer, intent(in) :: attr
    integer, intent(out) :: inheritsched

    integer :: info

    call thread_attr_getinheritsched(attr,inheritsched,info)
    forthread_attr_getinheritsched=info
  end function forthread_attr_getinheritsched

  integer function forthread_attr_setinheritsched(attr,inheritsched)
    integer, intent(in) :: attr
    integer, intent(in) :: inheritsched

    integer :: info

    call thread_attr_setinheritsched(attr,inheritsched,info)
    forthread_attr_setinheritsched=info
  end function forthread_attr_setinheritsched

  integer function forthread_attr_getschedparam(attr,param)
    integer, intent(in) :: attr
    type(sched_param), intent(out) :: param

    integer :: info

    call thread_attr_getschedparam(attr,param,info)
    forthread_attr_getschedparam=info
  end function forthread_attr_getschedparam

  integer function forthread_attr_setschedparam(attr,param)
    integer, intent(in) :: attr
    type(sched_param), intent(in) :: param

    integer :: info

    call thread_attr_setschedparam(attr,param,info)
    forthread_attr_setschedparam=info
  end function forthread_attr_setschedparam

  integer function forthread_attr_getschedpolicy(attr,policy)
    integer, intent(in) :: attr
    integer, intent(out) :: policy

    integer :: info

    call thread_attr_getschedpolicy(attr,policy,info)
    forthread_attr_getschedpolicy=info
  end function forthread_attr_getschedpolicy

  integer function forthread_attr_setschedpolicy(attr,policy)
    integer, intent(in) :: attr, policy

    integer :: info

    call thread_attr_setschedpolicy(attr,policy,info)
    forthread_attr_setschedpolicy=info
  end function forthread_attr_setschedpolicy

  integer function forthread_attr_getscope(attr,scope)
    integer, intent(in) :: attr
    integer, intent(out) :: scope

    integer :: info

    call thread_attr_getscope(attr,scope,info)
    forthread_attr_getscope=info
  end function forthread_attr_getscope

  integer function forthread_attr_setscope(attr,scope)
    integer, intent(in) :: attr, scope

    integer :: info

    call thread_attr_setscope(attr,scope,info)
    forthread_attr_setscope=info
  end function forthread_attr_setscope

  integer function forthread_attr_getstacksize(attr,stacksize)
    integer, intent(in) :: attr
    integer(size_t), intent(out) :: stacksize

    integer :: info

    call thread_attr_getstacksize(attr,stacksize,info)
    forthread_attr_getstacksize=info
  end function forthread_attr_getstacksize

  integer function forthread_attr_setstacksize(attr,stacksize)
    integer, intent(in) :: attr
    integer(size_t), intent(in) :: stacksize

    integer :: info

    call thread_attr_setstacksize(attr,stacksize,info)
    forthread_attr_setstacksize=info
  end function forthread_attr_setstacksize

  !*****************************************!
  !*       mutex attribute routines        *!
  !*****************************************!

  integer function forthread_mutexattr_destroy(attr)
    integer, intent(in) :: attr

    integer :: info

    call thread_mutexattr_destroy(attr,info)
    forthread_mutexattr_destroy=info
  end function forthread_mutexattr_destroy

  integer function forthread_mutexattr_init(attr)
    integer, intent(in) :: attr

    integer :: info

    call thread_mutexattr_init(attr,info)
    forthread_mutexattr_init=info
  end function forthread_mutexattr_init

  integer function forthread_mutexattr_getpshared(attr,pshared)
    integer, intent(in) :: attr
    integer, intent(out) :: pshared

    integer :: info

    call thread_mutexattr_getpshared(attr,pshared,info)
    forthread_mutexattr_getpshared=info
  end function forthread_mutexattr_getpshared

  integer function forthread_mutexattr_setpshared(attr,pshared)
    integer       , intent(in)      :: attr
    integer       , intent(in)      :: pshared
    integer           :: info

    call thread_mutexattr_setpshared(attr,pshared,info)
    forthread_mutexattr_setpshared=info
  end function forthread_mutexattr_setpshared

  integer function forthread_mutexattr_getprioceiling(attr,prioceiling)
    integer, intent(in) :: attr
    integer, intent(out) :: prioceiling

    integer :: info

    call thread_mutexattr_getprioceiling(attr,prioceiling,info)
    forthread_mutexattr_getprioceiling=info
  end function forthread_mutexattr_getprioceiling

  integer function forthread_mutexattr_setprioceiling(attr,prioceiling)
    integer, intent(in) :: attr, prioceiling

    integer :: info

    call thread_mutexattr_setprioceiling(attr,prioceiling,info)
    forthread_mutexattr_setprioceiling=info
  end function forthread_mutexattr_setprioceiling

  integer function forthread_mutexattr_getprotocol(attr,protocol)
    integer, intent(in) :: attr
    integer, intent(out) :: protocol

    integer :: info

    call thread_mutexattr_getprotocol(attr,protocol,info)
    forthread_mutexattr_getprotocol=info
  end function forthread_mutexattr_getprotocol

  integer function forthread_mutexattr_setprotocol(attr,protocol)
    integer, intent(in) :: attr, protocol

    integer :: info

    call thread_mutexattr_setprotocol(attr,protocol,info)
    forthread_mutexattr_setprotocol=info
  end function forthread_mutexattr_setprotocol

  integer function forthread_mutexattr_gettype(attr,mtype)
    integer, intent(in) :: attr
    integer, intent(out) :: mtype

    integer :: info

    call thread_mutexattr_gettype(attr,mtype,info)
    forthread_mutexattr_gettype=info
  end function forthread_mutexattr_gettype

  integer function forthread_mutexattr_settype(attr,mtype)
    integer, intent(in) :: attr, mtype

    integer :: info

    call thread_mutexattr_settype(attr,mtype,info)
    forthread_mutexattr_settype=info
  end function forthread_mutexattr_settype

  !*****************************************************!
  !*    condition attriubute variable routines         *!
  !*****************************************************!

  integer function forthread_condattr_destroy(attr)
    integer, intent(in) :: attr

    integer :: info

    call thread_condattr_destroy(attr,info)
    forthread_condattr_destroy=info
  end function forthread_condattr_destroy

  integer function forthread_condattr_init(attr)
    integer, intent(in) :: attr

    integer :: info

    call thread_condattr_init(attr,info)
    forthread_condattr_init=info
  end function forthread_condattr_init

  integer function forthread_condattr_getpshared(attr,pshared)
    integer, intent(in) :: attr
    integer, intent(out) :: pshared

    integer :: info

    call thread_condattr_getpshared(attr,pshared,info)
    forthread_condattr_getpshared=info
  end function forthread_condattr_getpshared

  integer function forthread_condattr_setpshared(attr,pshared)
    integer, intent(in) :: attr, pshared

    integer :: info

    call thread_condattr_setpshared(attr,pshared,info)
    forthread_condattr_setpshared=info
  end function forthread_condattr_setpshared

#ifndef __DARWIN
  integer function forthread_condattr_getclock(attr,clock_id)
    integer, intent(in) :: attr
    integer, intent(out) :: clock_id

    integer :: info

    call thread_condattr_getclock(attr,clock_id,info)
    forthread_condattr_getclock=info
  end function forthread_condattr_getclock

  integer function forthread_condattr_setclock(attr,clock_id)
    integer, intent(in) :: attr
    integer, intent(in) :: clock_id

    integer :: info

    call thread_condattr_setclock(attr,clock_id,info)
    forthread_condattr_setclock=info
  end function forthread_condattr_setclock

  !**************************************************!
  !*    barrier attribute variable routines         *!
  !**************************************************!

  integer function forthread_barrierattr_destroy(attr)
    integer, intent(in) :: attr

    integer :: info

    call thread_barrierattr_destroy(attr,info)
    forthread_barrierattr_destroy=info
  end function forthread_barrierattr_destroy

  integer function forthread_barrierattr_init(attr)
    integer, intent(in) :: attr

    integer :: info

    call thread_barrierattr_init(attr,info)
    forthread_barrierattr_init=info
  end function forthread_barrierattr_init

  integer function forthread_barrierattr_getpshared(attr,pshared)
    integer, intent(in) :: attr
    integer, intent(out) :: pshared

    integer :: info

    call thread_barrierattr_getpshared(attr,pshared,info)
    forthread_barrierattr_getpshared=info
  end function forthread_barrierattr_getpshared

  integer function forthread_barrierattr_setpshared(attr,pshared)
    integer, intent(in) :: attr, pshared
    integer :: info

    call thread_barrierattr_setpshared(attr,pshared,info)
    forthread_barrierattr_setpshared=info
  end function forthread_barrierattr_setpshared
#endif

  !**************************************************!
  !*    rwlock attribute variable routines         *!
  !**************************************************!

  integer function forthread_rwlockattr_destroy(attr)
    integer, intent(in) :: attr
    integer :: info

    call thread_rwlockattr_destroy(attr,info)
    forthread_rwlockattr_destroy=info
  end function forthread_rwlockattr_destroy

  integer function forthread_rwlockattr_init(attr)
    integer, intent(in) :: attr
    integer :: info

    call thread_rwlockattr_init(attr,info)
    forthread_rwlockattr_init=info
  end function forthread_rwlockattr_init

  integer function forthread_rwlockattr_getpshared(attr,pshared)
    integer, intent(in) :: attr
    integer, intent(out) :: pshared
    integer :: info

    call thread_rwlockattr_getpshared(attr,pshared,info)
    forthread_rwlockattr_getpshared=info
  end function forthread_rwlockattr_getpshared

  integer function forthread_rwlockattr_setpshared(attr,pshared)
    integer, intent(in) :: attr, pshared
    integer :: info

    call thread_rwlockattr_setpshared(attr,pshared,info)
    forthread_rwlockattr_setpshared=info
  end function forthread_rwlockattr_setpshared
end module forthread_mod

