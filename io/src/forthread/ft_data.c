#include "ft_data.h"

int is_initialized;
array_t *threads;
array_t *thread_attrs;
array_t *thread_keys;
array_t *mutexes;
array_t *mutex_attrs;
array_t *once_ctrls;
array_t *conds;
array_t *cond_attrs ;
array_t *barriers;
array_t *barrier_attrs;
varray_t *spinlocks;
array_t *rwlocks;
array_t *rwlock_attrs;

void thread_init_internal(int *info) {
  int i = 0;
  pthread_t stid;
  static int init = 0;
  *info = FT_OK;

  if (init) {
    *info = FT_EINIT;
    return;
  }
  threads = NULL;
  array_init(&threads,INIT_SIZE);
  thread_attrs = NULL;
  array_init(&thread_attrs,INIT_SIZE);
  thread_keys = NULL;
  array_init(&thread_keys,INIT_SIZE);
  once_ctrls = NULL;
  array_init(&once_ctrls,INIT_SIZE);
  mutexes = NULL;
  array_init(&mutexes,INIT_SIZE);
  mutex_attrs = NULL;
  array_init(&mutex_attrs,INIT_SIZE);
  conds = NULL;
  array_init(&conds,INIT_SIZE);
  cond_attrs = NULL;
  array_init(&cond_attrs,INIT_SIZE);
  barriers = NULL;
  array_init(&barriers,INIT_SIZE);
  barrier_attrs = NULL;
  array_init(&barrier_attrs,INIT_SIZE);
  spinlocks = NULL;
  varray_init(&spinlocks,INIT_SIZE);
  rwlocks = NULL;
  array_init(&rwlocks,INIT_SIZE);
  rwlock_attrs = NULL;
  array_init(&rwlock_attrs,INIT_SIZE);
  // allocate and store the thread master ID
  threads->data[0] = (pthread_t*) malloc(sizeof(pthread_t));
  stid = pthread_self();
  memcpy(threads->data[0],&stid,sizeof(pthread_t));
  threads->after++;
  
  init = 1;
  is_initialized = init;
}

void thread_destroy_internal(int* info) {
  int id;
  for(id = 1; id < threads->after; id++) {
    thread_cancel(&id,info);
  }
  array_delete(threads);
  array_delete(thread_attrs);
  array_delete(thread_keys);
  array_delete(once_ctrls);
  for(id = 0; id < mutexes->after; id++) {
    thread_mutex_destroy(&id,info);
  }
  array_delete(mutexes);
  array_delete(mutex_attrs);
  for(id = 0; id < conds->after; id++) {
    thread_cond_destroy(&id,info);
  }
  array_delete(conds);
  array_delete(cond_attrs);

#ifdef _POSIX_BARRIERS
  for(id = 0; id < barriers->after; id++) {
    thread_barrier_destroy(&id,info);
  }
  array_delete(barriers);
  array_delete(barrier_attrs);
#endif
#ifndef __DARWIN
  for(id = 0; id < spinlocks->after; id++) {
    thread_spin_destroy(&id,info);
  }
  varray_delete(spinlocks);
#endif
  for(id = 0; id < rwlocks->after; id++) {
    thread_rwlock_destroy(&id,info);
  }
  array_delete(rwlocks);
  array_delete(rwlock_attrs);
  *info = FT_OK;
}

/**
 * Initializes a given array. The argument array must be either
 * already allocated or a NULL pointer.
 **/
void array_init(array_t **array,int size) {
  int i;

  if (*array == NULL)
    *array = (array_t*)malloc(sizeof(array_t));
  pthread_mutex_init(&((*array)->mutex),NULL);
  (*array)->data = (void**)malloc(sizeof(void*)*size);
  for(i = 0; i < size; i++)
    (*array)->data[i] = NULL;
  (*array)->size = size;
  (*array)->after = 0;
}

/**
 * Initializes a given varray. The argument array must be either
 * already allocated or a NULL pointer.
 **/
void varray_init(varray_t **array,int size) {
  int i;

  if (*array == NULL)
    *array = (varray_t*)malloc(sizeof(varray_t));
  pthread_mutex_init(&((*array)->mutex),NULL);
  (*array)->data = (volatile void**)malloc(sizeof(void*)*size);
  for(i = 0; i < size; i++)
    (*array)->data[i] = NULL;
  (*array)->size = size;
  (*array)->after = 0;
}

/**
 * Resize array to size. We assume array to be NOT NULL
 **/
void array_resize(array_t **array,int size) {
  int i;

  (*array)->data = (void**)realloc((*array)->data,sizeof(void*)*size);
  (*array)->size = size;

  for(i = (*array)->after; i < size; i++)
    (*array)->data[i] = NULL;

}

/**
 * Resize varray to size. We assume varray to be NOT NULL
 **/
void varray_resize(varray_t **array,int size) {
  int i;

  (*array)->data = (volatile void**)realloc((*array)->data,sizeof(volatile void*)*size);
  (*array)->size = size;

  for(i = (*array)->after; i < size; i++)
    (*array)->data[i] = NULL;

}

/**
 * Free memory for array
 **/
void array_delete(array_t *array) {
  free(array->data);
  free(array);
}

/**
 * Free memory for varray
 **/
void varray_delete(varray_t *array) {
  free(array->data);
  free(array);
}

/**
 * A simple helper function to check whether an ID
 * is pointing to a valid element in arr. This assumes
 * that changes in arr are alwys done using this library.
 *
 * This function is not thread safe
 **/
int is_valid(array_t *arr, int id) {
  if ((id >= 0) && (id < arr->after) &&
      (arr->data[id] != NULL))
    return 1;
  else
    return 0;
}

/**
 * (varray version)
 * A simple helper function to check whether an ID
 * is pointing to a valid element in arr. This assumes
 * that changes in arr are alwys done using this library.
 *
 * This function is not thread safe
 **/
int vis_valid(varray_t *arr, int id) {
  if ((id >= 0) && (id < arr->after) && 
      (arr->data[id] != NULL))
    return 1;
  else
    return 0;
}
