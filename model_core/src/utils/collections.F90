!> Collection data structures.
!!
!! The collections utilities which provide common collection structures for different components of
!! MONC. Currently a list, map, stack and queue all with appropriate functionality are provided.
!! The core of all collections is currently a doubly linked list, this is abstracted to allow for the
!! internal structure of collections to change without requiring any user code modifications.
module collections_mod
  use datadefn_mod, only : STRING_LENGTH, DEFAULT_PRECISION
  use conversions_mod, only : conv_to_generic, generic_to_double_real, conv_to_integer, conv_to_string, conv_to_logical, &
       conv_single_real_to_double, conv_to_real
  use logging_mod, only : LOG_ERROR, log_log
  implicit none

#ifndef TEST_MODE
  private
#endif

  !> Number of entries in the hash table, this is a tradeoff - larger means more memory but smaller runtime
  !! assuming a hashing function with good distribution.
  integer, parameter, private :: hash_size  = 4993

  !> \private Private list node which holds the raw generic node data and pointers to next and previous list nodes
  type listnode_type
     type(listnode_type), pointer :: next=>null(),& !< The next nodes in the doubly linked list
          prev=>null() !< The next nodes in the doubly linked list
     !> Pointer to the generic data which is held by this node
     class(*), pointer :: data => null()
     logical :: memory_allocation_automatic
  end type listnode_type

  !> \private Private map key-value pair data structure
  type mapnode_type
     logical :: memory_allocation_automatic
     !> String key
     character(len=STRING_LENGTH) :: key
     !> Pointer to the generic data which is held by this key-value pair
     class(*), pointer :: value => null()
  end type mapnode_type

  !> \private Private set key structure
  type setnode_type
     !> String key
     character(len=STRING_LENGTH) :: key
  end type setnode_type

  type, public :: mapentry_type
     character(len=STRING_LENGTH) :: key
     class(*), pointer :: value => null()
  end type mapentry_type  

  type, public :: iterator_type
     type(listnode_type), pointer :: next_item
     type(list_type), dimension(:), pointer :: hash_structure
     integer :: hash_ptr
  end type iterator_type

  !> List data structure which implements a doubly linked list. This list will preserve its order
  !!
  !! Note that it holds a pointer (reference) to the data
  type, public :: list_type
     type(listnode_type), pointer, private :: head=>null(),& !< The head of the doubly linked list
          tail=>null()  !< The tail of the doubly linked list
     !> Number of elements in the list
     integer, private :: size=0
  end type list_type

  !> Queue (FIFO) data structure
  !!
  !! Note that it holds a pointer (reference) to the data
  type, public :: queue_type
     !> Underlying list data structure used to implement the queue
     type(list_type), private :: queue_ds
  end type queue_type

  !> Stack (FILO) data structure
  !!
  !! Note that it holds a pointer (reference) to the data
  type, public :: stack_type
     !> Underlying list data structure used to implement the stack
     type(list_type), private :: stack_ds
  end type stack_type

  !> Map data structure that holds string (length 20 maximum) key value pairs
  !!
  !! Note that it holds a pointer (reference) to the value data
  type, public :: map_type
     !> Underlying list data structure used to implement the map
     type(list_type), private :: map_ds
  end type map_type

  !> A hashmap structure, the same as a map but uses hashing for greatly improved performance when storing large numbers of
  !! entries, with the lookup complexisty being almost constant (assuming good hash distribution.) Note that this does 
  !! require more storage than a map and unlike a map does not preserve ordering
  type, public :: hashmap_type
     !> Underlying list data structure used to implement the map
     type(list_type), pointer, dimension(:), private :: map_ds => null()
     integer, private :: size=0
  end type hashmap_type

  !> Hashset structure which will store unique strings. The hashing aspect means that lookup is very fast but it does
  !! add extra memory overhead and the order is non-deterministic
  type, public :: hashset_type
     type(list_type), pointer, dimension(:), private :: set_ds => null()
     integer, private :: size=0
  end type hashset_type  

  !> Pushes a generic element onto the stack or queue
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific stack or queue involved
  !! @param data Pointer to the generic data to push onto the collection
  !! @param memory_allocation_automatic Whether the collections API should manage the freeing of memory
  interface c_push_generic
     module procedure stack_push_generic, queue_push_generic
  end interface c_push_generic

  !> Pushes an integer element onto the stack or queue
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific stack or queue involved
  !! @param data Integer data to push onto the collection
  interface c_push_integer
     module procedure stack_push_int, queue_push_int
  end interface c_push_integer

  !> Pushes a string element onto the stack or queue
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific stack or queue involved
  !! @param data String data to push onto the collection
  interface c_push_string
     module procedure stack_push_string, queue_push_string
  end interface c_push_string

  !> Pushes a double precision real element onto the stack or queue
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific stack or queue involved
  !! @param data Double precision real data to push onto the collection
  interface c_push_real
     module procedure stack_push_real, queue_push_real
  end interface c_push_real

  !> Pushes a logical element onto the stack or queue
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific stack or queue involved
  !! @param data Logical data to push onto the collection
  interface c_push_logical
     module procedure stack_push_logical, queue_push_logical
  end interface c_push_logical

  !> Pops a generic element off the stack or queue
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific stack or queue involved
  !! @returns Pointer to the generic data element which has been popped off the collection
  interface c_pop_generic
     module procedure stack_pop_generic, queue_pop_generic
  end interface c_pop_generic

  !> Pops an integer element off the stack or queue
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific stack or queue involved
  !! @returns Integer element which has been popped off the collection or raises and error if none is found
  interface c_pop_integer
     module procedure stack_pop_int, queue_pop_int
  end interface c_pop_integer

  !> Pops a string off the stack or queue
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific stack or queue involved
  !! @returns String which has been popped off the collection or raises and error if none is found
  interface c_pop_string
     module procedure stack_pop_string, queue_pop_string
  end interface c_pop_string

  !> Pops a double precision real element off the stack or queue
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific stack or queue involved
  !! @returns Double precision real element which has been popped off the collection or raises and error if none is found
  interface c_pop_real
     module procedure stack_pop_real, queue_pop_real
  end interface c_pop_real

  !> Pops a logical element off the stack or queue
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific stack or queue involved
  !! @returns Logical element which has been popped off the collection or raises and error if none is found
  interface c_pop_logical
     module procedure stack_pop_logical, queue_pop_logical
  end interface c_pop_logical

  !> Adds a generic element to the end of the list
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific list involved
  !! @param data Pointer to the raw generic data to add
  !! @param memory_allocation_automatic Whether the collections API should manage the freeing of memory
  interface c_add_generic
     module procedure list_add_generic
  end interface c_add_generic

  !> Adds an integer element to the end of the list
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific list involved
  !! @param data Integer data to add
  interface c_add_integer
     module procedure list_add_int
  end interface c_add_integer

  !> Adds a string to the end of the list
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific list involved
  !! @param data String data to add
  interface c_add_string
     module procedure list_add_string, hashset_add
  end interface c_add_string

  !> Adds a double precision real element to the end of the list
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific list involved
  !! @param data Double precision real data to add
  interface c_add_real
     module procedure list_add_real
  end interface c_add_real

  !> Adds a logical element to the end of the list
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific list involved
  !! @param data Logical data to add
  interface c_add_logical
     module procedure list_add_logical
  end interface c_add_logical

  !> Inserts a generic element into the list or places at the end if the index > list size
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific list involved
  !! @param data Pointer to the generic data to insert
  !! @param index The list index to insert to
  !! @param memory_allocation_automatic Whether the collections API should manage the freeing of memory
  interface c_insert_generic
     module procedure list_insert_generic
  end interface c_insert_generic

  !> Inserts an integer element into the list or places at the end if the index > list size
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific list involved
  !! @param data Integer data to insert
  !! @param index The list index to insert to
  interface c_insert_integer
     module procedure list_insert_int
  end interface c_insert_integer

  !> Inserts a string into the list or places at the end if the index > list size
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific list involved
  !! @param data String data to insert
  !! @param index The list index to insert to
  interface c_insert_string
     module procedure list_insert_string
  end interface c_insert_string

  !> Inserts a double precision real element into the list or places at the end if the index > list size
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific list involved
  !! @param data Double precision real data to insert
  !! @param index The list index to insert to
  interface c_insert_real
     module procedure list_insert_real
  end interface c_insert_real

  !> Inserts a logical element into the list or places at the end if the index > list size
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific list involved
  !! @param data Logical data to insert
  !! @param index The list index to insert to
  interface c_insert_logical
     module procedure list_insert_logical
  end interface c_insert_logical

  !> Puts a generic key-value pair into the map.
  !!
  !! If the key is not already held in the map then
  !! the key-value pair will be added, otherwise the existing key-value pair will be modified to
  !! hold this updated value (keys must be unique)
  !! This has a time complexity of O(n) due to key look up
  !! @param collection The specific map involved
  !! @param key The key to place in the map
  !! @param value Pointer to the generic data value to place in the map
  !! @param memory_allocation_automatic Whether the collections API should manage the freeing of memory
  interface c_put_generic
     module procedure map_put_generic, hashmap_put_generic
  end interface c_put_generic

  !> Puts an integer key-value pair into the map.
  !!
  !! If the key is not already held in the map then
  !! the key-value pair will be added, otherwise the existing key-value pair will be modified to
  !! hold this updated value (keys must be unique)
  !! This has a time complexity of O(n) due to key look up
  !! @param collection The specific map involved
  !! @param key The key to place in the map
  !! @param value Integer data value to place in the map
  interface c_put_integer
     module procedure map_put_int, hashmap_put_int
  end interface c_put_integer

  !> Puts a string key-value pair into the map.
  !!
  !! If the key is not already held in the map then
  !! the key-value pair will be added, otherwise the existing key-value pair will be modified to
  !! hold this updated value (keys must be unique)
  !! This has a time complexity of O(n) due to key look up
  !! @param collection The specific map involved
  !! @param key The key to place in the map
  !! @param value String data value to place in the map
  interface c_put_string
     module procedure map_put_string, hashmap_put_string
  end interface c_put_string

  !> Puts a double precision real key-value pair into the map.
  !!
  !! If the key is not already held in the map then
  !! the key-value pair will be added, otherwise the existing key-value pair will be modified to
  !! hold this updated value (keys must be unique)
  !! This has a time complexity of O(n) due to key look up
  !! @param collection The specific map involved
  !! @param key The key to place in the map
  !! @param value Double precision real data value to place in the map
  interface c_put_real
     module procedure map_put_real, hashmap_put_real
  end interface c_put_real

  !> Puts a logical key-value pair into the map.
  !!
  !! If the key is not already held in the map then
  !! the key-value pair will be added, otherwise the existing key-value pair will be modified to
  !! hold this updated value (keys must be unique)
  !! This has a time complexity of O(n) due to key look up
  !! @param collection The specific map involved
  !! @param key The key to place in the map
  !! @param value Logical data value to place in the map
  interface c_put_logical
     module procedure map_put_logical, hashmap_put_logical
  end interface c_put_logical

  !> Gets a specific generic element out of the list, stack, queue or map with the corresponding key
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific list, stack, queue or map involved
  !! @param key String look up key
  !! @returns Generic pointer to the value associated with the key or null if none exists
  interface c_get_generic
     module procedure list_get_generic, stack_get_generic, queue_get_generic, map_get_generic, hashmap_get_generic, &
          mapentry_get_generic
  end interface c_get_generic

  !> Gets a specific integer element out of the list, stack, queue or map with the corresponding key
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific list, stack, queue or map involved
  !! @param key String look up key
  !! @returns Integer value associated with the key or raises an error if none exists
  interface c_get_integer
     module procedure list_get_int, stack_get_int, queue_get_int, map_get_int, hashmap_get_int, mapentry_get_int
  end interface c_get_integer

  !> Gets a specific string element out of the list, stack, queue or map with the corresponding key
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific list, stack, queue or map involved
  !! @param key String look up key
  !! @returns String value associated with the key or raises an error if none exists
  interface c_get_string
     module procedure list_get_string, stack_get_string, queue_get_string, map_get_string, hashmap_get_string, &
          hashset_get_string, mapentry_get_string
  end interface c_get_string

  !> Gets a specific double precision real element out of the list, stack, queue or map with the corresponding key
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific list, stack, queue or map involved
  !! @param key String look up key
  !! @returns Double precision real value associated with the key or raises an error if none exists
  interface c_get_real
     module procedure list_get_real, stack_get_real, queue_get_real, map_get_real, hashmap_get_real, mapentry_get_real
  end interface c_get_real

  !> Gets a specific logical element out of the list, stack, queue or map with the corresponding key
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific list, stack, queue or map involved
  !! @param key String look up key
  !! @returns Logical value associated with the key or raises an error if none exists
  interface c_get_logical
     module procedure list_get_logical, stack_get_logical, queue_get_logical, map_get_logical, &
          hashmap_get_logical, mapentry_get_logical
  end interface c_get_logical

  !> Removes a specific element from the list or map
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific list or map involved
  !! @param identifier In the cast of a list the index to remove or for a map the key (of the key-value pair) to remove
  interface c_remove
     module procedure list_remove, map_remove, hashmap_remove, hashset_remove
  end interface c_remove

  !> Returns the number of elements in the collection
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific list, queue, stack or map involved
  !! @returns Number of elements in the collection
  interface c_size
     module procedure list_size, stack_size, queue_size, map_size, hashmap_size, hashset_size
  end interface c_size

  !> Returns whether a collection is empty
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific list, stack queue or map involved
  !! @returns Whether the collection is empty
  interface c_is_empty
     module procedure list_is_empty, stack_is_empty, queue_is_empty, map_is_empty, hashmap_is_empty, hashset_is_empty
  end interface c_is_empty

  !> Determines whether or not a map contains a specific key
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific map involved
  !! @param key The key we are looking up
  !! @returns Whether the map contains the key
  interface c_contains
     module procedure map_contains_key, hashmap_contains_key, hashset_contains
  end interface c_contains

  !> Retrieves the key currently being held at a specific index in the map or "" if the index > map elements
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific map involved
  !! @param index The index to retrieve the key from
  !! @returns The key string
  interface c_key_at
     module procedure map_key_at, hashmap_key_at
  end interface c_key_at

  !> Retrieves the generic value held at the specific map index or null if index > map elements
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific map involved
  !! @param index The index to get value from
  !! @returns Generic pointer to the value
  interface c_generic_at
     module procedure map_generic_at, hashmap_generic_at
  end interface c_generic_at

  !> Retrieves the integer value held at the specific map index or null if index > map elements
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific map involved
  !! @param index The index to get value from
  !! @returns Integer value or raises an error if none is found
  interface c_integer_at
     module procedure map_integer_at, hashmap_integer_at
  end interface c_integer_at

  !> Retrieves the string value held at the specific map index or null if index > map elements
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific map involved
  !! @param index The index to get value from
  !! @returns String value or raises an error if none is found
  interface c_string_at
     module procedure map_string_at, hashmap_string_at
  end interface c_string_at

  !> Retrieves the double precision real value held at the specific map index or null if index > map elements
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific map involved
  !! @param index The index to get value from
  !! @returns Double precision real value or raises an error if none is found
  interface c_real_at
     module procedure map_real_at, hashmap_real_at
  end interface c_real_at

  !> Retrieves the logical value held at the specific map index or null if index > map elements
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific map involved
  !! @param index The index to get value from
  !! @returns Logical value or raises an error if none is found
  interface c_logical_at
     module procedure map_logical_at, hashmap_logical_at
  end interface c_logical_at

  !> Retrieves a map entry at a specific index or null if index > map elements. This is more efficient than calling
  !! key at and then value at (or get with the key) as only requires one search for both the key and value
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific map involved
  !! @param index The index to get value from
  !! @param key The associated key
  !! @param value Generic pointer to corresponding value
  interface c_generic_entry_at
     module procedure map_generic_entry_at, hashmap_generic_entry_at
  end interface c_generic_entry_at

  !> Retrieves a map entry at a specific index. This is more efficient than calling
  !! key at and then value at (or get with the key) as only requires one search for both the key and value
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific map involved
  !! @param index The index to get value from
  !! @param key The associated key
  !! @param value Integer value or raises an error if none is found
  interface c_integer_entry_at
     module procedure map_integer_entry_at, hashmap_integer_entry_at
  end interface c_integer_entry_at

  !> Retrieves a map entry at a specific index. This is more efficient than calling
  !! key at and then value at (or get with the key) as only requires one search for both the key and value
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific map involved
  !! @param index The index to get value from
  !! @param key The associated key
  !! @param value String value or raises an error if none is found
  interface c_string_entry_at
     module procedure map_string_entry_at, hashmap_string_entry_at
  end interface c_string_entry_at

  !> Retrieves a map entry at a specific index. This is more efficient than calling
  !! key at and then value at (or get with the key) as only requires one search for both the key and value
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific map involved
  !! @param index The index to get value from
  !! @param key The associated key
  !! @param value Double precision real value or raises an error if none is found
  interface c_real_entry_at
     module procedure map_real_entry_at, hashmap_real_entry_at
  end interface c_real_entry_at

  !> Retrieves a map entry at a specific index. This is more efficient than calling
  !! key at and then value at (or get with the key) as only requires one search for both the key and value
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific map involved
  !! @param index The index to get value from
  !! @param key The associated key
  !! @param value Logical value or raises an error if none is found
  interface c_logical_entry_at
     module procedure map_logical_entry_at, hashmap_logical_entry_at
  end interface c_logical_entry_at

  !> Frees up all the allocatable, heap, memory associated with a list, stack, queue or map.
  !!
  !! This basically acts like a clear operation and once freed the data structure can be reused
  !! Note that it does not free the data memory at all as this might be referenced else where in the code
  !! This has a time complexity of O(n)
  !! @param collection The list, stack, queue or map to delete all members and free all allocated memory of
  interface c_free
     module procedure list_free, stack_free, queue_free, map_free, hashmap_free, hashset_free
  end interface c_free

  interface c_get_iterator
     module procedure list_get_iterator, map_get_iterator, hashmap_get_iterator, hashset_get_iterator, stack_get_iterator, &
          queue_get_iterator
  end interface c_get_iterator  

  interface c_has_next
     module procedure iteratior_has_next
  end interface c_has_next
  
  interface c_next_integer
     module procedure iterator_get_next_integer
  end interface c_next_integer

  interface c_next_string
     module procedure iterator_get_next_string
  end interface c_next_string

  interface c_next_real
     module procedure iterator_get_next_real
  end interface c_next_real

  interface c_next_logical
     module procedure iterator_get_next_logical
  end interface c_next_logical

  interface c_next_mapentry
     module procedure iterator_get_next_mapentry
  end interface c_next_mapentry

  interface c_next_generic
     module procedure iterator_get_next_generic
  end interface c_next_generic

  ! Explicit public interfaces and data items
  public c_push_generic, c_push_integer, c_push_string, c_push_real, c_push_logical, c_pop_generic, c_pop_integer, c_pop_real, &
       c_pop_string, c_pop_logical, c_add_generic, c_add_integer, c_add_string, c_add_real, c_add_logical, c_insert_generic, &
       c_insert_integer, c_insert_string, c_insert_real, c_insert_logical, c_put_generic, c_put_integer, c_put_string, &
       c_put_real, c_put_logical, c_get_generic, c_get_integer, c_get_string, c_get_real, c_get_logical, c_remove, c_size, &
       c_is_empty, c_contains, c_key_at, c_generic_at, c_integer_at, c_string_at, c_real_at, c_logical_at, c_generic_entry_at, &
       c_integer_entry_at, c_string_entry_at, c_real_entry_at, c_logical_entry_at, c_free, c_get_iterator, c_has_next, &
       c_next_integer, c_next_string, c_next_real, c_next_logical, c_next_mapentry, c_next_generic
contains

  !> Retrieves an iterator representation of the map, ready to access the first element
  !! @param specificmap Specific collection to base this iterator on
  !! @returns The iterator ready to access the first element
  type(iterator_type) function map_get_iterator(specificmap)
    type(map_type), intent(inout) :: specificmap

    map_get_iterator%next_item=>specificmap%map_ds%head
    map_get_iterator%hash_structure=>null()
    map_get_iterator%hash_ptr=0
  end function map_get_iterator 

  !> Puts a specific key-value pair into the map.
  !!
  !! If the key is not already held in the map then
  !! the key-value pair will be added, otherwise the existing key-value pair will be modified to
  !! hold this updated value (keys must be unique)
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key The key to place in the map
  !! @param data Integer value to place in the map
  subroutine map_put_int(specificmap, key, int_data)
    type(map_type), intent(inout) :: specificmap
    integer, intent(in) :: int_data
    character(len=*), intent(in) :: key

    class(*), pointer :: generic

    generic=>conv_to_generic(int_data, .true.)
    call map_put_generic(specificmap, key, generic, .true.)
  end subroutine map_put_int

  !> Puts a specific key-value pair into the map.
  !!
  !! If the key is not already held in the map then
  !! the key-value pair will be added, otherwise the existing key-value pair will be modified to
  !! hold this updated value (keys must be unique)
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key The key to place in the map
  !! @param data String value to place in the map
  subroutine map_put_string(specificmap, key, str_data)
    type(map_type), intent(inout) :: specificmap
    character(len=STRING_LENGTH), intent(in) :: str_data
    character(len=*), intent(in) :: key

    class(*), pointer :: generic

    generic=>conv_to_generic(str_data, .true.)
    call map_put_generic(specificmap, key, generic, .true.)
  end subroutine map_put_string
  
  !> Puts a specific key-value pair into the map.
  !!
  !! If the key is not already held in the map then
  !! the key-value pair will be added, otherwise the existing key-value pair will be modified to
  !! hold this updated value (keys must be unique)
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key The key to place in the map
  !! @param data Double precision real value to place in the map
  subroutine map_put_real(specificmap, key, real_data)
    type(map_type), intent(inout) :: specificmap
    real(kind=DEFAULT_PRECISION), intent(in) :: real_data
    character(len=*), intent(in) :: key

    class(*), pointer :: generic

    generic=>conv_to_generic(real_data, .true.)
    call map_put_generic(specificmap, key, generic, .true.)
  end subroutine map_put_real

  !> Puts a specific key-value pair into the map.
  !!
  !! If the key is not already held in the map then
  !! the key-value pair will be added, otherwise the existing key-value pair will be modified to
  !! hold this updated value (keys must be unique)
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key The key to place in the map
  !! @param data Logical value to place in the map
  subroutine map_put_logical(specificmap, key, logical_data)
    type(map_type), intent(inout) :: specificmap
    logical, intent(in) :: logical_data
    character(len=*), intent(in) :: key

    class(*), pointer :: generic

    generic=>conv_to_generic(logical_data, .true.)
    call map_put_generic(specificmap, key, generic, .true.)
  end subroutine map_put_logical
  
  !> Puts a specific key-value pair into the map.
  !!
  !! If the key is not already held in the map then
  !! the key-value pair will be added, otherwise the existing key-value pair will be modified to
  !! hold this updated value (keys must be unique)
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key The key to place in the map
  !! @param data Pointer to the generic data value to place in the map
  !! @param memory_allocation_automatic Whether the collections API should manage the freeing of memory
  subroutine map_put_generic(specificmap, key, data, memory_allocation_automatic)
    type(map_type), intent(inout) :: specificmap
    class(*), pointer, intent(in) :: data
    character(len=*), intent(in) :: key
    logical, intent(in) :: memory_allocation_automatic

    class(*), pointer :: raw_map_node, generic_map_node
    type(mapnode_type), pointer :: newmapnode

    ! Test to see if key already exists in the map
    raw_map_node=>map_getnode(specificmap, key)

    if (associated(raw_map_node)) then
      select type(raw_map_node)
      type is (mapnode_type)
        raw_map_node%value => data
      end select
    else
      allocate(newmapnode)
      newmapnode%value => data
      newmapnode%key = key
      newmapnode%memory_allocation_automatic=memory_allocation_automatic
      ! Clone and deallocate the newmapnode - this keeps GNU happy with passing the correct pointer and Cray
      ! doesn't link the generic pointer just pointing to the data structure hence we clone it
      allocate(generic_map_node, source=newmapnode)
      deallocate(newmapnode)
      call list_add_generic(specificmap%map_ds, generic_map_node, .false.)
    end if
  end subroutine map_put_generic

  !> Determines whether or not a map contains a specific key
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key The key we are looking up
  !! @returns Whether the map contains the key
  logical function map_contains_key(specificmap, key)
    type(map_type), intent(inout) :: specificmap
    character(len=*), intent(in) :: key

    integer :: key_location
    class(*), pointer :: raw_map_node

    raw_map_node => map_getnode(specificmap, key, key_location)
    map_contains_key = key_location .gt. 0
  end function map_contains_key

  !> Retrieves the key currently being held at a specific index in the map or "" if the index > map elements
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i The index to retrieve the key from
  !! @returns The key string
  character(len=STRING_LENGTH) function map_key_at(specificmap, i)
    type(map_type), intent(inout) :: specificmap
    integer, intent(in) :: i

    class(*), pointer :: raw_map_node
    integer :: the_map_size

    the_map_size = map_size(specificmap)
    if (i .le. the_map_size) then
      raw_map_node=>list_get_generic(specificmap%map_ds, i)
      if (associated(raw_map_node)) then
        select type(raw_map_node)
        type is(mapnode_type)
          map_key_at = raw_map_node%key
        end select
        return
      end if
    end if
    map_key_at=""
  end function map_key_at

  !> Retrieves the integer value held at the specific map index
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i Index to get value from
  !! @returns Integer value or raises an error if none is found
  function map_integer_at(specificmap, i)
    type(map_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    integer :: map_integer_at

    class(*), pointer :: generic

    generic=>map_generic_at(specificmap, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find integer at "//trim(conv_to_string(i)))
    map_integer_at=conv_to_integer(generic, .false.)
  end function map_integer_at

  !> Retrieves the string value held at the specific map index
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i Index to get value from
  !! @returns String value or raises an error if none is found
  function map_string_at(specificmap, i)
    type(map_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    character(len=STRING_LENGTH) :: map_string_at

    class(*), pointer :: generic

    generic=>map_generic_at(specificmap, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find string at "//trim(conv_to_string(i)))
    map_string_at=conv_to_string(generic, .false., STRING_LENGTH)
  end function map_string_at

  !> Retrieves the real value held at the specific map index. Converts between precision and int
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i Index to get value from
  !! @returns Double precision real value or raises an error if none is found
  function map_real_at(specificmap, i)
    type(map_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    real(kind=DEFAULT_PRECISION) :: map_real_at

    class(*), pointer :: generic

    generic=>map_generic_at(specificmap, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find real at "//trim(conv_to_string(i)))
    select type(vr=>generic)
    type is (real(kind=DEFAULT_PRECISION))
      map_real_at=vr
    type is (real)
      map_real_at=conv_single_real_to_double(vr)
    type is (integer)  
      map_real_at=conv_single_real_to_double(conv_to_real(vr))
    end select
  end function map_real_at

  !> Retrieves the logical value held at the specific map index
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i Index to get value from
  !! @returns Logical value or raises an error if none is found
  function map_logical_at(specificmap, i)
    type(map_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    logical :: map_logical_at

    class(*), pointer :: generic

    generic=>map_generic_at(specificmap, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find logical at "//trim(conv_to_string(i)))
    map_logical_at=conv_to_logical(generic, .false.)
  end function map_logical_at

  !> Retrieves the generic value held at the specific map index or null if index > map elements
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i Index to get value from
  !! @returns Pointer to the generic value
  function map_generic_at(specificmap, i)
    type(map_type), intent(inout) :: specificmap
    integer, intent(in) :: i

    class(*), pointer :: raw_map_node, map_generic_at
    integer :: the_map_size

    the_map_size = map_size(specificmap)
    if (i .le. the_map_size) then
      raw_map_node=>list_get_generic(specificmap%map_ds, i)
      if (associated(raw_map_node)) then
        select type(raw_map_node)
        type is (mapnode_type)
          map_generic_at => raw_map_node%value
        end select
        return
      end if
    end if
    map_generic_at=>null()
  end function map_generic_at

  !> Retrieves the entry at a specific map index or null if index > map elements
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i Index to get value from
  !! @param key The associated key
  !! @param value Integer value or raises an error if none is found
  logical function map_integer_entry_at(specificmap, i, key, int_val)
    type(map_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    character(len=*), intent(out) :: key
    integer, intent(out) :: int_val

    class(*), pointer :: generic

    map_integer_entry_at=map_generic_entry_at(specificmap, i, key, generic)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find integer entry with key '"//trim(key)//"'")
    int_val=conv_to_integer(generic, .false.)
  end function map_integer_entry_at

  !> Retrieves the entry at a specific map index or null if index > map elements
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i Index to get value from
  !! @param key The associated key
  !! @param value String value or raises an error if none is found
  logical function map_string_entry_at(specificmap, i, key, str_val)
    type(map_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    character(len=*), intent(out) :: key
    character(len=STRING_LENGTH), intent(out) :: str_val

    class(*), pointer :: generic

    map_string_entry_at=map_generic_entry_at(specificmap, i, key, generic)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find string entry with key '"//trim(key)//"'")
    str_val=conv_to_string(generic, .false., STRING_LENGTH)
  end function map_string_entry_at

  !> Retrieves the entry at a specific map index or null if index > map elements. This converts precision and from ints
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i Index to get value from
  !! @param key The associated key
  !! @param value Real value or raises an error if none is found
  logical function map_real_entry_at(specificmap, i, key, real_val)
    type(map_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    character(len=*), intent(out) :: key
    real(kind=DEFAULT_PRECISION), intent(out) :: real_val

    class(*), pointer :: generic

    map_real_entry_at=map_generic_entry_at(specificmap, i, key, generic)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find real entry with key '"//trim(key)//"'")
    select type(vr=>generic)
    type is (real(kind=DEFAULT_PRECISION))
      real_val=vr
    type is (real)
      real_val=conv_single_real_to_double(vr)
    type is (integer)  
      real_val=conv_single_real_to_double(conv_to_real(vr))
    end select
  end function map_real_entry_at

  !> Retrieves the entry at a specific map index or null if index > map elements
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i Index to get value from
  !! @param key The associated key
  !! @param value Logical value or raises an error if none is found
  logical function map_logical_entry_at(specificmap, i, key, logical_val)
    type(map_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    character(len=*), intent(out) :: key
    logical, intent(out) :: logical_val

    class(*), pointer :: generic

    map_logical_entry_at=map_generic_entry_at(specificmap, i, key, generic)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find logical entry with key '"//trim(key)//"'")
    logical_val=conv_to_logical(generic, .false.)
  end function map_logical_entry_at  

  !> Retrieves the entry at a specific map index or null if index > map elements
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i Index to get value from
  !! @param key The associated key
  !! @param value Generic pointer to corresponding value
  logical function map_generic_entry_at(specificmap, i, key, val)
    type(map_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    character(len=*), intent(out) :: key
    class(*), pointer, intent(out) :: val

    class(*), pointer :: raw_map_node
    integer :: the_map_size

    the_map_size = map_size(specificmap)
    if (i .le. the_map_size) then
      raw_map_node=>list_get_generic(specificmap%map_ds, i)
      if (associated(raw_map_node)) then
        select type(raw_map_node)
        type is (mapnode_type)
          val=>raw_map_node%value
          key=raw_map_node%key          
        end select
        map_generic_entry_at=.true.
        return
      end if
    end if
    val=>null()
    map_generic_entry_at=.false.
  end function map_generic_entry_at

  !> Removes a specific key-value pair from the map
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key Key of the key-value pair to remove from the map
  subroutine map_remove(specificmap, key)
    type(map_type), intent(inout) :: specificmap
    character(len=*), intent(in) :: key

    integer :: key_location
    class(*), pointer :: raw_map_node
    type(mapnode_type), pointer :: ptr

    raw_map_node=>map_getnode(specificmap, key, key_location)

    if (key_location .gt. 0) then
      select type (raw_map_node)
      type is (mapnode_type)
        if (raw_map_node%memory_allocation_automatic) then
          if (associated(raw_map_node%value)) deallocate(raw_map_node%value)
        end if
        ptr => raw_map_node
        deallocate(ptr)
      end select
      call list_remove(specificmap%map_ds, key_location)      
    end if
  end subroutine map_remove

  !> Gets a specific element out of the map with the corresponding key
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key Look up key
  !! @returns Integer value associated with the key or raises an error if none is found
  function map_get_int(specificmap, key)
    type(map_type), intent(inout) :: specificmap
    character(len=*), intent(in) :: key
    integer :: map_get_int

    class(*), pointer :: generic

    generic=>map_get_generic(specificmap, key)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find integer entry with key '"//trim(key)//"'")
    map_get_int=conv_to_integer(generic, .false.)
  end function map_get_int

  !> Gets a specific element out of the map with the corresponding key
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key Look up key
  !! @returns String value associated with the key or raises an error if none is found
  function map_get_string(specificmap, key)
    type(map_type), intent(inout) :: specificmap
    character(len=*), intent(in) :: key
    character(len=STRING_LENGTH) :: map_get_string

    class(*), pointer :: generic

    generic=>map_get_generic(specificmap, key)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find string entry with key '"//trim(key)//"'")
    map_get_string=conv_to_string(generic, .false., STRING_LENGTH)
  end function map_get_string

  !> Gets a specific element out of the map with the corresponding key. This converts between precision and from ints
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key Look up key
  !! @returns Real value associated with the key or raises an error if none is found
  function map_get_real(specificmap, key)
    type(map_type), intent(inout) :: specificmap
    character(len=*), intent(in) :: key
    real(kind=DEFAULT_PRECISION) :: map_get_real

    class(*), pointer :: generic

    generic=>map_get_generic(specificmap, key)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find real entry with key '"//trim(key)//"'")
    select type(vr=>generic)
    type is (real(kind=DEFAULT_PRECISION))
      map_get_real=vr
    type is (real)
      map_get_real=conv_single_real_to_double(vr)
    type is (integer)  
      map_get_real=conv_single_real_to_double(conv_to_real(vr))
    end select
  end function map_get_real

  !> Gets a specific element out of the map with the corresponding key
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key Look up key
  !! @returns Logical value associated with the key or raises an error if none is found
  function map_get_logical(specificmap, key)
    type(map_type), intent(inout) :: specificmap
    character(len=*), intent(in) :: key
    logical :: map_get_logical

    class(*), pointer :: generic

    generic=>map_get_generic(specificmap, key)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find logical entry with key '"//trim(key)//"'")
    map_get_logical=conv_to_logical(generic, .false.)
  end function map_get_logical

  !> Gets a specific element out of the map with the corresponding key
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key Look up key
  !! @returns Pointer to the generic value associated with the key or null if none exists
  function map_get_generic(specificmap, key)
    type(map_type), intent(inout) :: specificmap
    character(len=*), intent(in) :: key
    class(*), pointer :: map_get_generic, raw_map_node

    raw_map_node=>map_getnode(specificmap, key)
    if (associated(raw_map_node)) then
      select type (raw_map_node)
      type is (mapnode_type)
        map_get_generic => raw_map_node%value
      end select
      return
    end if
    map_get_generic => null()
  end function map_get_generic

  !> This gets the map node that the key represents (rather than the specific value)
  !!
  !! It allows us to twiddle with map properties but also grab the value out too
  !! Optionally provides us with the list index that the specific node was found at
  !! @param specificmap The specific map involved
  !! @param key Lookup key
  !! @param foundindex (Optional) index where the node is currently held
  !! @returns The pointer to the generic map node data structure
  function map_getnode(specificmap, key, foundindex)
    type(map_type), intent(inout) :: specificmap
    integer, intent(out), optional :: foundindex
    character(len=*), intent(in) :: key
    class(*), pointer :: raw_data, map_getnode

    integer :: i
    type(listnode_type), pointer :: node

    i=1
    node=>specificmap%map_ds%head
    if (associated(node)) then
      do while(1==1)
        raw_data=>node%data
        if (associated(raw_data)) then
          select type (raw_data)
          type is (mapnode_type)
            if (raw_data%key .eq. key) then
              map_getnode=>raw_data
              if (present(foundindex)) foundindex=i
              return
            end if
          end select
        end if
        node=>node%next
        i=i+1
        if (.not. associated(node)) exit
      end do
    end if
    map_getnode => null()
    if(present(foundindex)) foundindex = 0
  end function map_getnode

  !> Returns the number of elements in the map
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @returns Number of key-value pairs held in the map
  integer function map_size(specificmap)
    type(map_type), intent(inout) :: specificmap

    map_size = list_size(specificmap%map_ds)
  end function map_size

  !> Returns whether a map is empty
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @returns Whether the map is empty
  logical function map_is_empty(specificmap)
    type(map_type), intent(inout) :: specificmap

    map_is_empty = list_is_empty(specificmap%map_ds)
  end function map_is_empty

  !> Frees up all the allocatable, heap, memory associated with a specific map.
  !!
  !! This basically acts like a clear operation and once freed the data structure can be reused
  !! Note that it does not free the value at all as this might be referenced else where in the code
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The map to delete all members and free all allocated memory of
  subroutine map_free(specificmap)
    type(map_type), intent(inout) :: specificmap

    type(listnode_type), pointer :: node, previousnode

    node=>specificmap%map_ds%head
    previousnode=>null()

    if (associated(node)) then
      do while(1==1)
        previousnode=>node
        node=>node%next
        if (associated(previousnode%data)) then
          select type (n=>previousnode%data)
          type is (mapnode_type)
            if (n%memory_allocation_automatic) then
              if (associated(n%value)) deallocate(n%value)
            end if
          end select
          deallocate(previousnode%data) ! Free the mapnode data structure
        end if
        deallocate(previousnode)
        if (.not. associated(node)) exit
      end do
    end if

    specificmap%map_ds%tail=>null()
    specificmap%map_ds%head=>null()
    specificmap%map_ds%size=0
  end subroutine map_free

  !> Retrieves an iterator representation of the hashmap, ready to access the first element
  !! @param specificmap Specific collection to base this iterator on
  !! @returns The iterator ready to access the first element
  type(iterator_type) function hashmap_get_iterator(specificmap)
    type(hashmap_type), intent(inout) :: specificmap

    integer :: i

    hashmap_get_iterator%next_item=>null()
    if (associated(specificmap%map_ds)) then
      hashmap_get_iterator%hash_structure=>specificmap%map_ds

      do i=1, size(specificmap%map_ds)
        if (specificmap%map_ds(i)%size .gt. 0) then
          hashmap_get_iterator%next_item=>specificmap%map_ds(i)%head
          exit
        end if
      end do
      hashmap_get_iterator%hash_ptr=i+1
    end if
  end function hashmap_get_iterator 

  !> Puts a specific key-value pair into the hashmap.
  !!
  !! If the key is not already held in the hashmap then
  !! the key-value pair will be added, otherwise the existing key-value pair will be modified to
  !! hold this updated value (keys must be unique). This uses a hashing function for performance
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key The key to place in the map
  !! @param data Integer data value to place in the map
  subroutine hashmap_put_int(specificmap, key, int_data)
    type(hashmap_type), intent(inout) :: specificmap
    integer, intent(in) :: int_data
    character(len=*), intent(in) :: key

    class(*), pointer :: generic

    generic=>conv_to_generic(int_data, .true.)
    call hashmap_put_generic(specificmap, key, generic, .true.)
  end subroutine hashmap_put_int

  !> Puts a specific key-value pair into the hashmap.
  !!
  !! If the key is not already held in the hashmap then
  !! the key-value pair will be added, otherwise the existing key-value pair will be modified to
  !! hold this updated value (keys must be unique). This uses a hashing function for performance
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key The key to place in the map
  !! @param data String data value to place in the map
  subroutine hashmap_put_string(specificmap, key, str_data)
    type(hashmap_type), intent(inout) :: specificmap
    character(len=STRING_LENGTH), intent(in) :: str_data
    character(len=*), intent(in) :: key

    class(*), pointer :: generic

    generic=>conv_to_generic(str_data, .true.)
    call hashmap_put_generic(specificmap, key, generic, .true.)
  end subroutine hashmap_put_string

  !> Puts a specific key-value pair into the hashmap.
  !!
  !! If the key is not already held in the hashmap then
  !! the key-value pair will be added, otherwise the existing key-value pair will be modified to
  !! hold this updated value (keys must be unique). This uses a hashing function for performance
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key The key to place in the map
  !! @param data Double precision real data value to place in the map
  subroutine hashmap_put_real(specificmap, key, real_data)
    type(hashmap_type), intent(inout) :: specificmap
    real(kind=DEFAULT_PRECISION), intent(in) :: real_data
    character(len=*), intent(in) :: key

    class(*), pointer :: generic

    generic=>conv_to_generic(real_data, .true.)
    call hashmap_put_generic(specificmap, key, generic, .true.)
  end subroutine hashmap_put_real

  !> Puts a specific key-value pair into the hashmap.
  !!
  !! If the key is not already held in the hashmap then
  !! the key-value pair will be added, otherwise the existing key-value pair will be modified to
  !! hold this updated value (keys must be unique). This uses a hashing function for performance
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key The key to place in the map
  !! @param data Logical data value to place in the map
  subroutine hashmap_put_logical(specificmap, key, logical_data)
    type(hashmap_type), intent(inout) :: specificmap
    logical, intent(in) :: logical_data
    character(len=*), intent(in) :: key

    class(*), pointer :: generic

    generic=>conv_to_generic(logical_data, .true.)
    call hashmap_put_generic(specificmap, key, generic, .true.)
  end subroutine hashmap_put_logical

  !> Puts a specific key-value pair into the hashmap.
  !!
  !! If the key is not already held in the hashmap then
  !! the key-value pair will be added, otherwise the existing key-value pair will be modified to
  !! hold this updated value (keys must be unique). This uses a hashing function for performance
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key The key to place in the map
  !! @param data Pointer to the generic data value to place in the map
  !! @param memory_allocation_automatic Whether the collections API should manage the freeing of memory
  subroutine hashmap_put_generic(specificmap, key, data, memory_allocation_automatic)
    type(hashmap_type), intent(inout) :: specificmap
    class(*), pointer, intent(in) :: data
    character(len=*), intent(in) :: key
    logical, intent(in) :: memory_allocation_automatic

    class(*), pointer :: raw_map_node, generic_map_node
    type(mapnode_type), pointer :: newmapnode

    if (.not. associated(specificmap%map_ds)) allocate(specificmap%map_ds(hash_size))

    ! Test to see if key already exists in the map
    raw_map_node=>hashmap_getnode(specificmap, key)

    if (associated(raw_map_node)) then
      select type(raw_map_node)
      type is (mapnode_type)
        raw_map_node%value=>data
      end select
    else
      allocate(newmapnode)
      newmapnode%value=>data
      newmapnode%key=key
      newmapnode%memory_allocation_automatic=memory_allocation_automatic
      ! Clone and deallocate the newmapnode - this keeps GNU happy with passing the correct pointer and Cray
      ! doesn't link the generic pointer just pointing to the data structure hence we clone it
      allocate(generic_map_node, source=newmapnode)
      deallocate(newmapnode)
      call list_add_generic(specificmap%map_ds(get_hashkey(key)), generic_map_node, .false.)
      specificmap%size=specificmap%size+1
    end if
  end subroutine hashmap_put_generic

  !> Determines whether or not a hashmap contains a specific key
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key The key we are looking up
  !! @returns Whether the map contains the key
  logical function hashmap_contains_key(specificmap, key)
    type(hashmap_type), intent(inout) :: specificmap
    character(len=*), intent(in) :: key

    class(*), pointer :: raw_map_node

    raw_map_node=>hashmap_getnode(specificmap, key)
    hashmap_contains_key=associated(raw_map_node)
  end function hashmap_contains_key

  !> Retrieves the key currently being held at a specific index in the hashmap or "" if the index > map elements. Note
  !! that this is an expensive operation as it has to potentially process all internal hashed lists so avoid if can
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific hashmap involved
  !! @param i The index to retrieve the key from
  !! @returns The key string
  character(len=STRING_LENGTH) function hashmap_key_at(specificmap, i)
    type(hashmap_type), intent(inout) :: specificmap
    integer, intent(in) :: i

    class(*), pointer :: raw_map_node    

    raw_map_node=>hashmap_getnode_atindex(specificmap, i)
    if (associated(raw_map_node)) then
      select type(raw_map_node)
      type is(mapnode_type)
        hashmap_key_at = raw_map_node%key
      end select
      return
    else
      hashmap_key_at=""
    end if
  end function hashmap_key_at

  !> Retrieves the value held at the specific hashmap index. Note
  !! that this is an expensive operation as it has to potentially process all internal hashed lists so avoid if can
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific hashmap involved
  !! @param i Index to get value from
  !! @returns Integer value or raises an error if none is found
  function hashmap_integer_at(specificmap, i)
    type(hashmap_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    integer :: hashmap_integer_at

    class(*), pointer :: generic

    generic=>hashmap_generic_at(specificmap, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find integer at "//trim(conv_to_string(i)))
    hashmap_integer_at=conv_to_integer(generic, .false.)
  end function hashmap_integer_at

  !> Retrieves the value held at the specific hashmap index. Note
  !! that this is an expensive operation as it has to potentially process all internal hashed lists so avoid if can
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific hashmap involved
  !! @param i Index to get value from
  !! @returns String value or raises an error if none is found
  function hashmap_string_at(specificmap, i)
    type(hashmap_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    character(len=STRING_LENGTH) :: hashmap_string_at

    class(*), pointer :: generic

    generic=>hashmap_generic_at(specificmap, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find string at "//trim(conv_to_string(i)))
    hashmap_string_at=conv_to_string(generic, .false., STRING_LENGTH)
  end function hashmap_string_at

  !> Retrieves the value held at the specific hashmap index. Converts between precision and from int. Note
  !! that this is an expensive operation as it has to potentially process all internal hashed lists so avoid if can
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific hashmap involved
  !! @param i Index to get value from
  !! @returns Double precision real value or raises an error if none is found
  function hashmap_real_at(specificmap, i)
    type(hashmap_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    integer :: hashmap_real_at

    class(*), pointer :: generic

    generic=>hashmap_generic_at(specificmap, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find real at "//trim(conv_to_string(i)))
    select type(vr=>generic)
    type is (real(kind=DEFAULT_PRECISION))
      hashmap_real_at=vr
    type is (real)
      hashmap_real_at=conv_single_real_to_double(vr)
    type is (integer)  
      hashmap_real_at=conv_single_real_to_double(conv_to_real(vr))
    end select
  end function hashmap_real_at

  !> Retrieves the value held at the specific hashmap index. Note
  !! that this is an expensive operation as it has to potentially process all internal hashed lists so avoid if can
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific hashmap involved
  !! @param i Index to get value from
  !! @returns Logical value or raises an error if none is found
  function hashmap_logical_at(specificmap, i)
    type(hashmap_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    logical :: hashmap_logical_at

    class(*), pointer :: generic

    generic=>hashmap_generic_at(specificmap, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find logical at "//trim(conv_to_string(i)))
    hashmap_logical_at=conv_to_logical(generic, .false.)
  end function hashmap_logical_at

  !> Retrieves the value held at the specific hashmap index or null if index > map elements. Note
  !! that this is an expensive operation as it has to potentially process all internal hashed lists so avoid if can
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific hashmap involved
  !! @param i Index to get value from
  !! @returns Pointer to the generic value
  function hashmap_generic_at(specificmap, i)
    type(hashmap_type), intent(inout) :: specificmap
    integer, intent(in) :: i

    class(*), pointer :: raw_map_node, hashmap_generic_at    

    raw_map_node=>hashmap_getnode_atindex(specificmap, i)
    if (associated(raw_map_node)) then
      select type(raw_map_node)
      type is (mapnode_type)
        hashmap_generic_at=>raw_map_node%value
      end select
      return
    else
      hashmap_generic_at=>null()
    end if
  end function hashmap_generic_at

  !> Retrieves the entry at a specific map index
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i Index to get value from
  !! @param key The associated key
  !! @param value Integer value or raises an error if none is found
  logical function hashmap_integer_entry_at(specificmap, i, key, int_val)
    type(hashmap_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    character(len=*), intent(out) :: key
    integer, intent(out) :: int_val

    class(*), pointer :: generic

    hashmap_integer_entry_at=hashmap_generic_entry_at(specificmap, i, key, generic)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find integer entry with key '"//trim(key)//"'")
    int_val=conv_to_integer(generic, .false.)
  end function hashmap_integer_entry_at

  !> Retrieves the entry at a specific map index
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i Index to get value from
  !! @param key The associated key
  !! @param value String value or raises an error if none is found
  logical function hashmap_string_entry_at(specificmap, i, key, str_val)
    type(hashmap_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    character(len=*), intent(out) :: key
    character(len=STRING_LENGTH), intent(out) :: str_val

    class(*), pointer :: generic

    hashmap_string_entry_at=hashmap_generic_entry_at(specificmap, i, key, generic)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find string entry with key '"//trim(key)//"'")
    str_val=conv_to_string(generic, .false., STRING_LENGTH)
  end function hashmap_string_entry_at

  !> Retrieves the entry at a specific map index. This converts between precision and from int
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i Index to get value from
  !! @param key The associated key
  !! @param value Double precision realvalue or raises an error if none is found
  logical function hashmap_real_entry_at(specificmap, i, key, real_val)
    type(hashmap_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    character(len=*), intent(out) :: key
    real(kind=DEFAULT_PRECISION), intent(out) :: real_val

    class(*), pointer :: generic

    hashmap_real_entry_at=hashmap_generic_entry_at(specificmap, i, key, generic)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find real entry with key '"//trim(key)//"'")
    select type(vr=>generic)
    type is (real(kind=DEFAULT_PRECISION))
      real_val=vr
    type is (real)
      real_val=conv_single_real_to_double(vr)
    type is (integer)  
      real_val=conv_single_real_to_double(conv_to_real(vr))
    end select
  end function hashmap_real_entry_at

  !> Retrieves the entry at a specific map index
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i Index to get value from
  !! @param key The associated key
  !! @param value Logical value or raises an error if none is found
  logical function hashmap_logical_entry_at(specificmap, i, key, logical_val)
    type(hashmap_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    character(len=*), intent(out) :: key
    logical, intent(out) :: logical_val

    class(*), pointer :: generic

    hashmap_logical_entry_at=hashmap_generic_entry_at(specificmap, i, key, generic)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find logical entry with key '"//trim(key)//"'")
    logical_val=conv_to_logical(generic, .false.)
  end function hashmap_logical_entry_at  

  !> Retrieves the entry at a specific map index or null if index > map elements
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i Index to get value from
  !! @param key The associated key
  !! @param value Generic pointer to corresponding value
  logical function hashmap_generic_entry_at(specificmap, i, key, val)
    type(hashmap_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    character(len=*), intent(out) :: key
    class(*), pointer, intent(out) :: val

    class(*), pointer :: raw_map_node

    raw_map_node => hashmap_getnode_atindex(specificmap, i)
    if (associated(raw_map_node)) then
      select type(raw_map_node)
      type is (mapnode_type)
        val=>raw_map_node%value
        key=raw_map_node%key          
      end select
      hashmap_generic_entry_at=.true.
      return
    end if
    val=>null()
    hashmap_generic_entry_at=.false.
  end function hashmap_generic_entry_at

  !> Removes a specific key-value pair from the hashmap
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific hashmap involved
  !! @param key Key of the key-value pair to remove from the map
  subroutine hashmap_remove(specificmap, key)
    type(hashmap_type), intent(inout) :: specificmap
    character(len=*), intent(in) :: key

    integer :: key_location
    class(*), pointer :: raw_map_node
    type(mapnode_type), pointer :: ptr

    raw_map_node=>hashmap_getnode(specificmap, key, key_location)
    
    if (key_location .gt. 0) then
      select type (raw_map_node)
      type is (mapnode_type)
        if (raw_map_node%memory_allocation_automatic) then
          if (associated(raw_map_node%value)) deallocate(raw_map_node%value)
        end if
        ptr => raw_map_node
        deallocate(ptr)
      end select
      call list_remove(specificmap%map_ds(get_hashkey(key)), key_location)      
      specificmap%size=specificmap%size-1
    end if
  end subroutine hashmap_remove

  !> Gets a specific element out of the hashmap with the corresponding key
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific hashmap involved
  !! @param key Look up key
  !! @returns Integer value associated with the key or raises an error if none is found
  function hashmap_get_int(specificmap, key)
    type(hashmap_type), intent(inout) :: specificmap
    character(len=*), intent(in) :: key
    integer :: hashmap_get_int

    class(*), pointer :: generic

    generic=>hashmap_get_generic(specificmap, key)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find integer entry with key '"//trim(key)//"'")
    hashmap_get_int=conv_to_integer(generic, .false.)
  end function hashmap_get_int

  !> Gets a specific element out of the hashmap with the corresponding key
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific hashmap involved
  !! @param key Look up key
  !! @returns String value associated with the key or raises an error if none is found
  function hashmap_get_string(specificmap, key)
    type(hashmap_type), intent(inout) :: specificmap
    character(len=*), intent(in) :: key
    character(len=STRING_LENGTH) :: hashmap_get_string

    class(*), pointer :: generic

    generic=>hashmap_get_generic(specificmap, key)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find string entry with key '"//trim(key)//"'")
    hashmap_get_string=conv_to_string(generic, .false., STRING_LENGTH)
  end function hashmap_get_string

  !> Gets a specific element out of the hashmap with the corresponding key. Converts between precision and from int
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific hashmap involved
  !! @param key Look up key
  !! @returns Double precision real value associated with the key or raises an error if none is found
  function hashmap_get_real(specificmap, key)
    type(hashmap_type), intent(inout) :: specificmap
    character(len=*), intent(in) :: key
    real(kind=DEFAULT_PRECISION) :: hashmap_get_real

    class(*), pointer :: generic

    generic=>hashmap_get_generic(specificmap, key)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find real entry with key '"//trim(key)//"'")
    select type(vr=>generic)
    type is (real(kind=DEFAULT_PRECISION))
      hashmap_get_real=vr
    type is (real)
      hashmap_get_real=conv_single_real_to_double(vr)
    type is (integer)  
      hashmap_get_real=conv_single_real_to_double(conv_to_real(vr))
    end select
  end function hashmap_get_real

  !> Gets a specific element out of the hashmap with the corresponding key
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific hashmap involved
  !! @param key Look up key
  !! @returns Logical value associated with the key or raises an error if none is found
  function hashmap_get_logical(specificmap, key)
    type(hashmap_type), intent(inout) :: specificmap
    character(len=*), intent(in) :: key
    logical :: hashmap_get_logical

    class(*), pointer :: generic

    generic=>hashmap_get_generic(specificmap, key)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not find logical entry with key '"//trim(key)//"'")
    hashmap_get_logical=conv_to_logical(generic, .false.)
  end function hashmap_get_logical

  !> Gets a specific element out of the hashmap with the corresponding key
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific hashmap involved
  !! @param key Look up key
  !! @returns Pointer to the generic value associated with the key or null if none exists
  function hashmap_get_generic(specificmap, key)
    type(hashmap_type), intent(inout) :: specificmap
    character(len=*), intent(in) :: key
    class(*), pointer :: hashmap_get_generic, raw_map_node

    raw_map_node=>hashmap_getnode(specificmap, key)
    if (associated(raw_map_node)) then
      select type (raw_map_node)
      type is (mapnode_type)
        hashmap_get_generic=>raw_map_node%value
      end select
      return
    end if
    hashmap_get_generic=>null()
  end function hashmap_get_generic

  !> This gets the hashmap node that the key represents (rather than the specific value)
  !!
  !! It allows us to twiddle with map properties but also grab the value out too
  !! Optionally provides us with the list index that the specific node was found at
  !! @param specificmap The specific hashmap involved
  !! @param key Lookup key
  !! @param foundindex (Optional) index where the node is currently held (in that hash list)
  !! @returns The pointer to the generic map node data structure
  function hashmap_getnode(specificmap, key, key_location)
    type(hashmap_type), intent(inout) :: specificmap
    character(len=*), intent(in) :: key
    integer, intent(out), optional :: key_location
    class(*), pointer :: raw_data, hashmap_getnode

    integer :: i, hash
    type(listnode_type), pointer :: node

    hashmap_getnode=>null()
    if (present(key_location)) key_location=0

    if (.not. associated(specificmap%map_ds)) return

    hash=get_hashkey(key)
    
    i=1
    node=>specificmap%map_ds(hash)%head
    if (associated(node)) then
      do while(1==1)
        raw_data=>node%data
        if (associated(raw_data)) then
          select type (raw_data)
          type is (mapnode_type)
            if (raw_data%key .eq. key) then
              hashmap_getnode=>raw_data
              if (present(key_location)) key_location=i
              return
            end if
          end select
        end if
        node=>node%next
        i=i+1
        if (.not. associated(node)) exit
      end do
    end if
    if (present(key_location)) key_location=0
  end function hashmap_getnode

  !> This gets the hashmap node at a specific index, from the first hash linked list to the end.
  !!  
  !! @param specificmap The specific hashmap involved
  !! @param index Index to locate at
  !! @returns The pointer to the generic map node data structure
  function hashmap_getnode_atindex(specificmap, index)
    type(hashmap_type), intent(inout) :: specificmap
    integer, intent(in) :: index
    class(*), pointer :: hashmap_getnode_atindex

    integer :: i, current_size, prev

    hashmap_getnode_atindex=>null()
    if (.not. associated(specificmap%map_ds) .or. index .gt. specificmap%size) return

    current_size=0
    prev=0
    do i=1, hash_size
      current_size=current_size+list_size(specificmap%map_ds(i))
      if (current_size .ge. index) then
        hashmap_getnode_atindex=>list_get_generic(specificmap%map_ds(i), index-prev)        
        return
      end if
      prev=current_size
    end do
  end function hashmap_getnode_atindex

  !> Returns the number of elements in the hashmap
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific hashmap involved
  !! @returns Number of key-value pairs held in the map
  integer function hashmap_size(specificmap)
    type(hashmap_type), intent(inout) :: specificmap

    hashmap_size=specificmap%size
  end function hashmap_size

  !> Returns whether a hashmap is empty
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific hashmap involved
  !! @returns Whether the map is empty
  logical function hashmap_is_empty(specificmap)
    type(hashmap_type), intent(inout) :: specificmap

    hashmap_is_empty=(specificmap%size == 0)
  end function hashmap_is_empty

  !> Frees up all the allocatable, heap, memory associated with a specific hashmap.
  !!
  !! This basically acts like a clear operation and once freed the data structure can be reused
  !! Note that it does not free the value at all as this might be referenced else where in the code
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The hashmap to delete all members and free all allocated memory of
  subroutine hashmap_free(specificmap)
    type(hashmap_type), intent(inout) :: specificmap

    type(listnode_type), pointer :: node, previousnode
    integer :: i

    if (associated(specificmap%map_ds)) then
      do i=1, hash_size
        node=>specificmap%map_ds(i)%head
        previousnode=>null()

        if (associated(node)) then
          do while(1==1)
            previousnode=>node
            node=>node%next
            if (associated(previousnode%data)) then
              select type (n=>previousnode%data)
              type is (mapnode_type)
                if (n%memory_allocation_automatic) then
                  if (associated(n%value)) deallocate(n%value)
                end if
              end select
              deallocate(previousnode%data) ! Free the mapnode data structure
            end if
            deallocate(previousnode)
            if (.not. associated(node)) exit
          end do
        end if

        specificmap%map_ds(i)%tail=>null()
        specificmap%map_ds(i)%head=>null()
        specificmap%map_ds(i)%size=0
      end do
      specificmap%size=0
      deallocate(specificmap%map_ds)
    end if
  end subroutine hashmap_free

  !> Retrieves an iterator representation of the hashset, ready to access the first element
  !! @param specificset Specific collection to base this iterator on
  !! @returns The iterator ready to access the first element
  type(iterator_type) function hashset_get_iterator(specificset)
    type(hashset_type), intent(inout) :: specificset

    integer :: i

    hashset_get_iterator%next_item=>null()
    if (associated(specificset%set_ds)) then
      hashset_get_iterator%hash_structure=>specificset%set_ds

      do i=1, size(specificset%set_ds)
        if (specificset%set_ds(i)%size .gt. 0) then
          hashset_get_iterator%next_item=>specificset%set_ds(i)%head
          exit
        end if
      end do
      hashset_get_iterator%hash_ptr=i+1
    end if
  end function hashset_get_iterator 

  !> Adds a string to the hashset which stores unique strings, therefore if the string already exists then this is
  !! ignored
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificset The specific set involved
  !! @param key The string key to add to the set
  subroutine hashset_add(specificset, key)
    type(hashset_type), intent(inout) :: specificset
    character(len=*), intent(in) :: key

    class(*), pointer :: generic
    type(setnode_type), pointer :: newsetnode
    integer :: hash, location

    if (.not. associated(specificset%set_ds)) allocate(specificset%set_ds(hash_size))

    call hashset_getlocation(specificset, key, hash, location)

    if (hash .gt. 0 .and. location .eq. 0) then      
      allocate(newsetnode)
      newsetnode%key=key
      ! Clone and deallocate the newmapnode - this keeps GNU happy with passing the correct pointer and Cray
      ! doesn't link the generic pointer just pointing to the data structure hence we clone it
      allocate(generic, source=newsetnode)
      deallocate(newsetnode)
      call list_add_generic(specificset%set_ds(hash), generic, .true.)
      specificset%size=specificset%size+1
    end if
  end subroutine hashset_add

  !> Removes a string from the hashset
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificset The specific ste involved
  !! @param key The string key to remove
  subroutine hashset_remove(specificset, key)
    type(hashset_type), intent(inout) :: specificset
    character(len=*), intent(in) :: key

    integer :: location, hash

    call hashset_getlocation(specificset, key, hash, location)
    if (hash .gt. 0 .and. location .gt. 0) then
      call list_remove(specificset%set_ds(hash), location)
      specificset%size=specificset%size-1
    end if
  end subroutine hashset_remove

  !> Determines wheter the hashset contains a specific key or not
  !!
  !! @param specificset The specific set involved
  !! @param key The string key to test for
  !! @returns Whether the hashset contains this key or not
  logical function hashset_contains(specificset, key)
    type(hashset_type), intent(inout) :: specificset
    character(len=*), intent(in) :: key

    integer :: hash, key_location

    call hashset_getlocation(specificset, key, hash, key_location)
    hashset_contains= (hash .gt. 0 .and. key_location .gt. 0)
  end function hashset_contains  

  !> Determines the location and hash of a key within a specific hashset. The hash is set regardless of whether the key
  !! is found (assuming the set data structure is allocated) and corresponds to the entry in the hash table. The location
  !! is set if the key is found and is the entry in the linked list of a specific hash table entry. This is zero if no
  !! key is in the hash set
  !! @param specificset The specific set involved
  !! @param key The string key to locate
  !! @param hash The hash code of the key
  !! @param key_location Key location (relative to the hash table entry) or zero if no key is found
  subroutine hashset_getlocation(specificset, key, hash, key_location)
    type(hashset_type), intent(inout) :: specificset
    character(len=*), intent(in) :: key
    integer, intent(out) :: hash, key_location
    class(*), pointer :: raw_data

    integer :: i
    type(listnode_type), pointer :: node

    hash=0
    key_location=0

    if (.not. associated(specificset%set_ds)) return

    hash=get_hashkey(key)
    
    i=1
    node=>specificset%set_ds(hash)%head
    if (associated(node)) then
      do while(1==1)
        raw_data=>node%data
        if (associated(raw_data)) then
          select type (raw_data)
          type is (setnode_type)
            if (raw_data%key .eq. key) then            
              key_location=i
              return
            end if
          end select
        end if
        node=>node%next
        i=i+1
        if (.not. associated(node)) exit
      end do
    end if
    key_location=0
  end subroutine hashset_getlocation

  !> Determines whether or not the hashset is empty
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificset The specific set involved
  !! @returns Whether the set is empty or not
  logical function hashset_is_empty(specificset)
    type(hashset_type), intent(in) :: specificset

    hashset_is_empty = specificset%size == 0
  end function hashset_is_empty

  !> Retrieves the key at index i from the set or empty string if index < list size
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificset The specific set involved
  !! @param i Index to look up
  !! @returns The corresponding key at this location or empty string if none is found
  character(len=STRING_LENGTH) function hashset_get_string(specificset, index)
    type(hashset_type), intent(inout) :: specificset
    integer, intent(in) :: index
    class(*), pointer :: generic

    integer :: i, current_size, prev

    hashset_get_string=""
    if (.not. associated(specificset%set_ds) .or. index .gt. specificset%size) return

    current_size=0
    prev=0
    do i=1, hash_size
      current_size=current_size+list_size(specificset%set_ds(i))
      if (current_size .ge. index) then
        generic=>list_get_generic(specificset%set_ds(i), index-prev)        
        if (associated(generic)) then
          select type (generic)
          type is (setnode_type)
            hashset_get_string=generic%key
          end select
          return
        else
          call log_log(LOG_ERROR, "Can not find hashset entry at index "//trim(conv_to_string(index)))
        end if
      end if
      prev=current_size
    end do
  end function hashset_get_string

  !> Frees up all the allocatable, heap, memory associated with a specific set.
  !!
  !! This basically acts like a clear operation and once freed the data structure can be reused
  !! Note that it does not free the data memory at all as this might be referenced else where in the code
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificset The set to delete all members and free all allocated memory of
  subroutine hashset_free(specificset)
    type(hashset_type), intent(inout) :: specificset

    type(listnode_type), pointer :: node, previousnode
    integer :: i

    if (associated(specificset%set_ds)) then
      do i=1, hash_size
        node=>specificset%set_ds(i)%head
        previousnode=>null()

        if (associated(node)) then
          do while(1==1)        
            previousnode=>node
            node=>node%next
            if (associated(previousnode%data)) then            
              deallocate(previousnode%data) ! Free the mapnode data structure
            end if
            deallocate(previousnode)
            if (.not. associated(node)) exit
          end do
        end if

        specificset%set_ds(i)%tail=>null()
        specificset%set_ds(i)%head=>null()
        specificset%set_ds(i)%size=0
      end do
      specificset%size=0
      deallocate(specificset%set_ds)
    end if
  end subroutine hashset_free

  !> Returns the number of elements in a list
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @returns Number of list elements
  integer function hashset_size(specificset)
    type(hashset_type), intent(in) :: specificset

    hashset_size = specificset%size
  end function hashset_size

  !> Translates the string key into a hash from 1 to hash_size (inclusive.) This encoding is deterministic, so will result in the
  !! same hash for a key on multiple calls and is therefore used as the basis of our hashing collections
  !! @param key The key to find the hash code for
  !! @returns The corresponding hash code
  integer function get_hashkey(key)
    character(len=*), intent(in) :: key

    integer :: i

    get_hashkey=5381
    do i=1, len(trim(key))
      get_hashkey=(ishft(get_hashkey,5) + get_hashkey) + ichar(key(i:i))
    end do
    get_hashkey=abs(mod(get_hashkey, hash_size))+1
  end function get_hashkey

  !> Retrieves an iterator representation of the stack, ready to access the first element
  !! @param specificstack Specific collection to base this iterator on
  !! @returns The iterator ready to access the first element
  type(iterator_type) function stack_get_iterator(specificstack)
    type(stack_type), intent(inout) :: specificstack

    stack_get_iterator%next_item=>specificstack%stack_ds%head
    stack_get_iterator%hash_structure=>null()
    stack_get_iterator%hash_ptr=0
  end function stack_get_iterator 

  !> Pushes an element onto the stack (LIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @param data Integer data to push onto the stack
  subroutine stack_push_int(specificstack, int_data)
    type(stack_type), intent(inout) :: specificstack
    integer, intent(in) :: int_data

    class(*), pointer :: generic

    generic=>conv_to_generic(int_data, .true.)
    call stack_push_generic(specificstack, generic, .true.)
  end subroutine stack_push_int

  !> Pushes an element onto the stack (LIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @param data String data to push onto the stack
  subroutine stack_push_string(specificstack, str_data)
    type(stack_type), intent(inout) :: specificstack
    character(len=STRING_LENGTH), intent(in) :: str_data

    class(*), pointer :: generic

    generic=>conv_to_generic(str_data, .true.)
    call stack_push_generic(specificstack, generic, .true.)
  end subroutine stack_push_string

  !> Pushes an element onto the stack (LIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @param data Double precision real data to push onto the stack
  subroutine stack_push_real(specificstack, real_data)
    type(stack_type), intent(inout) :: specificstack
    real(kind=DEFAULT_PRECISION), intent(in) :: real_data

    class(*), pointer :: generic

    generic=>conv_to_generic(real_data, .true.)
    call stack_push_generic(specificstack, generic, .true.)
  end subroutine stack_push_real

  !> Pushes an element onto the stack (LIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @param data Logical data to push onto the stack
  subroutine stack_push_logical(specificstack, logical_data)
    type(stack_type), intent(inout) :: specificstack
    logical, intent(in) :: logical_data

    class(*), pointer :: generic

    generic=>conv_to_generic(logical_data, .true.)
    call stack_push_generic(specificstack, generic, .true.)
  end subroutine stack_push_logical

  !> Pushes an element onto the stack (LIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @param data Pointer to the generic data to push onto the stack
  !! @param memory_allocation_automatic Whether the collections API should manage the freeing of memory
  subroutine stack_push_generic(specificstack, data, memory_allocation_automatic)
    type(stack_type), intent(inout) :: specificstack
    class(*), pointer, intent(in) :: data
    logical, intent(in) :: memory_allocation_automatic

    call list_insert_generic(specificstack%stack_ds, data, 1, memory_allocation_automatic)
  end subroutine stack_push_generic

  !> Pops an element off the stack (LIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @returns Integer data of the last element added to the stack or raises an error if none is found
  function stack_pop_int(specificstack)
    type(stack_type), intent(inout) :: specificstack
    integer :: stack_pop_int

    class(*), pointer :: generic

    generic=>stack_pop_generic(specificstack)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not pop integer from stack")
    stack_pop_int=conv_to_integer(generic, .false.)
  end function stack_pop_int

  !> Pops an element off the stack (LIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @returns String data of the last element added to the stack or raises an error if none is found
  function stack_pop_string(specificstack)
    type(stack_type), intent(inout) :: specificstack
    character(len=STRING_LENGTH) :: stack_pop_string

    class(*), pointer :: generic

    generic=>stack_pop_generic(specificstack)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not pop string from stack")
    stack_pop_string=conv_to_string(generic, .false., STRING_LENGTH)
  end function stack_pop_string

  !> Pops an element off the stack (LIFO). Converts between precision and from int.
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @returns Double precision real data of the last element added to the stack or raises an error if none is found
  function stack_pop_real(specificstack)
    type(stack_type), intent(inout) :: specificstack
    real(kind=DEFAULT_PRECISION) :: stack_pop_real

    class(*), pointer :: generic

    generic=>stack_pop_generic(specificstack)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not pop real from stack")
    select type(vr=>generic)
    type is (real(kind=DEFAULT_PRECISION))
      stack_pop_real=vr
    type is (real)
      stack_pop_real=conv_single_real_to_double(vr)
    type is (integer)  
      stack_pop_real=conv_single_real_to_double(conv_to_real(vr))
    end select
  end function stack_pop_real

  !> Pops an element off the stack (LIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @returns Logical data of the last element added to the stack or raises an error if none is found
  function stack_pop_logical(specificstack)
    type(stack_type), intent(inout) :: specificstack
    logical :: stack_pop_logical

    class(*), pointer :: generic

    generic=>stack_pop_generic(specificstack)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not pop logical from stack")
    stack_pop_logical=conv_to_logical(generic, .false.)
  end function stack_pop_logical  

  !> Pops an element off the stack (LIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @returns Pointer to the generic data of the last element added to the stack
  function stack_pop_generic(specificstack)
    type(stack_type), intent(inout) :: specificstack
    class(*), pointer :: stack_pop_generic

    stack_pop_generic=>stack_get_generic(specificstack, 1)
    call list_remove(specificstack%stack_ds, 1)
  end function stack_pop_generic

  !> Gets a specific element from the stack at index specified
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @param i The index to retrieve the element at
  !! @returns Integer data or raises an error if none is found
  function stack_get_int(specificstack, i)
    type(stack_type), intent(inout) :: specificstack
    integer, intent(in) :: i
    integer :: stack_get_int

    class(*), pointer :: generic

    generic=>stack_get_generic(specificstack, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not get integer from stack at index "//trim(conv_to_string(i)))
    stack_get_int=conv_to_integer(generic, .false.)
  end function stack_get_int

  !> Gets a specific element from the stack at index specified
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @param i The index to retrieve the element at
  !! @returns String data or raises an error if none is found
  function stack_get_string(specificstack, i)
    type(stack_type), intent(inout) :: specificstack
    integer, intent(in) :: i
    character(len=STRING_LENGTH) :: stack_get_string

    class(*), pointer :: generic

    generic=>stack_get_generic(specificstack, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not get string from stack at index "//trim(conv_to_string(i)))
    stack_get_string=conv_to_string(generic, .false., STRING_LENGTH)
  end function stack_get_string

  !> Gets a specific element from the stack at index specified. Converts between precision and from int.
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @param i The index to retrieve the element at
  !! @returns Double precision real data or raises an error if none is found
  function stack_get_real(specificstack, i)
    type(stack_type), intent(inout) :: specificstack
    integer, intent(in) :: i
    real(kind=DEFAULT_PRECISION) :: stack_get_real

    class(*), pointer :: generic

    generic=>stack_get_generic(specificstack, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not get real from stack at index "//trim(conv_to_string(i)))
    select type(vr=>generic)
    type is (real(kind=DEFAULT_PRECISION))
      stack_get_real=vr
    type is (real)
      stack_get_real=conv_single_real_to_double(vr)
    type is (integer)  
      stack_get_real=conv_single_real_to_double(conv_to_real(vr))
    end select
  end function stack_get_real

  !> Gets a specific element from the stack at index specified
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @param i The index to retrieve the element at
  !! @returns Logical data or raises an error if none is found
  function stack_get_logical(specificstack, i)
    type(stack_type), intent(inout) :: specificstack
    integer, intent(in) :: i
    logical :: stack_get_logical

    class(*), pointer :: generic

    generic=>stack_get_generic(specificstack, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not get logical from stack at index "//trim(conv_to_string(i)))
    stack_get_logical=conv_to_logical(generic, .false.)
  end function stack_get_logical

  !> Gets a specific element from the stack at index specified or null if the index > stack size
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @param i The index to retrieve the element at
  !! @returns Pointer to the generic data
  function stack_get_generic(specificstack, i)
    type(stack_type), intent(inout) :: specificstack
    integer, intent(in) :: i
    class(*), pointer :: stack_get_generic

    stack_get_generic=>list_get_generic(specificstack%stack_ds, i)
  end function stack_get_generic

  !> Returns the number of elements held on the stack
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @returns Number of stack elements
  integer function stack_size(specificstack)
    type(stack_type), intent(inout) :: specificstack

    stack_size = list_size(specificstack%stack_ds)
  end function stack_size

  !> Returns whether a stack is empty
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @returns Whether the stack is empty
  logical function stack_is_empty(specificstack)
    type(stack_type), intent(inout) :: specificstack

    stack_is_empty = list_is_empty(specificstack%stack_ds)
  end function stack_is_empty

  !> Frees up all the allocatable, heap, memory associated with a specific stack.
  !!
  !! This basically acts like a clear operation and once freed the data structure can be reused
  !! Note that it does not free the data memory at all as this might be referenced else where in the code
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The stack to delete all members and free all allocated memory of
  subroutine stack_free(specificstack)
    type(stack_type), intent(inout) :: specificstack

    call list_free(specificstack%stack_ds)
  end subroutine stack_free

  !> Retrieves an iterator representation of the queue, ready to access the first element
  !! @param specificqueue Specific collection to base this iterator on
  !! @returns The iterator ready to access the first element
  type(iterator_type) function queue_get_iterator(specificqueue)
    type(queue_type), intent(inout) :: specificqueue

    queue_get_iterator%next_item=>specificqueue%queue_ds%head
    queue_get_iterator%hash_structure=>null()
    queue_get_iterator%hash_ptr=0
  end function queue_get_iterator 

  !> Adds an element to the end of the queue (FIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @param data Integer data to add to the queue
  subroutine queue_push_int(specificqueue, int_data)
    type(queue_type), intent(inout) :: specificqueue
    integer, intent(in) :: int_data

    class(*), pointer :: generic

    generic=>conv_to_generic(int_data, .true.)
    call queue_push_generic(specificqueue, generic, .true.)
  end subroutine queue_push_int

  !> Adds an element to the end of the queue (FIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @param data String data to add to the queue
  subroutine queue_push_string(specificqueue, str_data)
    type(queue_type), intent(inout) :: specificqueue
    character(len=STRING_LENGTH), intent(in) :: str_data

    class(*), pointer :: generic

    generic=>conv_to_generic(str_data, .true.)
    call queue_push_generic(specificqueue, generic, .true.)
  end subroutine queue_push_string

  !> Adds an element to the end of the queue (FIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @param data Double precision real data to add to the queue
  subroutine queue_push_real(specificqueue, real_data)
    type(queue_type), intent(inout) :: specificqueue
    real(kind=DEFAULT_PRECISION), intent(in) :: real_data

    class(*), pointer :: generic

    generic=>conv_to_generic(real_data, .true.)
    call queue_push_generic(specificqueue, generic, .true.)
  end subroutine queue_push_real

  !> Adds an element to the end of the queue (FIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @param data Logical data to add to the queue
  subroutine queue_push_logical(specificqueue, logical_data)
    type(queue_type), intent(inout) :: specificqueue
    logical, intent(in) :: logical_data

    class(*), pointer :: generic

    generic=>conv_to_generic(logical_data, .true.)
    call queue_push_generic(specificqueue, generic, .true.)
  end subroutine queue_push_logical

  !> Adds an element to the end of the queue (FIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @param data Pointer to the generic data to add to the queue
  !! @param memory_allocation_automatic Whether the collections API should manage the freeing of memory
  subroutine queue_push_generic(specificqueue, data, memory_allocation_automatic)
    type(queue_type), intent(inout) :: specificqueue
    class(*), pointer, intent(in) :: data
    logical, intent(in) :: memory_allocation_automatic

    call list_add_generic(specificqueue%queue_ds, data, memory_allocation_automatic)
  end subroutine queue_push_generic

  !> Pops the queue element off the head of the queue (FIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @returns Integer element at head of the queue or raises an error if there is none
  function queue_pop_int(specificqueue)
    type(queue_type), intent(inout) :: specificqueue
    integer :: queue_pop_int

    class(*), pointer :: generic

    generic=>queue_pop_generic(specificqueue)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not pop integer from queue")
    queue_pop_int=conv_to_integer(generic, .false.)
  end function queue_pop_int

  !> Pops the queue element off the head of the queue (FIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @returns String element at head of the queue or raises an error if there is none
  function queue_pop_string(specificqueue)
    type(queue_type), intent(inout) :: specificqueue
    character(len=STRING_LENGTH) :: queue_pop_string

    class(*), pointer :: generic

    generic=>queue_pop_generic(specificqueue)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not pop string from queue")
    queue_pop_string=conv_to_string(generic, .false., STRING_LENGTH)
  end function queue_pop_string

  !> Pops the queue element off the head of the queue (FIFO). Converts between precision and from int.
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @returns Double precision real element at head of the queue or raises an error if there is none
  function queue_pop_real(specificqueue)
    type(queue_type), intent(inout) :: specificqueue
    real(kind=DEFAULT_PRECISION) :: queue_pop_real

    class(*), pointer :: generic

    generic=>queue_pop_generic(specificqueue)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not pop real from queue")
    select type(vr=>generic)
    type is (real(kind=DEFAULT_PRECISION))
      queue_pop_real=vr
    type is (real)
      queue_pop_real=conv_single_real_to_double(vr)
    type is (integer)  
      queue_pop_real=conv_single_real_to_double(conv_to_real(vr))
    end select
  end function queue_pop_real

  !> Pops the queue element off the head of the queue (FIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @returns Logical element at head of the queue or raises an error if there is none
  function queue_pop_logical(specificqueue)
    type(queue_type), intent(inout) :: specificqueue
    logical :: queue_pop_logical

    class(*), pointer :: generic

    generic=>queue_pop_generic(specificqueue)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not pop logical from queue")
    queue_pop_logical=conv_to_logical(generic, .false.)
  end function queue_pop_logical

  !> Pops the queue element off the head of the queue (FIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @returns Pointer to the generic element at head of the queue or null if none
  function queue_pop_generic(specificqueue)
    type(queue_type), intent(inout) :: specificqueue
    class(*), pointer :: queue_pop_generic

    queue_pop_generic=>queue_get_generic(specificqueue, 1)
    call list_remove(specificqueue%queue_ds, 1)
  end function queue_pop_generic

  !> Returns a specific queue element at an index
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @param i The index to look up
  !! @returns Integer element at i or raises an error if none is found
  function queue_get_int(specificqueue, i)
    type(queue_type), intent(inout) :: specificqueue
    integer, intent(in) :: i
    integer :: queue_get_int

    class(*), pointer :: generic

    generic=>queue_get_generic(specificqueue, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not get integer from queue at index "//trim(conv_to_string(i)))
    queue_get_int=conv_to_integer(generic, .false.)
  end function queue_get_int

  !> Returns a specific queue element at an index
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @param i The index to look up
  !! @returns String element at i or raises an error if none is found
  function queue_get_string(specificqueue, i)
    type(queue_type), intent(inout) :: specificqueue
    integer, intent(in) :: i
    character(len=STRING_LENGTH) :: queue_get_string

    class(*), pointer :: generic

    generic=>queue_get_generic(specificqueue, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not get string from queue at index "//trim(conv_to_string(i)))
    queue_get_string=conv_to_string(generic, .false., STRING_LENGTH)
  end function queue_get_string

  !> Returns a specific queue element at an index. Converts between precision and from int.
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @param i The index to look up
  !! @returns Double precision real element at i or raises an error if none is found
  function queue_get_real(specificqueue, i)
    type(queue_type), intent(inout) :: specificqueue
    integer, intent(in) :: i
    real(kind=DEFAULT_PRECISION) :: queue_get_real

    class(*), pointer :: generic

    generic=>queue_get_generic(specificqueue, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not get real from queue at index "//trim(conv_to_string(i)))
    select type(vr=>generic)
    type is (real(kind=DEFAULT_PRECISION))
      queue_get_real=vr
    type is (real)
      queue_get_real=conv_single_real_to_double(vr)
    type is (integer)  
      queue_get_real=conv_single_real_to_double(conv_to_real(vr))
    end select
  end function queue_get_real

  !> Returns a specific queue element at an index
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @param i The index to look up
  !! @returns Logical element at i or raises an error if none is found
  function queue_get_logical(specificqueue, i)
    type(queue_type), intent(inout) :: specificqueue
    integer, intent(in) :: i
    logical :: queue_get_logical

    class(*), pointer :: generic

    generic=>queue_get_generic(specificqueue, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not get logical from queue at index "//trim(conv_to_string(i)))
    queue_get_logical=conv_to_logical(generic, .false.)
  end function queue_get_logical

  !> Returns a specific queue element at an index or null if index > queue size
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @param i The index to look up
  !! @returns Pointer to the generic queue element at i or null if the index does not exist
  function queue_get_generic(specificqueue, i)
    type(queue_type), intent(inout) :: specificqueue
    integer, intent(in) :: i
    class(*), pointer :: queue_get_generic

    queue_get_generic=>list_get_generic(specificqueue%queue_ds, i)
  end function queue_get_generic

  !> Returns the number of elements held in a queue
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @returns Number of elements in the queue
  integer function queue_size(specificqueue)
    type(queue_type), intent(inout) :: specificqueue

    queue_size = list_size(specificqueue%queue_ds)
  end function queue_size

  !> Returns whether a queue is empty
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @returns Whether the queue is empty
  logical function queue_is_empty(specificqueue)
    type(queue_type), intent(inout) :: specificqueue

    queue_is_empty = list_is_empty(specificqueue%queue_ds)
  end function queue_is_empty

  !> Frees up all the allocatable, heap, memory associated with a specific queue.
  !!
  !! This basically acts like a clear operation and once freed the data structure can be reused
  !! Note that it does not free the data memory at all as this might be referenced else where in the code
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The queue to delete all members and free all allocated memory of
  subroutine queue_free(specificqueue)
    type(queue_type), intent(inout) :: specificqueue

    call list_free(specificqueue%queue_ds)
  end subroutine queue_free

  !> Retrieves an iterator representation of the list, ready to access the first element
  !! @param specificlist Specific collection to base this iterator on
  !! @returns The iterator ready to access the first element
  type(iterator_type) function list_get_iterator(specificlist)
    type(list_type), intent(inout) :: specificlist

    list_get_iterator%next_item=>specificlist%head
    list_get_iterator%hash_structure=>null()
    list_get_iterator%hash_ptr=0
  end function list_get_iterator    

  !> Inserts an element into the list or places at the end if the index > list size
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param data Integer data to insert
  !! @param i The list index to insert to
  subroutine list_insert_int(specificlist, int_data, i)
    type(list_type), intent(inout) :: specificlist
    integer, intent(in) :: i, int_data

    class(*), pointer :: generic

    generic=>conv_to_generic(int_data, .true.)
    call list_insert_generic(specificlist, generic, i, .true.)
  end subroutine list_insert_int

  !> Inserts an element into the list or places at the end if the index > list size
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param data String data to insert
  !! @param i The list index to insert to
  subroutine list_insert_string(specificlist, str_data, i)
    type(list_type), intent(inout) :: specificlist
    integer, intent(in) :: i
    character(len=STRING_LENGTH), intent(in) :: str_data

    class(*), pointer :: generic

    generic=>conv_to_generic(str_data, .true.)
    call list_insert_generic(specificlist, generic, i, .true.)
  end subroutine list_insert_string

  !> Inserts an element into the list or places at the end if the index > list size
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param data Double precision real data to insert
  !! @param i The list index to insert to
  subroutine list_insert_real(specificlist, real_data, i)
    type(list_type), intent(inout) :: specificlist
    integer, intent(in) :: i
    real(kind=DEFAULT_PRECISION), intent(in) :: real_data

    class(*), pointer :: generic

    generic=>conv_to_generic(real_data, .true.)
    call list_insert_generic(specificlist, generic, i, .true.)
  end subroutine list_insert_real

  !> Inserts an element into the list or places at the end if the index > list size
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param data Logical data to insert
  !! @param i The list index to insert to
  subroutine list_insert_logical(specificlist, logical_data, i)
    type(list_type), intent(inout) :: specificlist
    integer, intent(in) :: i
    logical, intent(in) :: logical_data

    class(*), pointer :: generic

    generic=>conv_to_generic(logical_data, .true.)
    call list_insert_generic(specificlist, generic, i, .true.)
  end subroutine list_insert_logical

  !> Inserts an element into the list or places at the end if the index > list size
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param data Pointer to the generic data to insert
  !! @param i The list index to insert to
  subroutine list_insert_generic(specificlist, data, i, memory_allocation_automatic)
    type(list_type), intent(inout) :: specificlist
    integer, intent(in) :: i
    class(*), pointer, intent(in) :: data
    logical, intent(in) :: memory_allocation_automatic

    integer ::j
    type(listnode_type), pointer :: newnode, node

    allocate(newnode)
    newnode%data => data

    j=1
    node => specificlist%head
    if (associated(node)) then
      do while(j .lt. i)
        if (.not. associated(node%next)) exit
        node => node%next
        j=j+1
      end do
      if (j .eq. i) then
        ! Insert node
        newnode%next => node
        newnode%prev => node%prev
        newnode%memory_allocation_automatic=memory_allocation_automatic
        if (associated(node%prev)) node%prev%next=>newnode
        node%prev => newnode
        if (associated(node, target=specificlist%head)) specificlist%head=>newnode
      else
        ! Ran out of list nodes so add this one onto the end
        newnode%prev=>specificlist%tail
        if (associated(specificlist%tail)) then
          specificlist%tail%next => newnode
        end if
        specificlist%tail => newnode

        if (associated(specificlist%head) .eqv. .false.) then
          specificlist%head=>newnode
        end if
      end if
    else
      ! No current list data so set up the list with this node
      specificlist%head => newnode
      specificlist%tail => newnode
    end if
    specificlist%size=specificlist%size+1
  end subroutine list_insert_generic

  !> Adds an element to the end of the list
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param data Integer data to add
  subroutine list_add_int(specificlist, int_data)
    type(list_type), intent(inout) :: specificlist
    integer, intent(in) :: int_data

    class(*), pointer :: generic

    generic=>conv_to_generic(int_data, .true.)
    call list_add_generic(specificlist, generic, .true.)
  end subroutine list_add_int

  !> Adds an element to the end of the list
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param data String data to add
  subroutine list_add_string(specificlist, str_data)
    type(list_type), intent(inout) :: specificlist
    character(len=STRING_LENGTH), intent(in) :: str_data

    class(*), pointer :: generic

    generic=>conv_to_generic(str_data, .true.)
    call list_add_generic(specificlist, generic, .true.)
  end subroutine list_add_string

  !> Adds an element to the end of the list
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param data Double precision real data to add
  subroutine list_add_real(specificlist, real_data)
    type(list_type), intent(inout) :: specificlist
    real(kind=DEFAULT_PRECISION), intent(in) :: real_data

    class(*), pointer :: generic

    generic=>conv_to_generic(real_data, .true.)
    call list_add_generic(specificlist, generic, .true.)
  end subroutine list_add_real

  !> Adds an element to the end of the list
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param data Logical data to add
  subroutine list_add_logical(specificlist, logical_data)
    type(list_type), intent(inout) :: specificlist
    logical, intent(in) :: logical_data

    class(*), pointer :: generic

    generic=>conv_to_generic(logical_data, .true.)
    call list_add_generic(specificlist, generic, .true.)
  end subroutine list_add_logical

  !> Adds an element to the end of the list
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param data Pointer to the generic data to add
  !! @param memory_allocation_automatic Whether the collections API should manage the freeing of memory
  subroutine list_add_generic(specificlist, data, memory_allocation_automatic)
    type(list_type), intent(inout) :: specificlist
    class(*), pointer, intent(in) :: data
    logical, intent(in) :: memory_allocation_automatic

    type(listnode_type), pointer :: newnode

    allocate(newnode)
    newnode%data => data

    newnode%prev=>specificlist%tail
    newnode%memory_allocation_automatic=memory_allocation_automatic
    if (associated(specificlist%tail)) then
      specificlist%tail%next => newnode
    end if
    specificlist%tail => newnode

    if (associated(specificlist%head) .eqv. .false.) then
      specificlist%head=>newnode
    end if

    specificlist%size=specificlist%size+1
  end subroutine list_add_generic

  !> Removes an element from the list at a specific index
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param i Index to remove from
  subroutine list_remove(specificlist, i)
    type(list_type), intent(inout) :: specificlist
    integer, intent(in) :: i

    integer ::j
    type(listnode_type), pointer :: node

    j=1
    if (i .le. specificlist%size) then
      node => specificlist%head
      do while(j .lt. i)
        if (.not. associated(node)) exit
        node => node%next
        j=j+1
      end do
      if (associated(node)) then
        if (associated(node%prev)) node%prev%next => node%next
        if (associated(node%next)) node%next%prev => node%prev
        if (associated(node, target=specificlist%head)) specificlist%head => node%next
        if (associated(node, target=specificlist%tail)) specificlist%tail => node%prev
        if (node%memory_allocation_automatic) then
          if (associated(node%data)) deallocate(node%data)
        end if
        deallocate(node)
        specificlist%size = specificlist%size - 1
      end if
    end if
  end subroutine list_remove

  !> Determines whether or not the list is empty
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @returns Whether the list is empty or not
  logical function list_is_empty(specificlist)
    type(list_type), intent(in) :: specificlist

    list_is_empty = specificlist%size == 0
  end function list_is_empty

  !> Retrieves the element at index i from the list
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param i Index to look up
  !! @returns Integer data in the list or raises an error if none is found
  function list_get_int(specificlist, i)
    type(list_type), intent(inout) :: specificlist
    integer, intent(in) :: i
    integer :: list_get_int

    class(*), pointer :: generic

    generic=>list_get_generic(specificlist, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not get integer from list at index "//trim(conv_to_string(i)))
    list_get_int=conv_to_integer(generic, .false.)
  end function list_get_int

  !> Retrieves the element at index i from the list
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param i Index to look up
  !! @returns String data in the list or raises an error if none is found
  function list_get_string(specificlist, i)
    type(list_type), intent(inout) :: specificlist
    integer, intent(in) :: i
    character(len=STRING_LENGTH) :: list_get_string

    class(*), pointer :: generic

    generic=>list_get_generic(specificlist, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not get string from list at index "//trim(conv_to_string(i)))
    list_get_string=conv_to_string(generic, .false., STRING_LENGTH)
  end function list_get_string

  !> Retrieves the element at index i from the list. Converts between precision and from int.
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param i Index to look up
  !! @returns Double precision real data in the list or raises an error if none is found
  function list_get_real(specificlist, i)
    type(list_type), intent(inout) :: specificlist
    integer, intent(in) :: i
    real(kind=DEFAULT_PRECISION) :: list_get_real

    class(*), pointer :: generic

    generic=>list_get_generic(specificlist, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not get real from list at index "//trim(conv_to_string(i)))
    select type(vr=>generic)
    type is (real(kind=DEFAULT_PRECISION))
      list_get_real=vr
    type is (real)
      list_get_real=conv_single_real_to_double(vr)
    type is (integer)  
      list_get_real=conv_single_real_to_double(conv_to_real(vr))
    end select
  end function list_get_real

  !> Retrieves the element at index i from the list
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param i Index to look up
  !! @returns Logical data in the list or raises an error if none is found
  function list_get_logical(specificlist, i)
    type(list_type), intent(inout) :: specificlist
    integer, intent(in) :: i
    logical :: list_get_logical

    class(*), pointer :: generic

    generic=>list_get_generic(specificlist, i)
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not get logical from list at index "//trim(conv_to_string(i)))
    list_get_logical=conv_to_logical(generic, .false.)
  end function list_get_logical

  !> Retrieves the element at index i from the list or null if index < list size
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param i Index to look up
  !! @returns Pointer to the generic data in the list
  function list_get_generic(specificlist, i)
    type(list_type), intent(inout) :: specificlist
    integer, intent(in) :: i
    class(*), pointer :: list_get_generic

    integer :: j
    type(listnode_type), pointer :: node

    j=1
    if (specificlist%size .lt. i) then
      list_get_generic => null()
      return
    end if
    node => specificlist%head
    do while(j .lt. i)
      if (.not. associated(node)) exit
      node => node%next
      j=j+1
    end do
    list_get_generic => node%data
  end function list_get_generic

  !> Frees up all the allocatable, heap, memory associated with a specific list.
  !!
  !! This basically acts like a clear operation and once freed the data structure can be reused
  !! Note that it does not free the data memory at all as this might be referenced else where in the code
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The list to delete all members and free all allocated memory of
  subroutine list_free(specificlist)
    type(list_type), intent(inout) :: specificlist

    type(listnode_type), pointer :: node, previousnode

    node=>specificlist%head
    previousnode=>null()

    if (associated(node)) then
      do while(1==1)
        previousnode=>node
        node=>node%next
        if (previousnode%memory_allocation_automatic) then
          if (associated(previousnode%data)) deallocate(previousnode%data)
        end if
        deallocate(previousnode)
        if (.not. associated(node)) exit
      end do
    end if

    specificlist%tail=>null()
    specificlist%head=>null()
    specificlist%size=0
  end subroutine list_free

  !> Returns the number of elements in a list
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @returns Number of list elements
  integer function list_size(specificlist)
    type(list_type), intent(in) :: specificlist

    list_size = specificlist%size
  end function list_size

  !> Deduces whether an iterator has a next entry or not
  !! @param iterator The iterator to test upon
  !! @returns Whether there is a next entry to access
  logical function iteratior_has_next(iterator)
    type(iterator_type), intent(inout) :: iterator

    iteratior_has_next=associated(iterator%next_item)
  end function iteratior_has_next

  !> Returns the next integer referenced by the iterator and advanced it,  or an error if it has reached the end of iteration
  !! @param iterator The iterator of which to access and advance the next element
  !! @returns The next integer
  integer function iterator_get_next_integer(iterator)
    type(iterator_type), intent(inout) :: iterator

    class(*), pointer :: generic

    generic=>iterator_get_next_generic(iterator)
    if (associated(generic)) then
      iterator_get_next_integer=conv_to_integer(generic, .false.)
    else
      call log_log(LOG_ERROR, "Can not get next integer in iterator as iterator has reached end of collection")
    end if
  end function iterator_get_next_integer
  
  !> Returns the next string referenced by the iterator and advanced it,  or an error if it has reached the end of iteration
  !! @param iterator The iterator of which to access and advance the next element
  !! @returns The next string
  character(len=STRING_LENGTH) function iterator_get_next_string(iterator)
    type(iterator_type), intent(inout) :: iterator

    class(*), pointer :: generic

    generic=>iterator_get_next_generic(iterator)
    if (associated(generic)) then
      select type(generic)
        type is (setnode_type)
          iterator_get_next_string=generic%key
        type is (character(len=*))
          iterator_get_next_string = generic
        class default
          ! Intel compiler complains about the below line
          ! iterator_get_next_string=conv_to_string(generic, .false., STRING_LENGTH)
          ! Workaround to make Intel compiler happy
          iterator_get_next_string = ""
      end select      
    else
      call log_log(LOG_ERROR, "Can not get next string in iterator as iterator has reached end of collection")
    end if
  end function iterator_get_next_string

  !> Returns the next real (double precision) referenced by the iterator and advanced it,  
  !! or an error if it has reached the end of iteration
  !! @param iterator The iterator of which to access and advance the next element
  !! @returns The next double precision real, conversion between single and integers is done automatically
  real(kind=DEFAULT_PRECISION) function iterator_get_next_real(iterator)
    type(iterator_type), intent(inout) :: iterator

    class(*), pointer :: generic

    generic=>iterator_get_next_generic(iterator)
    if (associated(generic)) then
      select type(vr=>generic)
      type is (real(kind=DEFAULT_PRECISION))
        iterator_get_next_real=vr
      type is (real)
        iterator_get_next_real=conv_single_real_to_double(vr)
      type is (integer)  
        iterator_get_next_real=conv_single_real_to_double(conv_to_real(vr))
      end select
    else
      call log_log(LOG_ERROR, "Can not get next real in iterator as iterator has reached end of collection")
    end if
  end function iterator_get_next_real

  !> Returns the next logical referenced by the iterator and advanced it,  or an error if it has reached the end of iteration
  !! @param iterator The iterator of which to access and advance the next element
  !! @returns The next logical
  logical function iterator_get_next_logical(iterator)
    type(iterator_type), intent(inout) :: iterator

    class(*), pointer :: generic

    generic=>iterator_get_next_generic(iterator)
    if (associated(generic)) then
      iterator_get_next_logical=conv_to_logical(generic, .false.)
    else
      call log_log(LOG_ERROR, "Can not get next logical in iterator as iterator has reached end of collection")
    end if
  end function iterator_get_next_logical

  !> Returns the next mapentry referenced by the iterator and advanced it,  or an error if it has reached the end of iteration
  !! or the next item was not a mapentry
  !! @param iterator The iterator of which to access and advance the next element
  !! @returns The next map entry
  function iterator_get_next_mapentry(iterator)
    type(iterator_type), intent(inout) :: iterator
    type(mapentry_type) :: iterator_get_next_mapentry

    class(*), pointer :: generic

    generic=>iterator_get_next_generic(iterator)
    if (associated(generic)) then
      select type(generic)
        type is (mapnode_type)
          iterator_get_next_mapentry%key=generic%key
          iterator_get_next_mapentry%value=>generic%value
        class default
          call log_log(LOG_ERROR, "Next item in iterator is not a map entry")
      end select
    else
      call log_log(LOG_ERROR, "Can not get next map entry in iterator as iterator has reached end of collection")
    end if
  end function iterator_get_next_mapentry  
  
  !> Returns the next generic referenced by the iterator and advanced it, or null if it has reached the end of iteration
  !! @param iterator The iterator of which to access and advance the next element
  !! @returns The next generic or null if none is found
  function iterator_get_next_generic(iterator)
    type(iterator_type), intent(inout) :: iterator
    class(*), pointer :: iterator_get_next_generic

    integer :: i
    
    if (associated(iterator%next_item)) then
      iterator_get_next_generic=>iterator%next_item%data
      iterator%next_item=>iterator%next_item%next
      if (.not. associated(iterator%next_item) .and. associated(iterator%hash_structure) .and. &
           iterator%hash_ptr .le. size(iterator%hash_structure)) then
        do i=iterator%hash_ptr, size(iterator%hash_structure)
          if (iterator%hash_structure(i)%size .gt. 0) then
            iterator%next_item=>iterator%hash_structure(i)%head
            exit
          end if          
        end do
        iterator%hash_ptr=i+1
      end if      
    else
      iterator_get_next_generic=>null()
    end if
  end function iterator_get_next_generic

  !> Retrieves the integer value from a map entry
  !! @param mapentry_item The map entry to retrieve the integer value from
  function mapentry_get_int(mapentry_item)
    type(mapentry_type), intent(in) :: mapentry_item    
    integer :: mapentry_get_int

    class(*), pointer :: generic

    generic=>mapentry_item%value
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not get integer from map entry")
    mapentry_get_int=conv_to_integer(generic, .false.)
  end function mapentry_get_int

  !> Retrieves the string value from a map entry
  !! @param mapentry_item The map entry to retrieve the string value from
  function mapentry_get_string(mapentry_item)
    type(mapentry_type), intent(in) :: mapentry_item
    character(len=STRING_LENGTH) :: mapentry_get_string

    class(*), pointer :: generic

    generic=>mapentry_item%value
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not get string from map entry")
    mapentry_get_string=conv_to_string(generic, .false., STRING_LENGTH)
  end function mapentry_get_string

  !> Retrieves the double precision real value from a map entry
  !! @param mapentry_item The map entry to retrieve the double precision real value from (auto converts from single/ints)
  function mapentry_get_real(mapentry_item)
    type(mapentry_type), intent(in) :: mapentry_item
    real(kind=DEFAULT_PRECISION) :: mapentry_get_real

    class(*), pointer :: generic

    generic=>mapentry_item%value
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not get real from map entry")
    select type(vr=>generic)
    type is (real(kind=DEFAULT_PRECISION))
      mapentry_get_real=vr
    type is (real)
      mapentry_get_real=conv_single_real_to_double(vr)
    type is (integer)  
      mapentry_get_real=conv_single_real_to_double(conv_to_real(vr))
    end select
  end function mapentry_get_real

  !> Retrieves the logical value from a map entry
  !! @param mapentry_item The map entry to retrieve the logical value from
  function mapentry_get_logical(mapentry_item)
    type(mapentry_type), intent(in) :: mapentry_item
    logical :: mapentry_get_logical

    class(*), pointer :: generic

    generic=>mapentry_item%value
    if (.not. associated(generic)) call log_log(LOG_ERROR, "Can not get logical from map entry")
    mapentry_get_logical=conv_to_logical(generic, .false.)
  end function mapentry_get_logical

  !> Retrieves the generic value from a map entry
  !! @param mapentry_item The map entry to retrieve the generic value from
  function mapentry_get_generic(mapentry_item)
    type(mapentry_type), intent(in) :: mapentry_item
    class(*), pointer :: mapentry_get_generic
    
    mapentry_get_generic=>mapentry_item%value
  end function mapentry_get_generic  
end module collections_mod
