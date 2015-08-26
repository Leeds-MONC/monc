!> Collection data structures.
!!
!! The collections utilities which provide common collection structures for different components of
!! MONC. Currently a list, map, stack and queue all with appropriate functionality are provided.
!! The core of all collections is currently a doubly linked list, this is abstracted to allow for the
!! internal structure of collections to change without requiring any user code modifications.
module collections_mod
  use datadefn_mod, only : STRING_LENGTH
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
  end type listnode_type

  !> \private Private map key-value pair data structure
  type mapnode_type
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
     type(list_type), allocatable, dimension(:), private :: map_ds
     integer, private :: size=0
  end type hashmap_type

  !> Hashset structure which will store unique strings. The hashing aspect means that lookup is very fast but it does
  !! add extra memory overhead and the order is non-deterministic
  type, public :: hashset_type
     type(list_type), allocatable, dimension(:), private :: set_ds
     integer, private :: size=0
  end type hashset_type  

  !> Pushes an element onto the stack or queue
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific stack or queue involved
  !! @param data Pointer to the generic data to push onto the collection
  interface c_push
     module procedure stack_push, queue_push
  end interface c_push

  !> Pops an element off the stack or queue
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific stack or queue involved
  !! @returns Pointer to the generic data element which has been popped off the collection
  interface c_pop
     module procedure stack_pop, queue_pop
  end interface c_pop

  !> Adds an element to the end of the list
  !!
  !! This has a time complexity of O(1)
  !! @param collection The specific list involved
  !! @param data Pointer to the raw generic data to add
  interface c_add
     module procedure list_add, hashset_add
  end interface c_add

  !> Inserts an element into the list or places at the end if the index > list size
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific list involved
  !! @param data Pointer to the generic data to insert
  !! @param index The list index to insert to
  interface c_insert
     module procedure list_insert
  end interface c_insert

  !> Puts a specific key-value pair into the map.
  !!
  !! If the key is not already held in the map then
  !! the key-value pair will be added, otherwise the existing key-value pair will be modified to
  !! hold this updated value (keys must be unique)
  !! This has a time complexity of O(n) due to key look up
  !! @param collection The specific map involved
  !! @param key The key to place in the map
  !! @param value Pointer to the generic data value to place in the map
  interface c_put
     module procedure map_put, hashmap_put
  end interface c_put

  !> Gets a specific element out of the list, stack, queue or map with the corresponding key
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific list, stack, queue or map involved
  !! @param key String look up key
  !! @returns Generic pointer to the value associated with the key or null if none exists
  interface c_get
     module procedure list_get, stack_get, queue_get, map_get, hashmap_get, hashset_get
  end interface c_get

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

  !> Retrieves the value held at the specific map index or null if index > map elements
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific map involved
  !! @param index The index to get value from
  !! @returns Generic pointer to the value
  interface c_value_at
     module procedure map_value_at, hashmap_value_at
  end interface c_value_at

  !> Retrieves a map entry at a specific index or null if index > map elements. This is more efficient than calling
  !! key at and then value at (or get with the key) as only requires one search for both the key and value
  !!
  !! This has a time complexity of O(n)
  !! @param collection The specific map involved
  !! @param index The index to get value from
  !! @param key The associated key
  !! @param value Generic pointer to corresponding value
  interface c_entry_at
     module procedure map_entry_at, hashmap_entry_at
  end interface c_entry_at  

  !> Frees up all the allocatable, heap, memory associated with a list, stack, queue or map.
  !!
  !! This basically acts like a clear operation and once freed the data structure can be reused
  !! Note that it does not free the data memory at all as this might be referenced else where in the code
  !! This has a time complexity of O(n)
  !! @param collection The list, stack, queue or map to delete all members and free all allocated memory of
  interface c_free
     module procedure list_free, stack_free, queue_free, map_free, hashmap_free, hashset_free
  end interface c_free

  ! Explicit public interfaces and data items
  public c_push, c_pop, c_add, c_insert, c_put, c_get, c_remove, c_size, c_is_empty, c_contains, &
       c_key_at, c_value_at, c_entry_at, c_free

contains

  !> Puts a specific key-value pair into the map.
  !!
  !! If the key is not already held in the map then
  !! the key-value pair will be added, otherwise the existing key-value pair will be modified to
  !! hold this updated value (keys must be unique)
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key The key to place in the map
  !! @param data Pointer to the generic data value to place in the map
  subroutine map_put(specificmap, key, data)
    type(map_type), intent(inout) :: specificmap
    class(*), pointer, intent(in) :: data
    character(len=*), intent(in) :: key

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
      ! Clone and deallocate the newmapnode - this keeps GNU happy with passing the correct pointer and Cray
      ! doesn't link the generic pointer just pointing to the data structure hence we clone it
      allocate(generic_map_node, source=newmapnode)
      deallocate(newmapnode)
      call list_add(specificmap%map_ds, generic_map_node)
    end if
  end subroutine map_put

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
      raw_map_node => list_get(specificmap%map_ds, i)
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

  !> Retrieves the value held at the specific map index or null if index > map elements
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i Index to get value from
  !! @returns Pointer to the generic value
  function map_value_at(specificmap, i)
    type(map_type), intent(inout) :: specificmap
    integer, intent(in) :: i

    class(*), pointer :: raw_map_node, map_value_at
    integer :: the_map_size

    the_map_size = map_size(specificmap)
    if (i .le. the_map_size) then
      raw_map_node => list_get(specificmap%map_ds, i)
      if (associated(raw_map_node)) then
        select type(raw_map_node)
        type is (mapnode_type)
          map_value_at => raw_map_node%value
        end select
        return
      end if
    end if
    map_value_at=>null()
  end function map_value_at

  !> Retrieves the entry at a specific map index or null if index > map elements
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i Index to get value from
  !! @param key The associated key
  !! @param value Generic pointer to corresponding value
  logical function map_entry_at(specificmap, i, key, val)
    type(map_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    character(len=*), intent(out) :: key
    class(*), pointer, intent(out) :: val

    class(*), pointer :: raw_map_node
    integer :: the_map_size

    the_map_size = map_size(specificmap)
    if (i .le. the_map_size) then
      raw_map_node => list_get(specificmap%map_ds, i)
      if (associated(raw_map_node)) then
        select type(raw_map_node)
        type is (mapnode_type)
          val=>raw_map_node%value
          key=raw_map_node%key          
        end select
        map_entry_at=.true.
        return
      end if
    end if
    val=>null()
    map_entry_at=.false.
  end function map_entry_at

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

    raw_map_node => map_getnode(specificmap, key, key_location)

    if (key_location .gt. 0) then
      call list_remove(specificmap%map_ds, key_location)
    end if
  end subroutine map_remove

  !> Gets a specific element out of the map with the corresponding key
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key Look up key
  !! @returns Pointer to the generic value associated with the key or null if none exists
  function map_get(specificmap, key)
    type(map_type), intent(inout) :: specificmap
    character(len=*), intent(in) :: key
    class(*), pointer :: map_get, raw_map_node

    raw_map_node=>map_getnode(specificmap, key)
    if (associated(raw_map_node)) then
      select type (raw_map_node)
      type is (mapnode_type)
        map_get => raw_map_node%value
      end select
      return
    end if
    map_get => null()
  end function map_get

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
    do while(associated(node))
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
    end do
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

    do while(associated(node))
      previousnode=>node
      node=>node%next
      if (associated(previousnode%data)) then
        deallocate(previousnode%data) ! Free the mapnode data structure
      end if
      deallocate(previousnode)
    end do

    specificmap%map_ds%tail=>null()
    specificmap%map_ds%head=>null()
    specificmap%map_ds%size=0
  end subroutine map_free

  !> Puts a specific key-value pair into the hashmap.
  !!
  !! If the key is not already held in the hashmap then
  !! the key-value pair will be added, otherwise the existing key-value pair will be modified to
  !! hold this updated value (keys must be unique). This uses a hashing function for performance
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param key The key to place in the map
  !! @param data Pointer to the generic data value to place in the map
  subroutine hashmap_put(specificmap, key, data)
    type(hashmap_type), intent(inout) :: specificmap
    class(*), pointer, intent(in) :: data
    character(len=*), intent(in) :: key

    class(*), pointer :: raw_map_node, generic_map_node
    type(mapnode_type), pointer :: newmapnode

    if (.not. allocated(specificmap%map_ds)) allocate(specificmap%map_ds(hash_size))

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
      ! Clone and deallocate the newmapnode - this keeps GNU happy with passing the correct pointer and Cray
      ! doesn't link the generic pointer just pointing to the data structure hence we clone it
      allocate(generic_map_node, source=newmapnode)
      deallocate(newmapnode)
      call list_add(specificmap%map_ds(get_hashkey(key)), generic_map_node)
      specificmap%size=specificmap%size+1
    end if
  end subroutine hashmap_put

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
  !! that this is an expensive operation has it has to potentially process all internal hashed lists so avoid if can
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

  !> Retrieves the value held at the specific hashmap index or null if index > map elements. Note
  !! that this is an expensive operation has it has to potentially process all internal hashed lists so avoid if can
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific hashmap involved
  !! @param i Index to get value from
  !! @returns Pointer to the generic value
  function hashmap_value_at(specificmap, i)
    type(hashmap_type), intent(inout) :: specificmap
    integer, intent(in) :: i

    class(*), pointer :: raw_map_node, hashmap_value_at    

    raw_map_node=>hashmap_getnode_atindex(specificmap, i)
    if (associated(raw_map_node)) then
      select type(raw_map_node)
      type is (mapnode_type)
        hashmap_value_at=>raw_map_node%value
      end select
      return
    else
      hashmap_value_at=>null()
    end if
  end function hashmap_value_at

  !> Retrieves the entry at a specific map index or null if index > map elements
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific map involved
  !! @param i Index to get value from
  !! @param key The associated key
  !! @param value Generic pointer to corresponding value
  logical function hashmap_entry_at(specificmap, i, key, val)
    type(hashmap_type), intent(inout) :: specificmap
    integer, intent(in) :: i
    character(len=*), intent(out) :: key
    class(*), pointer, intent(out) :: val

    class(*), pointer :: raw_map_node
    integer :: the_map_size

    raw_map_node => hashmap_getnode_atindex(specificmap, i)
    if (associated(raw_map_node)) then
      select type(raw_map_node)
      type is (mapnode_type)
        val=>raw_map_node%value
        key=raw_map_node%key          
      end select
      hashmap_entry_at=.true.
      return
    end if
    val=>null()
    hashmap_entry_at=.false.
  end function hashmap_entry_at

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

    raw_map_node=>hashmap_getnode(specificmap, key, key_location)

    if (key_location .gt. 0) then
      call list_remove(specificmap%map_ds(get_hashkey(key)), key_location)
      specificmap%size=specificmap%size-1
    end if
  end subroutine hashmap_remove

  !> Gets a specific element out of the hashmap with the corresponding key
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific hashmap involved
  !! @param key Look up key
  !! @returns Pointer to the generic value associated with the key or null if none exists
  function hashmap_get(specificmap, key)
    type(hashmap_type), intent(inout) :: specificmap
    character(len=*), intent(in) :: key
    class(*), pointer :: hashmap_get, raw_map_node

    raw_map_node=>hashmap_getnode(specificmap, key)
    if (associated(raw_map_node)) then
      select type (raw_map_node)
      type is (mapnode_type)
        hashmap_get=>raw_map_node%value
      end select
      return
    end if
    hashmap_get=>null()
  end function hashmap_get

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

    if (.not. allocated(specificmap%map_ds)) return

    hash=get_hashkey(key)
    
    i=1
    node=>specificmap%map_ds(hash)%head
    do while(associated(node))
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
    end do
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
    class(*), pointer :: raw_data, hashmap_getnode_atindex

    integer :: i, current_size, prev

    hashmap_getnode_atindex=>null()
    if (.not. allocated(specificmap%map_ds) .or. index .gt. specificmap%size) return

    current_size=0
    prev=0
    do i=1, hash_size
      current_size=current_size+list_size(specificmap%map_ds(i))
      if (current_size .ge. index) then
        hashmap_getnode_atindex=>list_get(specificmap%map_ds(i), index-prev)        
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

    integer :: i, hash_size

    hashmap_size=specificmap%size
  end function hashmap_size

  !> Returns whether a hashmap is empty
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificmap The specific hashmap involved
  !! @returns Whether the map is empty
  logical function hashmap_is_empty(specificmap)
    type(hashmap_type), intent(inout) :: specificmap

    integer :: i

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

    if (allocated(specificmap%map_ds)) then
      do i=1, hash_size
        node=>specificmap%map_ds(i)%head
        previousnode=>null()

        do while(associated(node))
          previousnode=>node
          node=>node%next
          if (associated(previousnode%data)) then
            deallocate(previousnode%data) ! Free the mapnode data structure
          end if
          deallocate(previousnode)
        end do

        specificmap%map_ds(i)%tail=>null()
        specificmap%map_ds(i)%head=>null()
        specificmap%map_ds(i)%size=0
      end do
      specificmap%size=0
      deallocate(specificmap%map_ds)
    end if
  end subroutine hashmap_free

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

    if (.not. allocated(specificset%set_ds)) allocate(specificset%set_ds(hash_size))

    call hashset_getlocation(specificset, key, hash, location)

    if (hash .gt. 0 .and. location .eq. 0) then      
      allocate(newsetnode)
      newsetnode%key=key
      ! Clone and deallocate the newmapnode - this keeps GNU happy with passing the correct pointer and Cray
      ! doesn't link the generic pointer just pointing to the data structure hence we clone it
      allocate(generic, source=newsetnode)
      deallocate(newsetnode)
      call list_add(specificset%set_ds(hash), generic)
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

    if (.not. allocated(specificset%set_ds)) return

    hash=get_hashkey(key)
    
    i=1
    node=>specificset%set_ds(hash)%head
    do while(associated(node))
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
    end do
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
  character(len=STRING_LENGTH) function hashset_get(specificset, index)
    type(hashset_type), intent(inout) :: specificset
    integer, intent(in) :: index
    class(*), pointer :: generic

    integer :: i, current_size, prev

    hashset_get=""
    if (.not. allocated(specificset%set_ds) .or. index .gt. specificset%size) return

    current_size=0
    prev=0
    do i=1, hash_size
      current_size=current_size+list_size(specificset%set_ds(i))
      if (current_size .ge. index) then
        generic=>list_get(specificset%set_ds(i), index-prev)        
        if (associated(generic)) then
          select type (generic)
          type is (setnode_type)
            hashset_get=generic%key
          end select
          return
        end if
      end if
      prev=current_size
    end do
  end function hashset_get

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

    if (allocated(specificset%set_ds)) then
      do i=1, hash_size
        node=>specificset%set_ds(i)%head
        previousnode=>null()

        do while(associated(node))
          previousnode=>node
          node=>node%next
          if (associated(previousnode%data)) then
            deallocate(previousnode%data) ! Free the mapnode data structure
          end if
          deallocate(previousnode)
        end do

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

  !> Pushes an element onto the stack (LIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @param data Pointer to the generic data to push onto the stack
  subroutine stack_push(specificstack, data)
    type(stack_type), intent(inout) :: specificstack
    class(*), pointer, intent(in) :: data

    call list_insert(specificstack%stack_ds, data, 1)
  end subroutine stack_push

  !> Pops an element off the stack (LIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @returns Pointer to the generic data of the last element added to the stack
  function stack_pop(specificstack)
    type(stack_type), intent(inout) :: specificstack
    class(*), pointer :: stack_pop

    stack_pop => stack_get(specificstack, 1)
    call list_remove(specificstack%stack_ds, 1)
  end function stack_pop

  !> Gets a specific element from the stack at index specified or null if the index > stack size
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificstack The specific stack involved
  !! @param i The index to retrieve the element at
  !! @returns Pointer to the generic data
  function stack_get(specificstack, i)
    type(stack_type), intent(inout) :: specificstack
    class(*), pointer :: stack_get
    integer, intent(in) :: i

    stack_get => list_get(specificstack%stack_ds, i)
  end function stack_get

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

  !> Adds an element to the end of the queue (FIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @param data Pointer to the generic data to add to the queue
  subroutine queue_push(specificqueue, data)
    type(queue_type), intent(inout) :: specificqueue
    class(*), pointer, intent(in) :: data

    call list_add(specificqueue%queue_ds, data)
  end subroutine queue_push

  !> Pops the queue element off the head of the queue (FIFO)
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @returns Pointer to the generic element at head of the queue or null if none
  function queue_pop(specificqueue)
    type(queue_type), intent(inout) :: specificqueue
    class(*), pointer :: queue_pop

    queue_pop => queue_get(specificqueue, 1)
    call list_remove(specificqueue%queue_ds, 1)
  end function queue_pop

  !> Returns a specific queue element at an index or null if index > queue size
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificqueue The specific queue involved
  !! @param i The index to look up
  !! @returns Pointer to the generic queue element at i or null if the index does not exist
  function queue_get(specificqueue, i)
    type(queue_type), intent(inout) :: specificqueue
    integer, intent(in) :: i
    class(*), pointer :: queue_get

    queue_get => list_get(specificqueue%queue_ds, i)
  end function queue_get

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

  !> Inserts an element into the list or places at the end if the index > list size
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param data Pointer to the generic data to insert
  !! @param i The list index to insert to
  subroutine list_insert(specificlist, data, i)
    type(list_type), intent(inout) :: specificlist
    integer, intent(in) :: i
    class(*), pointer, intent(in) :: data

    integer ::j
    type(listnode_type), pointer :: newnode, node

    allocate(newnode)
    newnode%data => data

    j=1
    node => specificlist%head
    if (associated(node)) then
      do while(associated(node%next) .and. j .lt. i)
        node => node%next
        j=j+1
      end do
      if (j .eq. i) then
        ! Insert node
        newnode%next => node
        newnode%prev => node%prev
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
  end subroutine list_insert

  !> Adds an element to the end of the list
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param data Pointer to the generic data to add
  subroutine list_add(specificlist, data)
    type(list_type), intent(inout) :: specificlist
    class(*), pointer, intent(in) :: data

    type(listnode_type), pointer :: newnode

    allocate(newnode)
    newnode%data => data

    newnode%prev=>specificlist%tail
    if (associated(specificlist%tail)) then
      specificlist%tail%next => newnode
    end if
    specificlist%tail => newnode

    if (associated(specificlist%head) .eqv. .false.) then
      specificlist%head=>newnode
    end if

    specificlist%size=specificlist%size+1
  end subroutine list_add

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
      do while(associated(node) .and. j .lt. i)
        node => node%next
        j=j+1
      end do
      if (associated(node)) then
        if (associated(node%prev)) node%prev%next => node%next
        if (associated(node%next)) node%next%prev => node%prev
        if (associated(node, target=specificlist%head)) specificlist%head => node%next
        if (associated(node, target=specificlist%tail)) specificlist%tail => node%prev
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

  !> Retrieves the element at index i from the list or null if index < list size
  !!
  !! Do not call directly from external module, this is called via the appropriate interface
  !! @param specificlist The specific list involved
  !! @param i Index to look up
  !! @returns Pointer to the generic data in the list
  function list_get(specificlist, i)
    class(*), pointer :: list_get
    type(list_type), intent(inout) :: specificlist
    integer, intent(in) :: i

    integer :: j
    type(listnode_type), pointer :: node

    j=1
    if (specificlist%size .lt. i) then
      list_get => null()
      return
    end if
    node => specificlist%head
    do while(associated(node) .and. j .lt. i)
      node => node%next
      j=j+1
    end do
    list_get => node%data
  end function list_get

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

    do while(associated(node))
      previousnode=>node
      node=>node%next
      deallocate(previousnode)
    end do

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
end module collections_mod
