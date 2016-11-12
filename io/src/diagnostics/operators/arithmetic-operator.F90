!> The arithmetic operator which allows the user to define arithmetic formulas based on fields and constants 
!! which are then executed BDMAS style. This works by parsing the text forumula into an execution tree which is walked in order
!! to perform the final result. Building the tree is the potentially expensive aspect of this, so built trees are cached as it
!! is likely that the equation will be run many times.
!! Currently this is expressed in terms of scalars, and it will operate on all elements of the data
module arithmetic_operator_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : io_configuration_type, data_values_type, get_data_value_by_field_name
  use data_utils_mod, only : get_action_attribute_string
  use collections_mod, only : hashmap_type, list_type, map_type, hashmap_type, iterator_type, c_add_string, c_get_string, &
       c_get_generic, c_put_generic, c_is_empty, c_get_iterator, c_has_next, c_next_string
  use conversions_mod, only : conv_to_real, conv_is_real, conv_is_integer, conv_to_integer, conv_to_string
  use forthread_mod, only : forthread_rwlock_rdlock, forthread_rwlock_wrlock, forthread_rwlock_unlock, &
       forthread_rwlock_init, forthread_rwlock_destroy
  use threadpool_mod, only : check_thread_status
  use logging_mod, only : LOG_ERROR, log_log
  use optionsdatabase_mod, only : options_has_key, options_get_integer, options_get_string, options_get_real
  implicit none

#ifndef TEST_MODE
  private
#endif

  !< Represents the type of an operator
  integer, parameter :: TERMINAL_OP=0, ADD_OP=1, MINUS_OP=2, MUL_OP=3, DIV_OP=4, MOD_OP=5

  !< A specific node in the execution tree
  type arithmetic_execution_node
     character(len=STRING_LENGTH) :: variable
     integer :: operator
     type(arithmetic_execution_node), pointer :: left, right
  end type arithmetic_execution_node

  !< An equation cache item which olds the execution tree and fields required for that calculation
  type arithmetic_cache_item
     type(arithmetic_execution_node), pointer :: execution_tree
     type(list_type) :: required_fields
  end type arithmetic_cache_item

  type(hashmap_type), volatile :: equation_cache
  integer, volatile :: equation_cache_rwlock

  public initialise_arithmetic_operator, finalise_arithmetic_operator, arithmetic_operator_get_required_fields, &
       perform_arithmetic_operator
contains

  !> Initialises this operator
  subroutine initialise_arithmetic_operator()
    call check_thread_status(forthread_rwlock_init(equation_cache_rwlock, -1))
  end subroutine initialise_arithmetic_operator

  !> Finalises this opertor
  subroutine finalise_arithmetic_operator()
    call check_thread_status(forthread_rwlock_destroy(equation_cache_rwlock))
  end subroutine finalise_arithmetic_operator  

  !> Executes this arithmetic operator by attempting to retrieved the cached equation (and creates one if not found.) If there
  !! is no execution tree it then parses the equation into an execution tree and stores it. The stored execution tree is then
  !! executed and the real result returned
  !! @param io_configuration Configuration of the IO server  
  !! @param field_values The field values
  !! @param action_attributes Attributes associated with the running of this operator
  !! @returns The resulting value
  subroutine perform_arithmetic_operator(io_configuration, field_values, action_attributes, source_monc_location, &
       source_monc, operator_result_values)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(hashmap_type), intent(inout) :: field_values
    type(map_type), intent(inout) :: action_attributes
    integer, intent(in) :: source_monc_location, source_monc
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable, intent(inout) :: operator_result_values

    character(len=STRING_LENGTH) :: equation
    type(arithmetic_cache_item), pointer :: cached_equation
    integer :: data_size

    equation=get_action_attribute_string(action_attributes, "equation")
    cached_equation=>find_or_add_equation(equation)
    if (.not. associated(cached_equation%execution_tree)) then
      cached_equation%execution_tree=>build_equation_tree(io_configuration, equation)
    end if
    data_size=get_size_of_data_being_operated_on(cached_equation, field_values)
    allocate(operator_result_values(data_size))
    operator_result_values=execute_equation_tree(cached_equation%execution_tree, field_values, data_size)
  end subroutine perform_arithmetic_operator

  !> Retrieves the number of data elements that this will operate on. It will produce a log error if any variable lengths are
  !! inconsistent
  !! @param cached_equation The cached equation information
  !! @param field_values The variable value key value pair
  !! @returns The number of data elements being operated on
  integer function get_size_of_data_being_operated_on(cached_equation, field_values)
    type(arithmetic_cache_item), intent(inout) :: cached_equation
    type(hashmap_type), intent(inout) :: field_values

    type(data_values_type), pointer :: variable_data
    type(iterator_type) :: iterator
    integer :: temp_size, prev_size

    get_size_of_data_being_operated_on=-1
    if (.not. c_is_empty(cached_equation%required_fields)) then
      iterator=c_get_iterator(cached_equation%required_fields)
      do while (c_has_next(iterator))
        variable_data=>get_data_value_by_field_name(field_values, c_next_string(iterator))

        if (get_size_of_data_being_operated_on == -1) then
          get_size_of_data_being_operated_on=size(variable_data%values)
        else
          temp_size=size(variable_data%values)
          if (get_size_of_data_being_operated_on .ne. temp_size) then
            if (temp_size .gt. get_size_of_data_being_operated_on) then
              prev_size=get_size_of_data_being_operated_on
              get_size_of_data_being_operated_on=temp_size
              temp_size=prev_size
            end if            
            if (mod(get_size_of_data_being_operated_on, temp_size) .ne. 0) then
              call log_log(LOG_ERROR, &
                   "Can only perform arithmetic on variables with the same array sizes or sizes that divide evenly")
            end if
          end if
        end if
      end do
    end if
  end function get_size_of_data_being_operated_on

  !> Executes an equation tree by doing a post order traversal of the tree. If a node is a terminal then either the variable
  !! is looked up for its value or the constant is returned. If a node is a non terminal then the operator is performed on 
  !! the result value of its left and right subtrees and the value returned.
  !! @param equation_tree The equation tree to traverse
  !! @param field_values The variable value key value pair
  !! @returns A real result from traversing this equation tree
  recursive function execute_equation_tree(equation_tree, field_values, n) result(result_value)
    type(arithmetic_execution_node), pointer, intent(inout) :: equation_tree
    type(hashmap_type), intent(inout) :: field_values
    integer, intent(in) :: n
    real(kind=DEFAULT_PRECISION), dimension(n) :: result_value

    real(kind=DEFAULT_PRECISION), dimension(n) :: left_value, right_value
    type(data_values_type), pointer :: variable_data
    integer :: i

    if (equation_tree%operator==TERMINAL_OP) then
      if (conv_is_real(equation_tree%variable)) then
        result_value=conv_to_real(equation_tree%variable)
      else if (conv_is_integer(equation_tree%variable)) then
        result_value=conv_to_real(conv_to_integer(equation_tree%variable))
      else
        variable_data=>get_data_value_by_field_name(field_values, equation_tree%variable)
        if (size(variable_data%values) .lt. n) then
          do i=1, n, size(variable_data%values)
            result_value(i:i+size(variable_data%values)-1)=variable_data%values
          end do
        else
          result_value=variable_data%values
        end if        
      end if
    else
      left_value=execute_equation_tree(equation_tree%left, field_values, n)
      right_value=execute_equation_tree(equation_tree%right, field_values, n)
      if (equation_tree%operator == ADD_OP) then
        result_value(:)=left_value(:)+right_value(:)
      else if (equation_tree%operator == MINUS_OP) then
        result_value(:)=left_value(:)-right_value(:)
      else if (equation_tree%operator == MUL_OP) then
        result_value(:)=left_value(:)*right_value(:)
      else if (equation_tree%operator == DIV_OP) then
        result_value(:)=left_value(:)/right_value(:)
      else if (equation_tree%operator == MOD_OP) then
        do i=1, n
          result_value(i)=mod(left_value(i), right_value(i))
        end do
      end if
    end if    
  end function execute_equation_tree

  !> Builds the equation tree, this searches for the least significant operator and then splits the equation based upon that.
  !! Each sub equation is then passed to recursive calls of this function which return subtrees for these aspects. Some string
  !! manipulation is done to remove braces which would otherwise be included in the terminal characters
  !! @param io_configuration Configuration of the IO server  
  !! @param equation The equation to represent as a tree
  !! @returns The equation tree which represents the provided equation
  recursive function build_equation_tree(io_configuration, equation) result(equation_tree)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: equation
    type(arithmetic_execution_node), pointer :: equation_tree

    integer :: split_point

    allocate(equation_tree)
    split_point=get_location_of_least_significant_operator(equation)
    if (split_point .gt. 0) then
      equation_tree%operator=get_operator_representation(equation(split_point:split_point))
      equation_tree%left=>build_equation_tree(io_configuration, equation(:split_point-1))
      equation_tree%right=>build_equation_tree(io_configuration, equation(split_point+1:))
    else
      equation_tree%operator=TERMINAL_OP
      equation_tree%variable=equation
      call remove_character(equation_tree%variable, "(")
      call remove_character(equation_tree%variable, ")")
      equation_tree%variable=trim(adjustl(equation_tree%variable))      
      if (equation_tree%variable(1:1) .eq. "{" .and. &
           equation_tree%variable(len_trim(equation_tree%variable):len_trim(equation_tree%variable)) .eq. "}") then
        equation_tree%variable=conv_to_string(options_get_real(&
             io_configuration%options_database, equation_tree%variable(2:len_trim(equation_tree%variable)-1)))    
      end if
    end if
  end function build_equation_tree

  !> Removes all occurances of a character from a string in situ by replacing it with whitespace
  !! @param raw_string The string to process and remove characters from in place
  !! @param c The character to search for and remove (replace by whitespace)
  subroutine remove_character(raw_string, c)
    character(len=*), intent(inout) :: raw_string
    character, intent(in) :: c

    integer :: brace_index

    brace_index=index(raw_string, c)
    do while (brace_index .gt. 0)
      raw_string(brace_index:brace_index)=" "
      brace_index=index(raw_string, c)
    end do
  end subroutine remove_character  

  !> Given an equation this will retrieve the location of the least significant operator in that equation or 0 if no operator
  !! is found (i.e. the string is a terminal.) This takes account of parenthesis.
  integer function get_location_of_least_significant_operator(equation)
    character(len=*), intent(in) :: equation

    integer :: i, eq_len, location_value, op, op_val, brace_level

    get_location_of_least_significant_operator=0
    location_value=999999
    brace_level=0
    eq_len=len(trim(equation))

    do i=eq_len, 1, -1
      if (equation(i:i) == "(") brace_level=brace_level-1
      if (equation(i:i) == ")") brace_level=brace_level+1
      op=get_operator_representation(equation(i:i))
      if (op .gt. -1) then
        op_val=0
        if (op == DIV_OP) op_val=4
        if (op == MOD_OP) op_val=4
        if (op == MUL_OP) op_val=3
        if (op == ADD_OP) op_val=2
        if (op == MINUS_OP) op_val=1
        op_val=op_val + (brace_level*10)
        if (op_val .lt. location_value) then
          location_value=op_val
          get_location_of_least_significant_operator=i
        end if        
      end if      
    end do
  end function get_location_of_least_significant_operator

  !> Given a character representation of an operator this returns the internal numeric type representation of it
  !! @param op_char The character operator representation
  !! @returns The numeric internal type of this operator or -1 if none can be found
  integer function get_operator_representation(op_char)
    character, intent(in) :: op_char

    if (op_char .eq. "/") then 
      get_operator_representation=DIV_OP
    else if (op_char .eq. "*") then
      get_operator_representation=MUL_OP
    else if (op_char .eq. "-") then
      get_operator_representation=MINUS_OP
    else if (op_char .eq. "+") then
      get_operator_representation=ADD_OP
    else if (op_char .eq. "%") then
      get_operator_representation=MOD_OP
    else
      get_operator_representation=-1
    end if
  end function get_operator_representation  

  !> Retrieves the list of fields needed by this operator for a specific configuration
  !! @param action_attributes The attributes which configure the operator
  !! @returns A list of required fields before the operator can run
  type(list_type) function arithmetic_operator_get_required_fields(action_attributes)
    type(map_type), intent(inout) :: action_attributes

    character(len=STRING_LENGTH) :: equation
    type(arithmetic_cache_item), pointer :: cached_equation

    equation=get_action_attribute_string(action_attributes, "equation")
    cached_equation=>find_or_add_equation(equation)
    if (c_is_empty(cached_equation%required_fields)) then
      cached_equation%required_fields=process_equation_to_get_required_fields(equation)
    end if
    arithmetic_operator_get_required_fields=cached_equation%required_fields
  end function arithmetic_operator_get_required_fields

  !> Performs text processing on an equation to extract out all the required variable (fields) needed in order to
  !! run this equation and get the result. Note this ignores all values which are constants (reals or integers)
  !! @param equation Text equation to extract list of required fields (variables) from
  !! @returns A list of variables needed by this equation
  type(list_type) function process_equation_to_get_required_fields(equation)
    character(len=*), intent(in) :: equation

    character(len=STRING_LENGTH) :: str_to_write
    character :: c
    integer :: i, eq_length, starting_len

    eq_length=len(trim(equation))

    starting_len=1
    do i=1, eq_length
      c=equation(i:i)
      if (c .eq. "/" .or. c .eq. "*" .or. c .eq. "-" .or. c .eq. "+" .or. c .eq. "(" .or. c .eq. ")" .or. c .eq. "%") then
        if (starting_len .lt. i) then
          str_to_write=equation(starting_len: i-1)
          if (.not. (conv_is_real(str_to_write) .or. conv_is_integer(str_to_write) .or. str_to_write(1:1) .eq. "{")) then
            call c_add_string(process_equation_to_get_required_fields, str_to_write)
          end if
        end if
        starting_len=i+1
      end if      
    end do
    if (starting_len .le. eq_length) then
      str_to_write=equation(starting_len: i-1)
      if (.not. (conv_is_real(str_to_write) .or. conv_is_integer(str_to_write))) then
        call c_add_string(process_equation_to_get_required_fields, str_to_write)
      end if
    end if
  end function process_equation_to_get_required_fields  

  !> Locates an existing equation in the cache based upon the textual equation representation or creates a new entry and 
  !! returns this one
  !! @param equation Textual equation that we are looking up
  !! @returns Cached entry for the equation
  function find_or_add_equation(equation)
    character(len=*), intent(in) :: equation
    type(arithmetic_cache_item), pointer :: find_or_add_equation

    class(*), pointer :: generic

    find_or_add_equation=>find_equation(equation, .true.)
    if (.not. associated(find_or_add_equation)) then
      call check_thread_status(forthread_rwlock_wrlock(equation_cache_rwlock))
      find_or_add_equation=>find_equation(equation, .false.)
      if (.not. associated(find_or_add_equation)) then
        allocate(find_or_add_equation)
        find_or_add_equation%execution_tree=>null()
        generic=>find_or_add_equation
        call c_put_generic(equation_cache, equation, generic, .false.)
      end if
      call check_thread_status(forthread_rwlock_unlock(equation_cache_rwlock))
    end if
  end function find_or_add_equation

  !> Finds an equation in the cache based upon its textual equation representation or returns null if none is found
  !! @param equation Textual equation that we are looking up
  !! @param dolock Whether to issue a read lock which accessing the collection
  !! @returns Cached entry for the equation or null if none is found
  function find_equation(equation, dolock)
    character(len=*), intent(in) :: equation
    type(arithmetic_cache_item), pointer :: find_equation
    logical, intent(in) :: dolock

    class(*), pointer :: generic

    if (dolock) call check_thread_status(forthread_rwlock_rdlock(equation_cache_rwlock))
    generic=>c_get_generic(equation_cache, equation)
    if (dolock) call check_thread_status(forthread_rwlock_unlock(equation_cache_rwlock))
    if (associated(generic)) then
      select type(generic)
        type is (arithmetic_cache_item)
          find_equation=>generic
      end select
    else
      find_equation=>null()
    end if
  end function find_equation
end module arithmetic_operator_mod
