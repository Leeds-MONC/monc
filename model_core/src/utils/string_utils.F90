!> String utility functionality that is commonly used throughout MONC
module string_utils_mod
  implicit none

#ifndef TEST_MODE
  private
#endif

  public replace_character
contains

  !> Replaces all occurances of a character in a string with another character
  !! @param str The string to replace characters in (this is modified)
  !! @param src_char The source character to look for, which will be replaced
  !! @param tgt_char What to replace source characters with
  subroutine replace_character(str, src_char, tgt_char)
    character(len=*), intent(inout) :: str
    character, intent(in) :: src_char, tgt_char

    integer :: i, n, idx

    n=len_trim(str)
    i=1
    do while (i .le. n)
      idx=index(str(i:n), src_char)
      if (idx == 0) exit
      i=i+idx+1
      str(i-2:i-2)=tgt_char
    end do    
  end subroutine replace_character  
end module string_utils_mod

