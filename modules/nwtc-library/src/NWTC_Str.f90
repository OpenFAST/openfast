!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!    This file is part of the NWTC Subroutine Library.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!**********************************************************************************************************************************

!> This module contains string manipulation routines
MODULE NWTC_Str

   use Precision  ! ProgDesc and other types with copy and other routines for those types

   implicit  none

   interface is_numeric
      module procedure is_numericR4
      module procedure is_numericR8
   end interface

CONTAINS



!> Count number of occurence of a substring in an input string. 
function countsubstring(s1, s2) result(c)
   character(len=*), intent(in) :: s1 !< Input string
   character(len=*), intent(in) :: s2 !< string to be searched
   integer :: c !< number of substrings
   integer ::  p, posn
   c = 0
   if(len(s2) == 0) return
   p = 1
   do 
      posn = index(s1(p:), s2)
      if(posn == 0) return
      c = c + 1
      p = p + posn + len(s2)
   end do
end function countsubstring

!> split a string according to a delimiter of size 1
subroutine strsplit(String, StrArray, delimiter)
   character(len=*),intent(in)    :: String
   character(len=1),intent(in)    :: delimiter
   character(1024), intent(out), allocatable :: StrArray(:) ! Array of strings extracted from line
   ! Variable
   integer :: j, k, l, n, nmax
   logical :: EndOfLine
   ! Find number of occurences
   n    = countsubstring(String, delimiter)
   nmax = n+1
   ! Allocate substrings
   if(allocated(StrArray)) deallocate(StrArray)
   allocate(StrArray(nmax))
   StrArray(:)=''
   ! Loop on string and store splits
   n = 0
   k = 1
   l = len_trim(string)
   EndOfLine = l-k < 0
   do while (.not.EndOfLine)
      j = index(string(k:l),delimiter)
      if (j == 0) then
         j = l + 1
      else
         j = j + k - 1
      end if
      n = n + 1
      if(n==nmax) then
         StrArray(n) = String(k:len(String))
         EndOfLine = .true.
      else
         if (j /= k .and. len_trim(string(k:j-1)) /= 0) StrArray(n) = String(k:j-1)
         k = j + 1
         EndOfLine = l-k < 0
      endif
   end do
end subroutine strsplit

!> Return true if string is an integer, and also return the integer
logical function is_integer(string, x)
   character(len=*), intent(in ) :: string
   integer(IntKi),   intent(out) :: x
   integer :: e, n
   x = 0
   n=len_trim(string)
   if (n==0) then ! blank lines shouldn't be valid integers
      is_integer = .false.
      return
   end if
   read(string,*,iostat=e) x
   is_integer = e == 0
end function is_integer

logical function is_numericR4(string, x) result(is_numeric)
   character(len=*), intent(in ) :: string
   real(SiKi),       intent(out) :: x
   integer :: e,n
   character(len=12) :: fmt
   x = 0.0_ReKi
   n=len_trim(string)
   write(fmt,'("(F",I0,".0)")') n
   read(string,fmt,iostat=e) x
   is_numeric = e == 0
end function is_numericR4

logical function is_numericR8(string, x) result(is_numeric)
   character(len=*), intent(in ) :: string
   real(R8Ki),       intent(out) :: x
   integer :: e,n
   character(len=12) :: fmt
   x = 0.0_ReKi
   n=len_trim(string)
   write(fmt,'("(F",I0,".0)")') n
   read(string,fmt,iostat=e) x
   is_numeric = e == 0
end function is_numericR8

logical function is_logical(string, b)
   character(len=*), intent(in ) :: string
   logical,          intent(out) :: b
   integer :: e,n
   b = .false.
   n=len_trim(string)
   read(string,*,iostat=e) b
   is_logical = e == 0
end function is_logical

END MODULE NWTC_Str
