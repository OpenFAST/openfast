!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!    This file is part of NWTC_Library.
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
!
!**********************************************************************************************************************************
module YAML
   use Precision, only: IntKi, SiKi, R8Ki
   use NWTC_Base, only: ErrID_None, ErrID_Fatal
   use NWTC_IO, only: num2lstr

   implicit none

   integer(IntKi), parameter :: INDENT_SPACES = 2

   !> Write 1D or 2D array to file
   interface yaml_write_array
      module procedure yaml_write_array1I   ! Single dimension array (Ary) of IntKi
      module procedure yaml_write_array1R4  ! Single dimension array (Ary) of SiKi
      module procedure yaml_write_array2R4  ! Two dimension array of SiKi
      module procedure yaml_write_array1R8  ! Single dimension array (Ary) of R8Ki
      module procedure yaml_write_array2R8  ! Two dimension array of R8Ki
      module procedure yaml_write_array2I   ! Two dimension array of IntKi
   end interface

   !> Write list to file
   interface yaml_write_list
      module procedure yaml_write_listI   ! Single dimension array (Ary) of IntKi
      module procedure yaml_write_listR4  ! Single dimension array (Ary) of SiKi
      module procedure yaml_write_listR8  ! Single dimension array (Ary) of R8Ki
   end interface
   
   !> Write variable to file
   interface yaml_write_var
      module procedure yaml_write_varC   ! Character
      module procedure yaml_write_varI   ! IntKi
      module procedure yaml_write_varR4  ! SiKi
      module procedure yaml_write_varR8  ! R8Ki
   end interface
   private

   public :: yaml_write_comm
   public :: yaml_write_var
   public :: yaml_write_list
   public :: yaml_write_array

contains

!> Write comment to yaml file
subroutine yaml_write_comm(fid, comment, ErrStat, ErrMsg)
   integer(IntKi),             intent(in   ) :: fid     !< File Unit
   character(len=*),           intent(in   ) :: comment !< Comment
   integer,                    intent(  out) :: ErrStat !< A non-zero value indicates an error occurred
   character(len=*),           intent(  out) :: ErrMsg  !< Error message if errstat /= errid_none
   ErrStat = ErrID_None   
   ErrMsg  = ""
   write(fid, '("#",A)', iostat=ErrStat) comment
   if (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error writing comment `'//trim(comment)//'` to YAML file'
   endif
end subroutine yaml_write_comm
! --------------------------------------------------------------------------------
! --- Write variable 
! --------------------------------------------------------------------------------
!> Write simple key/variable to yaml file
subroutine yaml_write_varC(fid, key, val, VarFmt, ErrStat, ErrMsg, level, comment)
   integer(IntKi),             intent(in   ) :: fid     !< File Unit
   character(len=*),           intent(in   ) :: key     !< Variable name
   character(len=*),           intent(in   ) :: val     !< Value
   character(len=*),           intent(in   ) :: VarFmt  !< Format for printing real numbers  
   integer,                    intent(  out) :: ErrStat !< A non-zero value indicates an error occurred
   character(len=*),           intent(  out) :: ErrMsg  !< Error message if errstat /= errid_none
   integer(IntKi),   optional, intent(in   ) :: level   !< indentation level
   character(len=*), optional, intent(in   ) :: comment !< 
   character(256)     :: Fmt
   ErrStat = ErrID_None   
   ErrMsg  = ""
   Fmt = ''
   if (present(level)) Fmt = trim(Num2LStr(level*INDENT_SPACES))//'X,'
   if (present(comment)) then
      Fmt = '('//trim(Fmt)//'A,": ",'//trim(VarFmt)//', " # ",A)'   
      write(fid, Fmt, iostat=ErrStat) key, val, comment
   else
      Fmt = '('//trim(Fmt)//'A,": ",'//trim(VarFmt)//')'   
      write(fid, Fmt, iostat=ErrStat) key, val
   endif
   if (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error writing variable '//trim(key)//' to YAML file'
   endif
end subroutine yaml_write_varC

subroutine yaml_write_varI(fid, key, val, VarFmt, ErrStat, ErrMsg, level, comment)
   integer(IntKi),             intent(in   ) :: fid     !< File Unit
   character(len=*),           intent(in   ) :: key     !< Variable name
   integer(IntKi),             intent(in   ) :: val     !< Value
   character(len=*),           intent(in   ) :: VarFmt  !< Format for printing real numbers  
   integer,                    intent(  out) :: ErrStat !< A non-zero value indicates an error occurred
   character(len=*),           intent(  out) :: ErrMsg  !< Error message if errstat /= errid_none
   integer(IntKi),   optional, intent(in   ) :: level   !< indentation level
   character(len=*), optional, intent(in   ) :: comment !< 
   character(256)     :: Fmt
   ErrStat = ErrID_None   
   ErrMsg  = ""
   Fmt = ''
   if (present(level)) Fmt = trim(Num2LStr(level*INDENT_SPACES))//'X,'
   if (present(comment)) then
      Fmt = '('//trim(Fmt)//'A,": ",'//trim(VarFmt)//', " # ",A)'   
      write(fid, Fmt, iostat=ErrStat) key, val, comment
   else
      Fmt = '('//trim(Fmt)//'A,": ",'//trim(VarFmt)//')'   
      write(fid, Fmt, iostat=ErrStat) key, val
   endif
   if (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error writing variable '//trim(key)//' to YAML file'
   endif
end subroutine yaml_write_varI

subroutine yaml_write_varR4(fid, key, val, VarFmt, ErrStat, ErrMsg, level, comment)
   integer(IntKi),             intent(in   ) :: fid     !< File Unit
   character(len=*),           intent(in   ) :: key     !< Variable name
   real(SiKi),                 intent(in   ) :: val     !< Value
   character(len=*),           intent(in   ) :: VarFmt  !< Format for printing real numbers  
   integer,                    intent(  out) :: ErrStat !< A non-zero value indicates an error occurred
   character(len=*),           intent(  out) :: ErrMsg  !< Error message if errstat /= errid_none
   integer(IntKi),   optional, intent(in   ) :: level   !< indentation level
   character(len=*), optional, intent(in   ) :: comment !< 
   character(256)     :: Fmt
   ErrStat = ErrID_None   
   ErrMsg  = ""
   Fmt = ''
   if (present(level)) Fmt = trim(Num2LStr(level*INDENT_SPACES))//'X,'
   if (present(comment)) then
      Fmt = '('//trim(Fmt)//'A,": ",'//trim(VarFmt)//', " # ",A)'   
      write(fid, Fmt, iostat=ErrStat) key, val, comment
   else
      Fmt = '('//trim(Fmt)//'A,": ",'//trim(VarFmt)//')'   
      write(fid, Fmt, iostat=ErrStat) key, val
   endif
   if (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error writing variable '//trim(key)//' to YAML file'
   endif
end subroutine yaml_write_varR4

subroutine yaml_write_varR8(fid, key, val, VarFmt, ErrStat, ErrMsg, level, comment)
   integer(IntKi),             intent(in   ) :: fid     !< File Unit
   character(len=*),           intent(in   ) :: key     !< Variable name
   real(R8Ki),                 intent(in   ) :: val     !< Value
   character(len=*),           intent(in   ) :: VarFmt  !< Format for printing real numbers  
   integer,                    intent(  out) :: ErrStat !< A non-zero value indicates an error occurred
   character(len=*),           intent(  out) :: ErrMsg  !< Error message if errstat /= errid_none
   integer(IntKi),   optional, intent(in   ) :: level   !< indentation level
   character(len=*), optional, intent(in   ) :: comment !< 
   character(256)     :: Fmt
   ErrStat = ErrID_None   
   ErrMsg  = ""
   Fmt = ''
   if (present(level)) Fmt = trim(Num2LStr(level*INDENT_SPACES))//'X,'
   if (present(comment)) then
      Fmt = '('//trim(Fmt)//'A,": ",'//trim(VarFmt)//', " # ",A)'   
      write(fid, Fmt, iostat=ErrStat) key, val, comment
   else
      Fmt = '('//trim(Fmt)//'A,": ",'//trim(VarFmt)//')'   
      write(fid, Fmt, iostat=ErrStat) key, val
   endif
   if (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error writing variable '//trim(key)//' to YAML file'
   endif
end subroutine yaml_write_varR8


! --------------------------------------------------------------------------------}
! --- Write List  
! --------------------------------------------------------------------------------{
subroutine yaml_write_listI(fid, key, A, VarFmt, ErrStat, ErrMsg, level, comment)
   integer(IntKi),             intent(in   ) :: fid     !< File Unit
   character(len=*),           intent(in   ) :: key     !< list name
   integer(IntKi), dimension(:), intent(in   ) :: A     !< list
   character(len=*),           intent(in   ) :: VarFmt  !< Format for printing real numbers  
   integer,                    intent(  out) :: ErrStat !< A non-zero value indicates an error occurred
   character(len=*),           intent(  out) :: ErrMsg  !< Error message if errstat /= errid_none
   integer(IntKi),   optional, intent(in   ) :: level   !< indentation level
   character(len=*), optional, intent(in   ) :: comment !< 
   integer            :: nc  ! size (rows and columns) of A
   character(256)     :: Fmt
   ErrStat = ErrID_None   
   ErrMsg  = ""
   nc = size(A,1)
   Fmt = ''
   if (present(level)) Fmt = trim(Num2LStr(level*INDENT_SPACES))//'X,'
   if (nc==0) then
      if (present(comment)) then
         write(fid, '('//trim(Fmt)//'A,": []",1X,"#",1X,A)', iostat=ErrStat) trim(key), trim(comment)
      else
         write(fid, '('//trim(Fmt)//'A,": []")', iostat=ErrStat) trim(key)
      endif
   else
      if (nc==1) then
         Fmt = '('//trim(Fmt)//'A,": [",'//VarFmt//',"]"'
      else
         Fmt = '('//trim(Fmt)//'A,": [",'//trim(Num2LStr(nc-1))//'('//VarFmt//',","), '//VarFmt//'"]"' 
      endif
      if (present(comment)) then
         Fmt = trim(Fmt)//',1X,"#",1X,A)' 
         write(fid, Fmt, iostat=ErrStat)trim(key), A(:), trim(comment)
      else
         Fmt = trim(Fmt)//')' 
         write(fid, Fmt, iostat=ErrStat)trim(key), A(:)
      endif
   endif
   if (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error writing list '//trim(key)//' to YAML file'
   end if
end subroutine yaml_write_listI

subroutine yaml_write_listR4(fid, key, A, VarFmt, ErrStat, ErrMsg, level, comment)
   integer(IntKi),             intent(in   ) :: fid     !< File Unit
   character(len=*),           intent(in   ) :: key     !< list name
   real(SiKi), dimension(:),   intent(in   ) :: A       !< list
   character(len=*),           intent(in   ) :: VarFmt  !< Format for printing real numbers  
   integer,                    intent(  out) :: ErrStat !< A non-zero value indicates an error occurred
   character(len=*),           intent(  out) :: ErrMsg  !< Error message if errstat /= errid_none
   integer(IntKi),   optional, intent(in   ) :: level   !< indentation level
   character(len=*), optional, intent(in   ) :: comment !< 
   integer            :: nc  ! size (rows and columns) of A
   character(256)     :: Fmt
   ErrStat = ErrID_None   
   ErrMsg  = ""
   nc = size(A,1)
   Fmt = ''
   if (present(level)) Fmt = trim(Num2LStr(level*INDENT_SPACES))//'X,'
   if (nc==0) then
      if (present(comment)) then
         write(fid, '('//trim(Fmt)//'A,": []",1X,"#",1X,A)', iostat=ErrStat) trim(key), trim(comment)
      else
         write(fid, '('//trim(Fmt)//'A,": []")', iostat=ErrStat) trim(key)
      endif
   else
      if (nc==1) then
         Fmt = '('//trim(Fmt)//'A,": [",'//VarFmt//',"]"'
      else
         Fmt = '('//trim(Fmt)//'A,": [",'//trim(Num2LStr(nc-1))//'('//VarFmt//',","), '//VarFmt//'"]"' 
      endif
      if (present(comment)) then
         Fmt = trim(Fmt)//',1X,"#",1X,A)' 
         write(fid, Fmt, iostat=ErrStat)trim(key), A(:), trim(comment)
      else
         Fmt = trim(Fmt)//')' 
         write(fid, Fmt, iostat=ErrStat)trim(key), A(:)
      endif
   endif
   if (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error writing list '//trim(key)//' to YAML file'
   end if
end subroutine yaml_write_listR4

subroutine yaml_write_listR8(fid, key, A, VarFmt, ErrStat, ErrMsg, level, comment)
   integer(IntKi),             intent(in   ) :: fid     !< File Unit
   character(len=*),           intent(in   ) :: key     !< list name
   real(R8Ki), dimension(:),   intent(in   ) :: A       !< list
   character(len=*),           intent(in   ) :: VarFmt  !< Format for printing real numbers  
   integer,                    intent(  out) :: ErrStat !< A non-zero value indicates an error occurred
   character(len=*),           intent(  out) :: ErrMsg  !< Error message if errstat /= errid_none
   integer(IntKi),   optional, intent(in   ) :: level   !< indentation level
   character(len=*), optional, intent(in   ) :: comment !< 
   integer            :: nc  ! size (rows and columns) of A
   character(256)     :: Fmt
   ErrStat = ErrID_None   
   ErrMsg  = ""
   nc = size(A,1)
   Fmt = ''
   if (present(level)) Fmt = trim(Num2LStr(level*INDENT_SPACES))//'X,'
   if (nc==0) then
      if (present(comment)) then
         write(fid, '('//trim(Fmt)//'A,": []",1X,"#",1X,A)', iostat=ErrStat) trim(key), trim(comment)
      else
         write(fid, '('//trim(Fmt)//'A,": []")', iostat=ErrStat) trim(key)
      endif
   else
      if (nc==1) then
         Fmt = '('//trim(Fmt)//'A,": [",'//VarFmt//',"]"'
      else
         Fmt = '('//trim(Fmt)//'A,": [",'//trim(Num2LStr(nc-1))//'('//VarFmt//',","), '//VarFmt//'"]"' 
      endif
      if (present(comment)) then
         Fmt = trim(Fmt)//',1X,"#",1X,A)' 
         write(fid, Fmt, iostat=ErrStat)trim(key), A(:), trim(comment)
      else
         Fmt = trim(Fmt)//')' 
         write(fid, Fmt, iostat=ErrStat)trim(key), A(:)
      endif
   endif
   if (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error writing list '//trim(key)//' to YAML file'
   end if
end subroutine yaml_write_listR8


! --------------------------------------------------------------------------------
! --- Write 1D array
! --------------------------------------------------------------------------------
subroutine yaml_write_array1I(fid, key, A, VarFmt, ErrStat, ErrMsg, level, comment)
   integer(IntKi),             intent(in   ) :: fid     !< File Unit
   character(len=*),           intent(in   ) :: key     !< Array name
   integer(IntKi), dimension(:), intent(in   ) :: A     !< Array
   character(len=*),           intent(in   ) :: VarFmt  !< Format for printing real numbers  
   integer,                    intent(  out) :: ErrStat !< A non-zero value indicates an error occurred
   character(len=*),           intent(  out) :: ErrMsg  !< Error message if errstat /= errid_none
   integer(IntKi),   optional, intent(in   ) :: level   !< indentation level
   character(len=*), optional, intent(in   ) :: comment !< 
   integer            :: nc  ! size (rows and columns) of A
   character(256)     :: Fmt
   ErrStat = ErrID_None   
   ErrMsg  = ""
   nc = size(A,1)
   Fmt = ''
   if (present(level)) Fmt = trim(Num2LStr(level*INDENT_SPACES))//'X,'
   if (present(comment)) then
      write(fid, '('//trim(Fmt)//'A,": # ",I0," x ",I0,1X,A)', iostat=ErrStat ) trim(key), 1, nc, trim(comment)
   else
      write(fid, '('//trim(Fmt)//'A,": # ",I0," x ",I0)'     , iostat=ErrStat ) trim(key), 1, nc
   end if
   if (present(level)) then
      Fmt = trim(Num2LStr((level+1)*INDENT_SPACES))//'X,'
   else
      Fmt = trim(Num2LStr(INDENT_SPACES))//'X,'
   endif
   if (nc==0) then
      write(fid, '('//trim(Fmt)//'"- [ ]")', iostat=ErrStat) 
   else
      Fmt = '('//trim(Fmt)//'"- [",'//trim(Num2LStr(nc))//'('//VarFmt//',","),"]")'   
      write(fid, Fmt, iostat=ErrStat) A(:)
   endif
   if (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error writing array '//trim(key)//' to YAML file'
   end if
end subroutine yaml_write_array1I

subroutine yaml_write_array1R4(fid, key, A, VarFmt, ErrStat, ErrMsg, level, comment)
   integer(IntKi),             intent(in   ) :: fid     !< File Unit
   character(len=*),           intent(in   ) :: key     !< Array name
   real(SiKi), dimension(:),   intent(in   ) :: A       !< Array
   character(len=*),           intent(in   ) :: VarFmt  !< Format for printing real numbers  
   integer,                    intent(  out) :: ErrStat !< A non-zero value indicates an error occurred
   character(len=*),           intent(  out) :: ErrMsg  !< Error message if errstat /= errid_none
   integer(IntKi),   optional, intent(in   ) :: level   !< indentation level
   character(len=*), optional, intent(in   ) :: comment !< 
   integer            :: nc  ! size (rows and columns) of A
   character(256)     :: Fmt
   ErrStat = ErrID_None   
   ErrMsg  = ""
   nc = size(A,1)
   Fmt = ''
   if (present(level)) Fmt = trim(Num2LStr(level*INDENT_SPACES))//'X,'
   if (present(comment)) then
      write(fid, '('//trim(Fmt)//'A,": # ",I0," x ",I0,1X,A)', iostat=ErrStat ) trim(key), 1, nc, trim(comment)
   else
      write(fid, '('//trim(Fmt)//'A,": # ",I0," x ",I0)'     , iostat=ErrStat ) trim(key),1,nc
   end if
   if (present(level)) then
      Fmt = trim(Num2LStr((level+1)*INDENT_SPACES))//'X,'
   else
      Fmt = trim(Num2LStr(INDENT_SPACES))//'X,'
   endif
   if (nc==0) then
      write(fid, '('//trim(Fmt)//'"- [ ]")', iostat=ErrStat) 
   else
      Fmt = '('//trim(Fmt)//'"- [",'//trim(Num2LStr(nc))//'('//VarFmt//',","),"]")'   
      write(fid, Fmt, iostat=ErrStat) A(:)
   endif
   if (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error writing array '//trim(key)//' to YAML file'
   end if
end subroutine yaml_write_array1R4

subroutine yaml_write_array1R8(fid, key, A, VarFmt, ErrStat, ErrMsg, level, comment)
   integer(IntKi),             intent(in   ) :: fid     !< File Unit
   character(len=*),           intent(in   ) :: key     !< Array name
   real(R8Ki), dimension(:),   intent(in   ) :: A       !< Array
   character(len=*),           intent(in   ) :: VarFmt  !< Format for printing real numbers  
   integer,                    intent(  out) :: ErrStat !< A non-zero value indicates an error occurred
   character(len=*),           intent(  out) :: ErrMsg  !< Error message if errstat /= errid_none
   integer(IntKi),   optional, intent(in   ) :: level   !< indentation level
   character(len=*), optional, intent(in   ) :: comment !< 
   integer            :: nc  ! size (rows and columns) of A
   character(256)     :: Fmt
   ErrStat = ErrID_None   
   ErrMsg  = ""
   nc = size(A,1)
   Fmt = ''
   if (present(level)) Fmt = trim(Num2LStr(level*INDENT_SPACES))//'X,'
   if (present(comment)) then
      write(fid, '('//trim(Fmt)//'A,": # ",I0," x ",I0,1X,A)', iostat=ErrStat ) trim(key), 1, nc, trim(comment)
   else
      write(fid, '('//trim(Fmt)//'A,": # ",I0," x ",I0)'     , iostat=ErrStat ) trim(key), 1, nc
   end if
   if (present(level)) then
      Fmt = trim(Num2LStr((level+1)*INDENT_SPACES))//'X,'
   else
      Fmt = trim(Num2LStr(INDENT_SPACES))//'X,'
   endif
   if (nc==0) then
      write(fid, '('//trim(Fmt)//'"- [ ]")', iostat=ErrStat) 
   else
      Fmt = '('//trim(Fmt)//'"- [",'//trim(Num2LStr(nc))//'('//VarFmt//',","),"]")'   
      write(fid, Fmt, iostat=ErrStat) A(:)
   endif
   if (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error writing array '//trim(key)//' to YAML file'
   end if
end subroutine yaml_write_array1R8


! --------------------------------------------------------------------------------}
! --- Write Array 2D 
! --------------------------------------------------------------------------------{
subroutine yaml_write_array2I(fid, key, A, VarFmt, ErrStat, ErrMsg, level, comment, label)
   integer(IntKi),             intent(in   ) :: fid     !< File Unit
   character(len=*),           intent(in   ) :: key     !< Array name
   integer(IntKi), dimension(:,:), intent(in   ) :: A   !< Array
   character(len=*),           intent(in   ) :: VarFmt  !< Format for printing real numbers  
   integer,                    intent(  out) :: ErrStat !< A non-zero value indicates an error occurred
   character(len=*),           intent(  out) :: ErrMsg  !< Error message if errstat /= errid_none
   integer(IntKi),   optional, intent(in   ) :: level   !< indentation level
   character(len=*), optional, intent(in   ) :: comment !< 
   logical,          optional, intent(in   ) :: label   !< If present, add a index label at end of line
   integer            :: nr, nc, i  ! size (rows and columns) of A
   character(256)     :: Fmt
   ErrStat = ErrID_None   
   ErrMsg  = ""
   nr = size(A,1)
   nc = size(A,2)
   Fmt=''

   ! Indent and Key
   if (present(level)) Fmt = trim(Num2LStr(level*INDENT_SPACES))//'X,'
   if (present(comment)) then
      write(fid, '('//trim(Fmt)//'A,": # ",I0," x ",I0,1X,A)', iostat=ErrStat ) trim(key), nr, nc, trim(comment)
   else
      write(fid, '('//trim(Fmt)//'A,": # ",I0," x ",I0)'     , iostat=ErrStat ) trim(key),nr,nc
   end if
   if (present(level)) then
      Fmt = trim(Num2LStr((level+1)*INDENT_SPACES))//'X,'
   else
      Fmt = trim(Num2LStr(INDENT_SPACES))//'X,'
   endif
   ! Data
   if (nr==0) then
      write(fid, '('//trim(Fmt)//'"- [ ]")', iostat=ErrStat) 
   else
      ! YAML Line format
      if (present(label)) then
         if (nc==1) then
            Fmt = '('//trim(Fmt)//'"- [", '//VarFmt//', "] # ",I0)'   
         else
            Fmt = '('//trim(Fmt)//'"- [", '//trim(Num2LStr(nc-1))//'('//VarFmt//', ","), '//VarFmt//', "] # ",I0)'   
         endif
      else
         if (nc==1) then
            Fmt = '('//trim(Fmt)//'"- [", '//trim(Num2LStr(nc))//'('//VarFmt//'), "]")'   
         else
            Fmt = '('//trim(Fmt)//'"- [", '//trim(Num2LStr(nc-1))//'('//VarFmt//', ","), '//VarFmt//', "]")'   
         endif
      endif
      ! Write line by line
      do i=1,nr
         if (present(label)) then
            write(fid, Fmt, iostat=ErrStat) A(i,:), i
         else
            write(fid, Fmt, iostat=ErrStat) A(i,:)
         endif
      end do
   endif
   if (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error writing array '//trim(key)//' to YAML file'
   end if
end subroutine yaml_write_array2I

subroutine yaml_write_array2R4(fid, key, A, VarFmt, ErrStat, ErrMsg, level, comment, AllFmt)
   integer(IntKi),             intent(in   ) :: fid     !< File Unit
   character(len=*),           intent(in   ) :: key     !< Array name
   real(SiKi), dimension(:,:), intent(in   ) :: A       !< Array
   character(len=*),           intent(in   ) :: VarFmt  !< Format for printing real numbers  
   integer,                    intent(  out) :: ErrStat !< A non-zero value indicates an error occurred
   character(len=*),           intent(  out) :: ErrMsg  !< Error message if errstat /= errid_none
   integer(IntKi),   optional, intent(in   ) :: level   !< indentation level
   character(len=*), optional, intent(in   ) :: comment !< 
   character(len=*), optional, intent(in   ) :: AllFmt  !< Format for printing a line
   integer            :: nr, nc, i  ! size (rows and columns) of A
   character(256)     :: Fmt
   ErrStat = ErrID_None   
   ErrMsg  = ""
   nr = size(A,1)
   nc = size(A,2)
   Fmt=''

   ! Indent and Key
   if (present(level)) Fmt = trim(Num2LStr(level*INDENT_SPACES))//'X,'
   if (present(comment)) then
      write(fid, '('//trim(Fmt)//'A,": # ",I0," x ",I0,1X,A)', iostat=ErrStat ) trim(key), nr, nc, trim(comment)
   else
      write(fid, '('//trim(Fmt)//'A,": # ",I0," x ",I0)'     , iostat=ErrStat ) trim(key),nr,nc
   end if
   if (present(level)) then
      Fmt = trim(Num2LStr((level+1)*INDENT_SPACES))//'X,'
   else
      Fmt = trim(Num2LStr(INDENT_SPACES))//'X,'
   endif
   ! Data
   if (nr==0) then
      write(fid, '('//trim(Fmt)//'"- [ ]")', iostat=ErrStat) 
   else
      ! YAML Line format
      if (present(AllFmt)) then
         Fmt = '('//trim(Fmt)//'"- [",'//trim(AllFmt)//'"]")'   
      elseif (nc==1) then
         Fmt = '('//trim(Fmt)//'"- [", '//trim(Num2LStr(nc))//'('//VarFmt//'), "]")'   
      else
         Fmt = '('//trim(Fmt)//'"- [", '//trim(Num2LStr(nc-1))//'('//VarFmt//', ","), '//VarFmt//', "]")'   
      endif
      ! Write line by line
      do i=1,nr
         write(fid, Fmt, iostat=ErrStat) A(i,:)
      end do
   endif
   if (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error writing array '//trim(key)//' to YAML file'
   end if
end subroutine yaml_write_array2R4

subroutine yaml_write_array2R8(fid, key, A, VarFmt, ErrStat, ErrMsg, level, comment, AllFmt)
   integer(IntKi),             intent(in   ) :: fid     !< File Unit
   character(len=*),           intent(in   ) :: key     !< Array name
   real(R8Ki), dimension(:,:), intent(in   ) :: A       !< Array
   character(len=*),           intent(in   ) :: VarFmt  !< Format for printing real numbers  
   integer,                    intent(  out) :: ErrStat !< A non-zero value indicates an error occurred
   character(len=*),           intent(  out) :: ErrMsg  !< Error message if errstat /= errid_none
   integer(IntKi),   optional, intent(in   ) :: level   !< indentation level
   character(len=*), optional, intent(in   ) :: comment !< 
   character(len=*), optional, intent(in   ) :: AllFmt  !< Format for printing a line
   integer            :: nr, nc, i  ! size (rows and columns) of A
   character(256)     :: Fmt
   ErrStat = ErrID_None   
   ErrMsg  = ""
   nr = size(A,1)
   nc = size(A,2)
   Fmt=''

   ! Indent and Key
   if (present(level)) Fmt = trim(Num2LStr(level*INDENT_SPACES))//'X,'
   if (present(comment)) then
      write(fid, '('//trim(Fmt)//'A,": # ",I0," x ",I0,1X,A)', iostat=ErrStat ) trim(key), nr, nc, trim(comment)
   else
      write(fid, '('//trim(Fmt)//'A,": # ",I0," x ",I0)'     , iostat=ErrStat ) trim(key),nr,nc
   end if
   if (present(level)) then
      Fmt = trim(Num2LStr((level+1)*INDENT_SPACES))//'X,'
   else
      Fmt = trim(Num2LStr(INDENT_SPACES))//'X,'
   endif
   ! Data
   if (nr==0) then
      write(fid, '('//trim(Fmt)//'"- [ ]")', iostat=ErrStat) 
   else
      ! YAML Line format
      if (present(AllFmt)) then
         Fmt = '('//trim(Fmt)//'"- [",'//trim(AllFmt)//'"]")'   
      elseif (nc==1) then
         Fmt = '('//trim(Fmt)//'"- [", '//trim(Num2LStr(nc))//'('//VarFmt//'), "]")'   
      else
         Fmt = '('//trim(Fmt)//'"- [", '//trim(Num2LStr(nc-1))//'('//VarFmt//', ","), '//VarFmt//', "]")'   
      endif
      ! Write line by line
      do i=1,nr
         write(fid, Fmt, iostat=ErrStat) A(i,:)
      end do
   endif
   if (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error writing array '//trim(key)//' to YAML file'
   end if
end subroutine yaml_write_array2R8

end module YAML
