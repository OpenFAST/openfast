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
module JSON
   use Precision, only: IntKi, SiKi, R8Ki
   use NWTC_Base, only: ErrID_None, ErrID_Fatal
   use NWTC_IO, only: num2lstr

   implicit none

   !> Write 2D array to file in JSON format
   interface json_write_array
      module procedure json_write_array2R4  ! Two dimension array of SiKi
      module procedure json_write_array2I   ! Two dimension array of IntKi
      module procedure json_write_array2R8  ! Two dimension array of R8Ki
   end interface

   private

   public :: json_write_array

contains


! --------------------------------------------------------------------------------}
! --- Write Array 2D 
! --------------------------------------------------------------------------------{
subroutine json_write_array2I(fid, key, A, VarFmt, ErrStat, ErrMsg)
   integer(IntKi),             intent(in   ) :: fid     !< File Unit
   character(len=*),           intent(in   ) :: key     !< Array name
   integer(IntKi), dimension(:,:), intent(in   ) :: A   !< Array
   character(len=*),           intent(in   ) :: VarFmt  !< Format for printing real numbers  
   integer,                    intent(  out) :: ErrStat !< A non-zero value indicates an error occurred
   character(len=*),           intent(  out) :: ErrMsg  !< Error message if errstat /= errid_none
   integer            :: nr, nc, i  ! size (rows and columns) of A
   character(256)     :: Fmt
   ErrStat = ErrID_None   
   ErrMsg  = ""
   nr = size(A,1)
   nc = size(A,2)

   ! Key
   write(fid, '(A,": [")', iostat=ErrStat, advance='no') trim(key)
   if (nr==0) then
      write(fid, '("[]]")', iostat=ErrStat) 
   else
      ! JSON Line format
      if (nc==1) then
         Fmt = '("[",'//VarFmt//',"],")'   
      else
         Fmt = '("[",'//trim(Num2LStr(nc-1))//'('//VarFmt//',","),'//VarFmt//',"]")'   
      endif
      ! Write line by line
      do i=1,nr
         write(fid, Fmt, iostat=ErrStat, advance='no') A(i,:)
         if (i<nr) then
            write(fid, '(A)', iostat=ErrStat, advance='no') ','
         endif
      enddo
      write(fid, '("]")', iostat=ErrStat, advance='no') 
   endif
   if (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error writing array '//trim(key)//' to JSON file'
   end if
end subroutine json_write_array2I

subroutine json_write_array2R4(fid, key, A, VarFmt, ErrStat, ErrMsg, AllFmt)
   integer(IntKi),             intent(in   ) :: fid     !< File Unit
   character(len=*),           intent(in   ) :: key     !< Array name
   real(SiKi), dimension(:,:), intent(in   ) :: A       !< Array
   character(len=*),           intent(in   ) :: VarFmt  !< Format for printing real numbers  
   integer,                    intent(  out) :: ErrStat !< A non-zero value indicates an error occurred
   character(len=*),           intent(  out) :: ErrMsg  !< Error message if errstat /= errid_none
   character(len=*), optional, intent(in   ) :: AllFmt  !< Format for printing a line
   integer            :: nr, nc, i  ! size (rows and columns) of A
   character(256)     :: Fmt
   ErrStat = ErrID_None   
   ErrMsg  = ""
   nr = size(A,1)
   nc = size(A,2)

   ! Key
   write(fid, '(A,": [")', iostat=ErrStat, advance='no') trim(key)
   if (nr==0) then
      write(fid, '("[]]")', iostat=ErrStat) 
   else
      ! JSON Line format
      if (present(AllFmt)) then
         Fmt = '("[",'//trim(AllFmt)//'"]")'   
      elseif (nc==1) then
         Fmt = '("[",'//VarFmt//',"],")'   
      else
         Fmt = '("[",'//trim(Num2LStr(nc-1))//'('//VarFmt//',","),'//VarFmt//',"]")'   
      endif
      ! Write line by line
      do i=1,nr
         write(fid, Fmt, iostat=ErrStat, advance='no') A(i,:)
         if (i<nr) then
            write(fid, '(A)', iostat=ErrStat, advance='no') ','
         endif
      enddo
      write(fid, '("]")', iostat=ErrStat, advance='no') 
   endif
   if (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error writing array '//trim(key)//' to JSON file'
   end if
end subroutine json_write_array2R4

subroutine json_write_array2R8(fid, key, A, VarFmt, ErrStat, ErrMsg, AllFmt)
   integer(IntKi),             intent(in   ) :: fid     !< File Unit
   character(len=*),           intent(in   ) :: key     !< Array name
   real(R8Ki), dimension(:,:), intent(in   ) :: A       !< Array
   character(len=*),           intent(in   ) :: VarFmt  !< Format for printing real numbers  
   integer,                    intent(  out) :: ErrStat !< A non-zero value indicates an error occurred
   character(len=*),           intent(  out) :: ErrMsg  !< Error message if errstat /= errid_none
   character(len=*), optional, intent(in   ) :: AllFmt  !< Format for printing a line
   integer            :: nr, nc, i  ! size (rows and columns) of A
   character(256)     :: Fmt
   ErrStat = ErrID_None   
   ErrMsg  = ""
   nr = size(A,1)
   nc = size(A,2)

   ! Key
   write(fid, '(A,": [")', iostat=ErrStat, advance='no') trim(key)
   if (nr==0) then
      write(fid, '("[]]")', iostat=ErrStat) 
   else
      ! YSON Line format
      if (present(AllFmt)) then
         Fmt = '("[",'//trim(AllFmt)//'"]")'   
      elseif (nc==1) then
         Fmt = '("[",'//VarFmt//',"],")'   
      else
         Fmt = '("[",'//trim(Num2LStr(nc-1))//'('//VarFmt//',","),'//VarFmt//',"]")'   
      endif
      ! Write line by line
      do i=1,nr
         write(fid, Fmt, iostat=ErrStat, advance='no') A(i,:)
         if (i<nr) then
            write(fid, '(A)', iostat=ErrStat, advance='no') ','
         endif
      enddo
      write(fid, '("]")', iostat=ErrStat, advance='no') 
   endif
   if (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error writing array '//trim(key)//' to JSON file'
   end if
end subroutine json_write_array2R8

end module JSON
