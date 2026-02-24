module amrex_utils
use iso_c_binding
use NWTC_Library
use amrex_base_module, only: amrex_init, amrex_finalize

implicit none

private

interface
   subroutine amrex_read_header_c(dir_path, t, dims, dx, origin, err_stat, err_msg, err_msg_len) bind(c)
      use iso_c_binding
      import
      implicit none
      character(kind=c_char), intent(in)    :: dir_path(*)
      real(kind=c_double), intent(out)      :: t
      integer(kind=c_int), intent(out)      :: dims(3)
      real(kind=c_double), intent(out)      :: dx(3)
      real(kind=c_double), intent(out)      :: origin(3)
      integer(kind=c_int), intent(out)      :: err_stat
      character(kind=c_char), intent(out)   :: err_msg(*)
      integer(kind=c_int), intent(in)       :: err_msg_len
   end subroutine

   subroutine amrex_read_data_c(dir_path, data, err_stat, err_msg, err_msg_len) bind(c)
      use iso_c_binding
      import
      implicit none
      character(kind=c_char), intent(in)    :: dir_path(*)
      real(kind=c_float), intent(out)       :: data(*)
      integer(kind=c_int), intent(out)      :: err_stat
      character(kind=c_char), intent(out)   :: err_msg(*)
      integer(kind=c_int), intent(in)       :: err_msg_len
   end subroutine
end interface

public :: amrex_init, amrex_finalize
public :: amrex_read_header, amrex_read_data

contains

subroutine amrex_read_header(DirPath, time, nXYZ, dXYZ, oXYZ, ErrStat, ErrMsg)
   character(*), intent(in)      :: DirPath
   real(DbKi), intent(out)       :: Time
   Integer(IntKi), intent(out)   :: nXYZ(3)
   real(ReKi), intent(out)       :: dXYZ(3)
   real(ReKi), intent(out)       :: oXYZ(3)
   integer(IntKi), intent(out)   :: ErrStat
   character(*), intent(out)     :: ErrMsg

   character(c_char), allocatable :: dir_path(:)
   integer(c_int)    :: err_stat_c
   character(c_char) :: err_msg_c(ErrMsgLen)
   integer(c_int)    :: dims(3)
   real(c_double)    :: t, origin(3), gridSpacing(3)
   integer(IntKi)    :: i

   ! Convert directory path to C type
   allocate(dir_path(len_trim(DirPath) + 1))
   dir_path = transfer(DirPath, dir_path) // c_null_char

   ! Call C++ function to read header
   call amrex_read_header_c(dir_path, t, dims, gridSpacing, origin, err_stat_c, err_msg_c, ErrMsgLen)

   ! Transfer outputs back to fortran types
   time = transfer(t, time)
   nXYZ = transfer(dims, nXYZ)
   dXYZ = transfer(gridSpacing, dXYZ)
   oXYZ = transfer(origin, oXYZ)
   ErrStat = transfer(err_stat_c, ErrStat)
   ErrMsg = transfer(err_msg_c, ErrMsg)
   i = index(ErrMsg, c_null_char)
   if (i > 0) ErrMsg = ErrMsg(:i-1)
   
end subroutine

subroutine amrex_read_data(DirPath, gridData, ErrStat, ErrMsg)
   character(*), intent(in)      :: DirPath
   real(SiKi), intent(out)       :: gridData(:,:,:,:)
   integer(IntKi), intent(out)   :: ErrStat
   character(*), intent(out)     :: ErrMsg

   character(c_char), allocatable :: dir_path(:)
   integer(c_int)    :: err_stat_c
   character(c_char) :: err_msg_c(ErrMsgLen)
   integer(IntKi)    :: i

   ! Convert directory path to C type
   allocate(dir_path(len_trim(DirPath) + 1))
   dir_path = transfer(DirPath, dir_path) // c_null_char

   ! Call C++ function to read header
   call amrex_read_data_c(dir_path, gridData, err_stat_c, err_msg_c, ErrMsgLen)

   ! Transfer outputs back to fortran types
   ErrStat = transfer(err_stat_c, ErrStat)
   ErrMsg = transfer(err_msg_c, ErrMsg)
   i = index(ErrMsg, c_null_char)
   if (i > 0) ErrMsg = ErrMsg(:i-1)

end subroutine

end module
