module read_amrex_subdomain_module
  use iso_c_binding
  implicit none
  private
  public :: ReadAMReXHeader, ReadAMReXSubdomain

  integer(c_int), parameter :: error_stat_noerror = 0
  integer(c_int), parameter :: error_stat_severe = 3  ! Points outside domain.
  integer(c_int), parameter :: error_stat_fatal = 4   ! fatal error

  interface
     subroutine set_error_code (noerror, severe, fatal) bind(c)
       import
       implicit none
       integer(kind=c_int), value, intent(in) :: noerror, severe, fatal
     end subroutine set_error_code

     subroutine read_amrex_header (name,dims,origin,dx,t,err,msg,nmsg) bind(c)
       import
       implicit none
       character(kind=c_char), intent(in) :: name(*)
       integer(kind=c_int), intent(out) :: dims(3)
       real(kind=c_double), intent(out) :: origin(3)
       real(kind=c_double), intent(out) :: dx(3)
       real(kind=c_double), intent(out) :: t
       integer(kind=c_int), intent(out) :: err
       character(kind=c_char), intent(out) :: msg(*)
       integer(kind=c_int), value, intent(in) :: nmsg
     end subroutine read_amrex_header

     subroutine read_amrex_subdomain (a,lo,hi,err,msg,nmsg) bind(c)
       import
       implicit none
       integer(kind=c_int), intent(in) :: lo(3), hi(3)
       real(kind=c_double), intent(out) :: a(*)
       integer(kind=c_int), intent(out) :: err
       character(kind=c_char), intent(out) :: msg(*)
       integer(kind=c_int), value, intent(in) :: nmsg
     end subroutine read_amrex_subdomain
  end interface

contains

  function string_f_to_c (fstr) result(cstr)
    character(*), intent(in) :: fstr
    character(kind=c_char) :: cstr(len_trim(fstr)+1)
    integer :: i, n
    n = len_trim(fstr)
    do i = 1, n
       cstr(i) = fstr(i:i)
    end do
    cstr(n+1) = c_null_char
  end function string_f_to_c

  function string_c_to_f (cstr) result(fstr)
    character(kind=c_char), intent(in) :: cstr(:)
    character(len=size(cstr)-1) :: fstr
    integer :: i, n
    n = size(cstr)-1   ! skip the null character
    fstr = ""
    do i = 1, n
       if (cstr(i) == c_null_char) exit
       fstr(i:i) = cstr(i)
    enddo
  end function string_c_to_f

  ! This subroutine is not thread safe! Do NOT call this in omp parallel region!!!
  ! If we want, we could add an assertion.
  subroutine ReadAMReXHeader (FileName, dims, origin, gridSpacing, time, ErrorStat, ErrorMsg)
    character(1024), intent(in   )  :: FileName        ! Name of Header file
    integer(c_int),  intent(  out)  :: dims(3)         ! Dimension of the subdomain grid (nx,ny,nz)
    real(c_double),  intent(  out)  :: origin(3)       ! The lower-left corner of the subdomain
    real(c_double),  intent(  out)  :: gridSpacing(3)  ! Spacing between grid points (dx,dy,dz)
    real(c_double),  intent(  out)  :: time
    integer(c_int),  intent(  out)  :: ErrorStat
    character(1024), intent(  out)  :: ErrorMsg
    character(kind=c_char) :: error_msg(len(ErrorMsg)+1)
    call set_error_code(error_stat_noerror, error_stat_severe, error_stat_fatal)
    call read_amrex_header(string_f_to_c(FileName), dims, origin, gridSpacing, time, &
         &                 ErrorStat, error_msg, size(error_msg,kind=c_int))
    ErrorMsg = string_c_to_f(error_msg)
  end subroutine ReadAMReXHeader

  ! The natural ordering in amrex is (nx,ny,nz,3).  We will transpose it to
  ! (3,nx,ny,nz).  The starting index of the entire subdomain data is 0.
  ! This subroutine reads part of the subdomain specified by lo and hi.  The
  ! center of cell (i,j,k) is at origin + ((i,j,k)+0.5)*gridSpacing.
  ! ErrorStat is set to error_stat_severe, if some part of the array cannot
  ! be filled. This subroutine is thread safe. But using too many threads
  ! may not be good for performance.
  subroutine ReadAMReXSubdomain (a, lo, hi, ErrorStat, ErrorMsg)
    integer(c_int), intent(in) :: lo(3), hi(3)
    real(c_double), intent(out) :: a(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    integer(c_int),  intent(  out)  :: ErrorStat
    character(1024), intent(  out)  :: ErrorMsg
    character(kind=c_char) :: error_msg(len(ErrorMsg)+1)
    call read_amrex_subdomain(a, lo, hi, ErrorStat, error_msg, size(error_msg,kind=c_int))
    ErrorMsg = string_c_to_f(error_msg)
  end subroutine ReadAMReXSubdomain

end module read_amrex_subdomain_module
