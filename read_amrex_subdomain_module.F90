module read_amrex_subdomain_module
  use iso_c_binding
  implicit none
  private
  public :: ReadAMReXHeader, ReadAMReXSubdomain

  interface
     subroutine read_amrex_header (name,dims,origin,dx,t) bind(c)
       import
       implicit none
       character(kind=c_char), intent(in) :: name(*)
       integer(kind=c_int), intent(out) :: dims(3)
       real(kind=c_double), intent(out) :: origin(3)
       real(kind=c_double), intent(out) :: dx(3)
       real(kind=c_double), intent(out) :: t
     end subroutine read_amrex_header

     subroutine read_amrex_subdomain (a,lo,hi) bind(c)
       import
       implicit none
       integer(kind=c_int), intent(in) :: lo(3), hi(3)
       real(kind=c_double), intent(out) :: a(*)
     end subroutine read_amrex_subdomain
  end interface

contains

  ! This subroutine is not thread safe! Do NOT call this in omp parallel region!!!
  ! If we want, we could add an assertion.
  subroutine ReadAMReXHeader (FileName, dims, origin, gridSpacing, time)
    character(1024), intent(in   )  :: FileName        ! Name of Header file 
    integer(c_int),  intent(  out)  :: dims(3)         ! Dimension of the full grid (nx,ny,nz)
    real(c_double),  intent(  out)  :: origin(3)       ! The lower-left corner of the full grid (x0,y0,z0)
    real(c_double),  intent(  out)  :: gridSpacing(3)  ! Spacing between grid points (dx,dy,dz)
    real(c_double),  intent(  out)  :: time

    character(kind=c_char) :: cstr(len_trim(FileName)+1)
    integer :: i, n
    n = len_trim(FileName)
    do i = 1, n
       cstr(i) = FileName(i:i)
    end do
    cstr(n+1) = c_null_char
    call read_amrex_header(cstr, dims, origin, gridSpacing, time)
  end subroutine ReadAMReXHeader

  ! The natural ordering in amrex is (nx,ny,nz,3).
  ! We will transpose it to (3,nx,ny,nz)
  ! The center of cell (i,j,k) is at origins + ((i,j,k)+0.5)*gridSpacing.
  subroutine ReadAMReXSubdomain (a, lo, hi)
    integer(c_int), intent(in) :: lo(3), hi(3)
    real(c_double), intent(out) :: a(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    call read_amrex_subdomain(a, lo, hi)
  end subroutine ReadAMReXSubdomain

end module read_amrex_subdomain_module
