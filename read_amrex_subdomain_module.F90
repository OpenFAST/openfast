module read_amrex_subdomain_module
  use iso_c_binding
  implicit none
  private
  public :: ReadAMReXHeader, ReadAMReXSubdomain

  integer, public, parameter :: IntKi = c_int
  integer, public, parameter :: ReKi = c_double
  integer, public, parameter :: DbKi = c_double
  integer, public, parameter :: SiKi = c_double

  interface
     subroutine read_amrex_header (name,dims,origin,dx,t) bind(c)
       import
       implicit none
       character(kind=c_char), intent(in) :: name(*)
       integer(kind=c_int), intent(inout) :: dims(3)
       real(kind=c_double), intent(inout) :: origin(3)
       real(kind=c_double), intent(inout) :: dx(3)
       real(kind=c_double), intent(inout) :: t
     end subroutine read_amrex_header
  end interface

contains

  ! This subroutine is not thread safe! Do NOT call this in omp parallel region!!!
  ! If we want, we could add an assertion.
  subroutine ReadAMReXHeader (FileName, dims, origin, gridSpacing, time)
    character(1024), intent(in   )  :: FileName        ! Name of Header file 
    integer(IntKi),  intent(  out)  :: dims(3)         ! Dimension of the full grid (nx,ny,nz)
    real(ReKi),      intent(  out)  :: origin(3)       ! The lower-left corner of the full grid (x0,y0,z0)
    real(ReKi),      intent(  out)  :: gridSpacing(3)  ! Spacing between grid points (dx,dy,dz)
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

  subroutine ReadAMReXSubdomain
  end subroutine ReadAMReXSubdomain

end module read_amrex_subdomain_module
