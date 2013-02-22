module qsort_c_module

implicit none
public :: QsortC
private :: Partition

contains

recursive subroutine QsortC(A)
  integer, intent(in out), dimension(:,:) :: A
  integer :: iq

  if(size(A)/2 > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1, :))
     call QsortC(A(iq:, :))
  endif
end subroutine QsortC

subroutine Partition(A, marker)
  integer, intent(in out), dimension(:,:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  integer :: temp
  integer :: x      ! pivot point
  x = A(1, 1)
  i= 0
  j= size(A)/2 + 1

  do
     j = j-1
     do
        if (A(j,1) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i,1) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i,1)
        A(i,1) = A(j,1)
        A(j,1) = temp

        temp = A(i,2)
        A(i,2) = A(j,2)
        A(j,2) = temp

     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module qsort_c_module