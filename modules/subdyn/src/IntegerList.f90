!> Module providing suport for an integer list stored as an array and  not a chained list
!! Used since registry does not support pointer with recursive types.
module IntegerList
   use SubDyn_Types, only: IList
   use NWTC_Library, only: IntKi, ReKi, AllocAry, ErrID_None, ErrID_Fatal, num2lstr

   implicit none

   public :: IList

   public :: init_list
   public :: destroy_list
   public :: len
   public :: append
   public :: pop
   public :: get
   public :: find
   public :: sort
   public :: reverse
   interface pop
      module procedure pop_last
      module procedure pop_item
   end interface
   interface init_list
      module procedure init_list_n_def
      module procedure init_list_vect
   end interface
   interface find
      module procedure find_list
      module procedure find_intarray
   end interface
   interface unique
      module procedure unique_list
      module procedure unique_intarray
      module procedure unique_intarray_in_place
   end interface
contains

   !> Concatenate lists: I3=[I1,I2]
   subroutine concatenate_lists(I1,I2,I3, ErrStat, ErrMsg)
      integer(intki), intent(in)    :: i1(:), i2(:)
      integer(intki), intent(out)   :: i3(:)
      integer(IntKi), intent(  out) :: ErrStat   !< Error status of the operation
      character(*),   intent(  out) :: ErrMsg    !< Error message if ErrStat /    = ErrID_None
      ErrStat=ErrID_None
      ErrMsg=''
      I3(1:size(I1)) = I1 
      I3(size(I1)+1:size(I1)+size(I2)) = I2
   endsubroutine
   subroutine concatenate_3lists(I1,I2,I3,I4, ErrStat, ErrMsg)
      integer(intki), intent(in)    :: i1(:), i2(:), i3(:)
      integer(intki), intent(out)   :: i4(:)
      integer(IntKi), intent(  out) :: ErrStat   !< Error status of the operation
      character(*),   intent(  out) :: ErrMsg    !< Error message if ErrStat /    = ErrID_None
      ErrStat=ErrID_None
      ErrMsg=''
      I4(                  1:size(I1)                  ) = I1 
      I4(size(I1)         +1:size(I1)+size(I2)         ) = I2
      I4(size(I1)+size(I2)+1:size(I1)+size(I2)+size(I3)) = I3
   endsubroutine

   !> Set difference: I3=I1-I2 (assumes I1 is biggger than I2), elements of I1 not in I2
   subroutine lists_difference(I1, I2, I3, ErrStat, ErrMsg)
      integer(IntKi), intent(in)    :: I1(:), I2(:)
      integer(IntKi), intent(out)   :: I3(:)
      integer(IntKi), intent(  out) :: ErrStat      !< Error status of the operation
      character(*),   intent(  out) :: ErrMsg       !< Error message if ErrStat /= ErrID_None
      integer(IntKi) :: I
      logical, dimension(:), allocatable :: bUnique
      ErrStat = ErrID_None
      ErrMsg  = ""
      allocate(bUnique(1:size(I1)))
      ! Then, remove DOFs on the boundaries:
      DO i = 1, size(I1)  !Boundary DOFs (Interface + Constraints)
         if (find(I2,I1(i))>0) then
            bUnique(I) = .false.
         else
            bUnique(I) = .true.
         endif
      ENDDO
      if (count(bUnique) /= size(I3)) then
         ErrStat=ErrID_Fatal; ErrMsg='Storage for list difference is of wrong size'; return
      endif
      I3 = pack(I1, bUnique)
      deallocate(bUnique)
   endsubroutine

   !> Initialize an integer list
   subroutine init_list_n_def(L,n,default_val,ErrStat,ErrMsg)
      type(IList), intent(inout)                  :: L !< List
      integer(IntKi), intent(in)                  :: n !< number of initial values
      integer(IntKi), intent(in)                  :: default_val !< default values
      integer(IntKi),               intent(  out) :: ErrStat     !< Error status of the operation
      character(*),                 intent(  out) :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      ErrStat = ErrID_None
      ErrMsg  = ""
      call AllocAry(L%List, n, 'L%List', ErrStat, ErrMsg) 
      if (ErrStat/=ErrID_None) return
      L%List(1:n) = default_val
   end subroutine init_list_n_def

   subroutine init_list_vect(L,vect,ErrStat,ErrMsg)
      type(IList), intent(inout)               :: L    !< List
      integer(IntKi), dimension(:), intent(in) :: vect !< number of initial values
      integer(IntKi),               intent(  out) :: ErrStat     !< Error status of the operation
      character(*),                 intent(  out) :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      ErrStat = ErrID_None
      ErrMsg  = ""
      call AllocAry(L%List, size(vect), 'L%List', ErrStat, ErrMsg) 
      if (ErrStat/=ErrID_None) return
      L%List = vect
   end subroutine init_list_vect

   !> Deallocate list
   subroutine destroy_list(L,ErrStat,ErrMsg)
      type(IList), intent(inout)                  :: L !< List
      integer(IntKi),               intent(  out) :: ErrStat     !< Error status of the operation
      character(*),                 intent(  out) :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      ErrStat = ErrID_None
      ErrMsg  = ""
      if (allocated(L%List)) deallocate(L%List)
   end subroutine destroy_list

   !> Returns list length
   integer function len(L)
      type(IList), intent(in) :: L
      if (allocated(L%List)) then
         len=size(L%List)
      else
         len=0
      endif
   end function len

   !> Append element to list
   subroutine append(L,e, ErrStat, ErrMsg)
      type(IList),     intent(inout) :: L
      integer(IntKi),  intent(in)    :: e
      integer(IntKi),  intent(  out) :: ErrStat !< Error status of the operation
      character(*),    intent(  out) :: ErrMsg  !< Error message if ErrStat /    = ErrID_None
      ErrStat = ErrID_None
      ErrMsg  = ""
      if (allocated(L%List)) then
         call resize_array(L%List,len(L)+1,e)
      else
         call init_list(L, 1, e, ErrStat, ErrMsg)
      endif
   end subroutine append

   !> Get element i from list
   integer function get(L,i, ErrStat, ErrMsg)
      type(IList), intent(inout)     :: L
      integer(IntKi), intent(in)     :: i
      integer(IntKi),  intent(  out) :: ErrStat !< Error status of the operation
      character(*),    intent(  out) :: ErrMsg  !< Error message if ErrStat /    = ErrID_None
      if ((i<=0).or.(i>len(L))) then 
         ErrStat=ErrID_Fatal
         ErrMsg="Index out of bound "//trim(num2lstr(i))//", list length is "//trim(num2lstr(len(L)))
         get=-9999
      else
         ErrStat = ErrID_None
         ErrMsg  = ""
         get = L%List(i) ! No error handling, throws "index array out of bound", like a regular array
      endif
   end function get

   !> Pop last element of the list and reduce list size by 1
   integer function pop_last(L,ErrStat,ErrMsg)
      type(IList), intent(inout)     :: L
      integer(IntKi),  intent(  out) :: ErrStat !< Error status of the operation
      character(*),    intent(  out) :: ErrMsg  !< Error message if ErrStat /    = ErrID_None
      integer(IntKi) :: n
      ErrStat = ErrID_None
      ErrMsg  = ""
      n=len(L)
      pop_last = get(L, n, ErrStat, ErrMsg) ! index array out of bound will be thrown
      call resize_array(L%List,n-1,0)
   end function pop_last

   !> Pop element i from the list and reduce the size of the list by 1
   integer function pop_item(L,i,ErrStat,ErrMsg)
      type(IList), intent(inout)     :: L
      integer(IntKi), intent(in)     :: i
      integer(IntKi),  intent(  out) :: ErrStat !< Error status of the operation
      character(*),    intent(  out) :: ErrMsg  !< Error message if ErrStat /    = ErrID_None
      integer(IntKi) :: n
      ErrStat = ErrID_None
      ErrMsg  = ""
      n=len(L)
      pop_item = get(L, i, ErrStat, ErrMsg) ! index array out of bound will be thrown
      L%List(i:n-1)=L%List(i+1:n)
      call resize_array(L%List,n-1,0)
   end function pop_item

   !> Sort list
   subroutine sort(L, ErrStat, ErrMsg)
      type(IList),     intent(inout) :: L
      integer(IntKi),  intent(  out) :: ErrStat !< Error status of the operation
      character(*),    intent(  out) :: ErrMsg  !< Error message if ErrStat /    = ErrID_None
      ErrStat = ErrID_None
      ErrMsg  = ""
      if (allocated(L%List)) then
         call sort_in_place(L%List)
      else
         ErrStat=ErrID_Fatal
         ErrMsg="Cannot sort a list not allocated"
      endif
   end subroutine sort

   !> Reverse list
   subroutine reverse(L, ErrStat, ErrMsg)
      type(IList),     intent(inout) :: L
      integer(IntKi),  intent(  out) :: ErrStat !< Error status of the operation
      character(*),    intent(  out) :: ErrMsg  !< Error message if ErrStat /    = ErrID_None
      integer(IntKi) :: i
      integer(IntKi) :: n
      ErrStat = ErrID_None
      ErrMsg  = ""
      n=len(L)
      do i =1,int(n/2)
         call swap(i, n-i+1)
      enddo
      contains 
         subroutine swap(i,j)
            integer(IntKi), intent(in) :: i,j
            integer(IntKi) :: tmp
            tmp=L%List(i)
            L%List(i) = L%List(j)
            L%List(j) = tmp
         end subroutine
   end subroutine reverse


   !> Returns index of element e in L, returns 0 if not found
   !! NOTE: list but be sorted to call this function
   integer(IntKi) function find_list(L, e, ErrStat, ErrMsg)
      type(IList),     intent(inout) :: L  
      integer(IntKi),  intent(in   ) :: e  
      integer(IntKi),  intent(  out) :: ErrStat !< Error status of the operation
      character(*),    intent(  out) :: ErrMsg  !< Error message if ErrStat /    = ErrID_None
      ErrStat = ErrID_None
      ErrMsg  = ""
      if (len(L)>0) then
         find_list = binary_search(L%List, e) ! Binary search returns index for inequality List(i)<=e
         if (find_list>0) then
            if (L%List(find_list)/=e) then
               find_list=-1
            endif
         endif
      else
         find_list=-1
      endif
   end function find_list

   !> Unique, in place
   subroutine unique_list(L, ErrStat, ErrMsg)
      type(IList),     intent(inout) :: L  
      integer(IntKi),  intent(  out) :: ErrStat !< Error status of the operation
      character(*),    intent(  out) :: ErrMsg  !< Error message if ErrStat /    = ErrID_None
      ErrStat = ErrID_None
      ErrMsg  = ""
      if (len(L)>0) then
         call unique_intarray_in_place(L%List)
      endif
   end subroutine

   !> Print
   subroutine print_list(L, varname, u_opt)
      type(IList),     intent(in)          :: L
      character(len=*),intent(in)          :: varname
      integer(IntKi),  intent(in),optional :: u_opt
      ! 
      character(len=*),parameter :: IFMT='I7.0' !<
      integer(IntKi) :: u
      integer(IntKi) :: n
      character(len=20) :: fmt
      ! Optional args
      if (present(u_opt)) then
         u=u_opt
      else
         u=6
      endif
      n=len(L)
      if (n>0) then
         write(fmt,*) n
         write(u,"(A,A,"// adjustl(fmt)//IFMT//",A)") varname,"=[",L%List,"];"
      else
         write(u,'(A,A)') varname,'=[];'
      endif
   end subroutine print_list

   ! --------------------------------------------------------------------------------
   ! --- Generic helper functions (should be part of NWTC library)
   ! --------------------------------------------------------------------------------
   !> Sort integer array in place 
   pure subroutine sort_in_place(a)
      integer(IntKi), intent(inout), dimension(:) :: a
      integer(IntKi) :: temp
      integer(IntKi) :: i, j
      do i = 2, size(a)
         j = i - 1
         temp = a(i)
         do while (j>=1 .and. a(j)>temp)
            a(j+1) = a(j)
            j = j - 1
            if (j<1) then
               exit
            endif
         end do
         a(j+1) = temp
      end do
   end subroutine sort_in_place

    !> Performs binary search and return the largest index such that x(i) <= x0 
    !! allows equlity 
    Integer(IntKi) function binary_search(x, x0) result(i_inf)
        ! Arguments declarations 
        integer(IntKi), dimension(:),intent(in) :: x  !<  x *sorted* vector 
        integer(IntKi), intent(in) :: x0  !<  
        ! Variable declarations 
        integer(IntKi) :: i_sup  !<  
        integer(IntKi) :: mid  !<  
        i_inf=1
        i_sup=size(x)
        ! Safety test 
        if (x0<x(1)) then 
            i_inf=-1
            return
        end if 
        if (x0>=x(i_sup)) then 
            i_inf=i_sup
            return
        end if 
        ! We loop until we narrow down to one index 
        do while (i_inf+1<i_sup) 
            mid=(int((i_inf+i_sup)/2))
            if (x(mid)<=x0) then 
                i_inf=mid
            else
                i_sup=mid
            end if 
        end do 
    end function binary_search 

   !> Returns index of val in Array (val is an integer!)
   ! NOTE: in the future use intrinsinc function findloc
   function find_intarray(Array, Val) result(i)
      integer(IntKi), dimension(:), intent(in) :: Array !< Array to search in
      integer(IntKi), intent(in)               :: val   !< Val
      integer(IntKi)                           :: i     !< Index of joint in joint table
      i = 1
      do while ( i <= size(Array) )
         if ( Val == Array(i) ) THEN
            return ! Exit when found
         else
            i = i + 1
         endif
      enddo
      i=-1
   end function

    !> return in res the unique values of v
    subroutine unique_intarray(v,res) 
       ! Arguments
       integer(IntKi),dimension(:),intent(in) :: v
       integer(IntKi),dimension(:),allocatable::res
       !
       integer(IntKi),dimension(:),pointer::tmp
       integer :: k    !< number of unique elements
       integer :: i, j
       if (allocated(res)) deallocate(res)
       allocate(tmp(1:size(v)))
       k = 1
       tmp(1) = v(1)
       outer: do i=2,size(v)
          do j=1,k
             if (tmp(j) == v(i)) then
                ! Found a match so start looking again
                cycle outer
             end if
          end do
          ! No match found so add it to the output
          k = k + 1
          tmp(k) = v(i)
       end do outer
       allocate(res(1:k))
       res(1:k)=tmp(1:k)
       deallocate(tmp)
    end subroutine

    subroutine unique_intarray_in_place(v) 
        integer(IntKi),dimension(:),allocatable :: v
        integer(IntKi),dimension(:),allocatable::res
        integer :: k    !< number of unique elements
        integer :: i, j
        allocate(res(1:size(v)))
        k = 1
        res(1) = v(1)
        outer: do i=2,size(v)
            do j=1,k
                if (res(j) == v(i)) then
                    ! Found a match so start looking again
                    cycle outer
                end if
            end do
            ! No match found so add it to the output
            k = k + 1
            res(k) = v(i)
        end do outer
        deallocate(v)
        allocate(v(1:k))
        v(1:k)=res(1:k)
        deallocate(res)
    end subroutine

   !> Resize integer array of dimension 1
   subroutine resize_array(array,nNewSize,default_val)
      integer(IntKi),dimension(:),allocatable,intent(inout) :: array
      integer(IntKi) , intent(in) :: nNewSize
      integer(IntKi), intent(in) :: default_val
      ! Local variables
      integer(IntKi),dimension(:),allocatable :: tmp !< backup of input
      integer(IntKi) :: nDimTmp
      integer(IntKi) :: AllocateStatus
      ! To save memory, if nNewSize is below second dim, we take the min
      nDimTmp= min(size(array,1),nNewSize)
      ! Making of copy of the input
      allocate(tmp(1:nDimTmp), STAT = AllocateStatus)
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      tmp(1:nDimTmp)=array(1:nDimTmp)
      ! Reallocating the array 
      deallocate(array)
      allocate(array(1:nNewSize), STAT = AllocateStatus)
      if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      ! We copy the original data into it
      array(1:nDimTmp)=tmp(1:nDimTmp)
      if(nDimTmp+1<=nNewSize) array(nDimTmp+1:nNewSize)=default_val
   end subroutine 

   !> Append two integer arrays of dimension 1
   subroutine append_arrays(array1,n1,array2,n2)
      integer(IntKi), dimension(:), allocatable :: array1
      integer(IntKi), dimension(:)              :: array2
      integer(IntKi), intent(inout)             :: n1     !< SIDE EFFECTS
      integer(IntKi), intent(in)                :: n2
      ! Local variables
      integer  :: nNew
      nNew=n1+n2
      ! --- Making enough space if needed
      if(nNew>size(array1,1)) then
         call resize_array(array1,nNew,0)
      endif
      ! --- Appending
      array1((n1+1):(n1+n2))=array2(1:n2)
      ! updating n1
      n1=n1+n2;
   end subroutine 

end module IntegerList
