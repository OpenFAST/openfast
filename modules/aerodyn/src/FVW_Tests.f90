module FVW_Tests 

   use NWTC_Library

   use FVW_Types
   use FVW_Subs
   use FVW_VortexTools
   use FVW_Wings
   use FVW_IO
   use FVW_BiotSavart
   use VTK, only : VTK_Misc

   implicit none

   public :: FVW_RunTests
   private

   interface test_equal; module procedure &
         test_equal_i1, &
         test_equal_i0
   end interface
    interface test_almost_equal; module procedure &
          test_almost_equal_0, &
          test_almost_equal_1, &
          test_almost_equal_2
    end interface
contains
   ! --------------------------------------------------------------------------------
   ! --- Helper functions (should be part of NWTC library)
   ! --------------------------------------------------------------------------------
    subroutine test_success(testname,info,bPrint_in)
        character(len=*), intent(in) :: testname
        character(len=*), intent(in) :: info
        logical, intent(in), optional  ::  bPrint_in
        if(present(bPrint_in)) then
            if(bPrint_in) then
                write(*,'(A)')'[ OK ] '//trim(testname)//': '//trim(Info)
            endif
        else
            write(*,'(A)')'[ OK ] '//trim(testname)//': '//trim(Info)
        endif
    end subroutine

    subroutine test_fail(testname,info,bPrint_in,bStop_in)
        character(len=*), intent(in) :: testname
        character(len=*), intent(in) :: info
        logical, intent(in), optional  ::  bPrint_in
        logical, intent(in), optional  ::  bStop_in
        if(present(bPrint_in)) then
            if(bPrint_in) then
                write(*,'(A)')'[FAIL] '//trim(testname)//': '//trim(Info)
            endif
        else
            write(*,'(A)')'[FAIL] '//trim(testname)//': '//trim(Info)
        endif
        if(present(bStop_in)) then
            if(bStop_in) then
                STOP 
            endif
        else
            STOP
        endif
    end subroutine

    subroutine  test_equal_i0(testname,Var,iTry,iRef)
        ! Arguments
        character(len=*), intent(in) :: testname
        character(len=*), intent(in) :: Var
        integer, intent(in) :: iTry         !< 
        integer, intent(in) :: iRef         !< 
        ! Variables
        character(len=255) :: InfoAbs
        if(iRef/=iTry) then
            write(InfoAbs,'(A,I0,A,I0)') trim(Var),iRef,'/',iTry
            call test_fail(testname,InfoAbs)
            STOP 
        else
            write(InfoAbs,'(A,A,I0)') trim(Var),' ok ',iRef
            call test_success(testname,InfoAbs)
        endif
    end subroutine

    subroutine  test_equal_i1(testname,Var,VecTry,VecRef,bTest,bPrintOnly,bPassed)
        ! Arguments
        character(len=*), intent(in) :: testname
        character(len=*), intent(in) :: Var
        integer, dimension(:), intent(in) :: VecTry         !< 
        integer, dimension(:), intent(in) :: VecRef         !< 
        logical, intent(in) :: bTest
        logical, intent(in) :: bPrintOnly
        logical, intent(out),optional :: bPassed
        ! Variables
        character(len=255) :: InfoAbs
        integer :: i,cpt
        ! 
        cpt=0
        do i=1,size(VecRef)
            if(VecRef(i)/=VecTry(i)) then
                cpt=cpt+1
            endif
        enddo
        if(cpt>0) then
            write(InfoAbs,'(A,I0)') trim(Var)//' Elements different: ',cpt
            if(present(bPassed)) then
                bPassed=.false.
            endif
        else
            write(InfoAbs,'(A)') trim(Var)//' reproduced to identity'
            if(present(bPassed)) then
                bPassed=.true.
            endif
        endif
        if(bPrintOnly) then
            print'(A)',trim(InfoAbs)
        endif
        if(bTest) then
            if(cpt>0) then
                call test_fail(testname,InfoAbs)
                STOP 
            else
                call test_success(testname,InfoAbs)
            endif
        endif
    end subroutine

    subroutine  test_almost_equal_0(testname,Var,Ref,Try,MINNORM,bStop,bPrint,bPassed)
        ! Arguments
        character(len=*), intent(in) :: testname
        character(len=*), intent(in) :: Var
        real(ReKi), intent(in) :: Ref         !< 
        real(ReKi), intent(in) :: Try         !< 
        real(ReKi), intent(in) :: MINNORM
        logical, intent(in) :: bStop
        logical, intent(in) :: bPrint
        logical, intent(out),optional :: bPassed
        ! Variables
        character(len=255) :: InfoAbs
        real(ReKi) :: delta
        integer :: cpt
        ! 
        cpt=0
        delta=abs(Ref-Try)
        if(delta>MINNORM) then
            write(InfoAbs,'(A,ES8.1E2,A,ES8.1E2,A,I0)') trim(Var)//' tol: ',MINNORM,', mean: ',delta,' - Failed:',cpt
            call test_fail(testname,InfoAbs,bPrint,bStop)
        else
            write(InfoAbs,'(A,ES8.1E2,A,ES8.1E2)') trim(Var)//' tol: ',MINNORM,', mean: ',delta
            call test_success(testname,InfoAbs,bPrint)
        endif
        if(present(bPassed)) then
            bPassed=delta>MINNORM
        endif
    end subroutine
    subroutine  test_almost_equal_1(testname,Var,VecRef,VecTry,MINNORM,bStop,bPrint,bPassed)
        ! Arguments
        character(len=*), intent(in) :: testname
        character(len=*), intent(in) :: Var
        real(ReKi), dimension(:), intent(in) :: VecRef         !< 
        real(ReKi), dimension(:), intent(in) :: VecTry         !< 
        real(ReKi), intent(in) :: MINNORM
        logical, intent(in) :: bStop
        logical, intent(in) :: bPrint
        logical, intent(out),optional :: bPassed
        ! Variables
        character(len=255) :: InfoAbs
        integer :: i,cpt
        real(ReKi) :: delta
        real(ReKi) :: delta_cum
        ! 
        cpt=0
        delta_cum=0.0_ReKi
        do i=1,size(VecRef,1)
            delta=abs(VecRef(i)-VecTry(i))
            delta_cum=delta_cum+delta
            if(delta>MINNORM) then
                cpt=cpt+1
            endif
        enddo
        delta_cum=delta_cum/size(VecRef)

        if(cpt>0) then
            write(InfoAbs,'(A,ES8.1E2,A,ES8.1E2,A,I0)') trim(Var)//' tol: ',MINNORM,', mean: ',delta_cum,' - Failed:',cpt
            call test_fail(testname,InfoAbs,bPrint,bStop)
        else
            write(InfoAbs,'(A,ES8.1E2,A,ES8.1E2)') trim(Var)//' tol: ',MINNORM,', mean: ',delta_cum
            call test_success(testname,InfoAbs,bPrint)
        endif
        if(present(bPassed)) then
            bPassed=(cpt==0)
        endif
    end subroutine
    subroutine  test_almost_equal_2(testname,Var,VecRef,VecTry,MINNORM,bStop,bPrint,bPassed)
        ! Arguments
        character(len=*), intent(in) :: testname
        character(len=*), intent(in) :: Var
        real(ReKi), dimension(:,:), intent(in) :: VecRef         !< 
        real(ReKi), dimension(:,:), intent(in) :: VecTry         !< 
        real(ReKi), intent(in) :: MINNORM
        logical, intent(in) :: bStop
        logical, intent(in) :: bPrint
        logical, intent(out),optional :: bPassed
        ! Variables
        real(ReKi), dimension(:),allocatable :: VecRef2    !< 
        real(ReKi), dimension(:),allocatable :: VecTry2   !<
        integer :: p, i,j,n1,n2,nCPs
        ! 
        n1 = size(VecRef,1); n2 = size(VecRef,2); nCPs=n1*n2
        allocate ( VecRef2 (n1*n2)  ) ; allocate ( VecTry2 (n1*n2)  ) 
        p=0
        do j=1,n2; do i=1,n1
            p=p+1
            VecRef2(p)=VecRef(i,j)
            VecTry2(p)=VecTry(i,j)
        enddo; enddo;
        call  test_almost_equal(testname,Var,VecRef2,VecTry2,MINNORM,bStop,bPrint,bPassed)
    end subroutine

   ! --------------------------------------------------------------------------------}
   ! --- Specific FVW tests 
   ! --------------------------------------------------------------------------------{
   !>
   subroutine Test_BiotSavart_Sgmt(testname, ErrStat, ErrMsg)
      character(len=*),     intent(in)  :: testname
      integer(IntKi)      , intent(out) :: ErrStat !< Error status of the operation
      character(ErrMsgLen), intent(out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
      real(ReKi), dimension(3) :: P1,P2,P3,CP
      real(ReKi), dimension(3) :: U1
      real(ReKi) :: SegGamma1 !< Circulation  [m^2/s]
      real(ReKi) :: RegParam1 !< 
      integer(IntKi) :: i1,i2
      integer(IntKi) :: RegFunction 
      integer(IntKi), parameter :: nSegTot  = 2
      integer(IntKi), parameter :: nSegPTot = 3
      integer(IntKi), parameter :: nCPsTot  = 1
      real(ReKi),     dimension(3,nCPsTot) :: CPs            !< Control points
      real(ReKi),     dimension(3,nSegPTot) :: SegPoints     !< Segment points
      integer(IntKi), dimension(2,nSegTot) :: SegConnct      !< Connectivity, indices of segments points iSeg1, iSeg2
      real(ReKi),     dimension(nSegTot)   :: SegGamma       !< Segment circulation
      real(ReKi),     dimension(nSegTot)   :: RegParam       !< Regularization parameter
      real(ReKi),     dimension(3,nCPsTot) :: Uind_out       !< Induced velocity vector - Side effects!!!
      real(ReKi),     dimension(3,4) :: CPs_test   !< 
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""
      ! --- Test that the two functions return the same values 
      P1=(/0.  ,0.,-1./)
      P2=(/0.  ,0., 1./)
      CPs_test(:,1) = (/ 0.0,  0., 0.0  /) ! Middle
      CPs_test(:,2) = P1                   ! Extremity
      CPs_test(:,3) = (/ 0.05, 0., -0.5 /) ! Close
      CPs_test(:,4) = (/ 10.,  0., 0.0  /) ! Far
      do i2 = 1, size(CPs_test,2)
         ! Segment param
         CP=CPs_test(:,i2)
         SegGamma1=1
         RegParam1=0.5
         ! One segment param
         SegConnct(:,1)=(/1,2/)
         SegPoints(:,1) = P1
         SegPoints(:,2) = P2
         SegGamma(:) = SegGamma1
         RegParam(:) = RegParam1
         CPs (:,1)   = CP
         do i1=1,5
            RegFunction = idRegVALID(i1)
            ! Method 1
            Uind_out =0.0_ReKi
            call ui_seg(1, 1, CPs, &
                  1, 1, SegPoints, SegConnct, SegGamma,  &
                  RegFunction, RegParam, Uind_out)
            ! Method 2
            call ui_seg_11(CP-P1, CP-P2, SegGamma1, RegFunction, RegParam1, U1)
            ! Test
            !print*,'Reg function', RegFunction, 'CP',CP
            !print*,'Uind_out',Uind_out
            !print*,'U1      ',U1
            call test_almost_equal(testname,'Uind method1/2', U1, Uind_out(:,1), 1e-4_ReKi, .true.,.true.)
            !call test_almost_equal('Uind method1/2', U1, Uind_out(:,1), 1e-4, .false.,.true.)
         enddo
      enddo

      ! --- Test that the two segments or one segment returns the same value
      P1=(/0.  ,0.,-1./)
      P2=(/0.  ,0., 1./)
      P3=(/0.  ,0., 0./)
      CPs_test(:,1) = (/ 0.0,  0., 0.0  /) ! Middle
      CPs_test(:,2) = P1                   ! Extremity
      CPs_test(:,3) = (/ 0.05, 0., -0.5 /) ! Close
      CPs_test(:,4) = (/ 100.,  0., -0.5  /) ! Far
      do i2 = 1,size(CPs_test,2)
         ! Segment param
         CP=CPs_test(:,i2)
         SegGamma1=1
         RegParam1=0.5
         ! One segment param
         SegConnct(:,1)=(/1,2/)
         SegConnct(:,2)=(/2,3/)
         SegPoints(:,1) = P1
         SegPoints(:,2) = P3
         SegPoints(:,3) = P2
         SegGamma(:) = SegGamma1
         RegParam(:) = RegParam1
         CPs (:,1)   = CP
         do i1=1,4 ! NOTE stopping at 4 since Offset is not linear
            RegFunction = idRegVALID(i1)
            ! Method 1
            Uind_out =0.0_ReKi
            call ui_seg(1, 1, CPs, &
                  1, 2, SegPoints, SegConnct, SegGamma,  &
                  RegFunction, RegParam, Uind_out)
            ! Method 2
            call ui_seg_11(CP-P1, CP-P2, SegGamma1, RegFunction, RegParam1, U1)
            !print*,'Reg function', RegFunction, 'CP',CP
            !print*,'Uind_out',Uind_out
            !print*,'U1      ',U1
            call test_almost_equal(testname,'Uind 1seg/2seg', U1, Uind_out(:,1), 1e-4_ReKi, .true.,.true.)
         enddo
      enddo
   end subroutine

   !>
   subroutine Test_BiotSavart_Part(testname, ErrStat, ErrMsg)
      character(len=*),     intent(in)  :: testname
      integer(IntKi)      , intent(out) :: ErrStat !< Error status of the operation
      character(ErrMsgLen), intent(out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
      real(ReKi), dimension(3) :: P1,CP
      real(ReKi), dimension(3) :: U1
      real(ReKi), dimension(3) :: PartAlpha1 !< Particle intensity alpha=om.dV [m^3/s]
      real(ReKi) :: RegParam1 !< 
      integer(IntKi) :: i1,i2
      integer(IntKi) :: RegFunction 
      integer(IntKi), parameter :: nPart  = 1
      integer(IntKi), parameter :: nCPs  = 1
      real(ReKi),     dimension(3,nCPs) :: CPs             !< Control points
      real(ReKi),     dimension(3,nPart):: PartPoints     !< Particle points
      real(ReKi),     dimension(3,nPart):: PartAlpha    !< Particle circulation
      real(ReKi),     dimension(nPart)  :: RegParam    !< Regularization parameter
      real(ReKi),     dimension(3,nCPs) :: Uind_out     !< Induced velocity vector - Side effects!!!
      real(ReKi),     dimension(3,4) :: CPs_test   !< 
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""
      ! --- Test that the two functions return the same values 
      P1=(/0.0, 0.0, -1.0 /)
      CPs_test(:,1) = (/ 0.0,  0., 0.0  /) ! Middle
      CPs_test(:,2) = P1    ! Extremity
      CPs_test(:,3) = (/ 0.01, 0.01, -0.9 /) ! Close
      CPs_test(:,4) = (/ 10.,  0., 0.0  /) ! Far
      do i1=1,3
         do i2 = 1, size(CPs_test,2)
            ! Segment param
            CP                = CPs_test(:,i2)
            PartAlpha1(1:2) = 0
            PartAlpha1(3  ) = 2
            RegParam1         = 0.5
            ! One segment param
            PartPoints(:,1) = P1
            PartAlpha(:,1)  = PartAlpha1
            RegParam(:)     = RegParam1
            CPs (:,1)       = CP
            RegFunction = idRegPartVALID(i1)
            ! Method 1
            Uind_out =0.0_ReKi
            call ui_part_nograd(nCPS, CPs, nPart, PartPoints, PartAlpha, RegFunction, RegParam, Uind_out)
            ! Method 2
            call ui_part_nograd_11(CP-P1, PartAlpha1, RegFunction, RegParam1, U1)
            ! Test
            !print*,'Reg function', RegFunction, 'CP',CP
            !print*,'Uind_out',Uind_out
            !print*,'U1      ',U1
            call test_almost_equal(testname,'Uind part method1/2', U1, Uind_out(:,1), 1e-4_ReKi, .true.,.true.)
         enddo
      enddo
   end subroutine Test_BiotSavart_Part

   !> This test compares calls using the tree algorithm and the direct N^2 evaluation
   subroutine Test_BiotSavart_PartTree(testname, ErrStat, ErrMsg)
      character(len=*),     intent(in)  :: testname
      integer(IntKi)      , intent(out) :: ErrStat !< Error status of the operation
      character(ErrMsgLen), intent(out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
      type(T_Tree) :: Tree
      real(ReKi), dimension(3) :: U_ref
      integer(IntKi) :: i1,i2,i3,k, iCP
      integer(IntKi) :: RegFunction
      integer(IntKi) :: nPart  = 1
      integer(IntKi) :: nCPs  = 1
      real(ReKi), dimension(:,:), allocatable :: CPs        !< Control points
      real(ReKi), dimension(:,:), allocatable :: PartPoints !< Particle points
      real(ReKi), dimension(:,:), allocatable :: PartAlpha  !< Particle circulation
      real(ReKi), dimension(:)  , allocatable :: RegParam   !< Regularization parameter
      real(ReKi), dimension(:,:), allocatable :: Uind1      !< Induced velocity vector - Side effects!!!
      real(ReKi), dimension(:,:), allocatable :: Uind2      !< Induced velocity vector - Side effects!!!
      real(ReKi) :: BranchFactor, BranchSmall
      real(ReKi),     dimension(3,5) :: CPs_test   !< 
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""
      BranchFactor = 2.0_ReKi !< Should be above1
      BranchSmall  = 0.0_ReKi
      RegFunction = 1

      ! --- Test with 0 particle
      nPart=0; nCPs= 1
      call alloc(nPart,nCPs)
      CPs(:,1) = (/0.0,0.0,0.0/)
      Uind1 =0.0_ReKi
      Uind2 =0.0_ReKi
      U_ref =0.0_ReKi
      call grow_tree_part(Tree, nPart, PartPoints, PartAlpha, RegFunction, RegParam, 0)
      !call print_tree(Tree)
      call ui_tree_part(Tree, nCPs, CPs, BranchFactor, BranchSmall,  Uind2, ErrStat, ErrMsg)
      call ui_part_nograd(nCPS, CPs, nPart, PartPoints, PartAlpha, RegFunction, RegParam, Uind1)
      ! Test
      call test_almost_equal(testname,'Uind tree 0 part', U_ref, Uind2(:,1), 1e-4_ReKi, .true.,.true.)
      call cut_tree(Tree)
      call dealloc()


      ! --- Test with 1 particle
      nPart=1; nCPs= 1
      call alloc(nPart,nCPs)
      CPs(:,1) = (/0.0,0.0,0.0/)
      PartPoints(1:3,1) = (/1.0,0.0,0.0/)
      U_ref =0.0_ReKi
      call grow_tree_part(Tree, nPart, PartPoints, PartAlpha, RegFunction, RegParam, 0)
      !call print_tree(Tree)
      call ui_tree_part(Tree, nCPs, CPs, BranchFactor, BranchSmall,  Uind2, ErrStat, ErrMsg)
      call ui_part_nograd(nCPS, CPs, nPart, PartPoints, PartAlpha, RegFunction, RegParam, Uind1)
      ! Test
      call test_almost_equal(testname,'Uind tree 1 part', Uind1, Uind2, 1e-4_ReKi, .true.,.true.)
      call cut_tree(Tree)
      !call print_tree(Tree)
      call dealloc()

      ! --- Test with 81 particles on different CPs, inside and outside the distribution of particles
      nPart=3*3**3; nCPs= 1
      call alloc(nPart,nCPs)
      k=0
      do i1 = -1,1,1
         do i2 = -1,1,1
            do i3 = -1,1,1
               ! NOTE: here we purposely duplicate a point, since since is a challenging case
               k=k+1; PartPoints(1:3,k) = (/ i1, i2, i3  /)
               k=k+1; PartPoints(1:3,k) = (/ i1, i2, i3  /)
               k=k+1; PartPoints(1:3,k) = (/ i1*1.2, i2*1.3, i3*1.1  /)
            enddo
         enddo
      enddo
      CPs_test(:,1) = (/ 0.0,  0., 0.0  /) ! Middle
      CPs_test(:,2) = (/ 1.0, 1.0, 1.0  /) ! Close to a cell center
      CPs_test(:,3) = PartPoints(:,5)      ! On a particle point
      CPs_test(:,4) = (/ 2.0, 2.0, 2.0 /)  ! Starts to be far from most points
      CPs_test(:,5) = (/ 10., 10., 10.0  /) ! Far from all

      call grow_tree_part(Tree, nPart, PartPoints, PartAlpha, RegFunction, RegParam, 0)
      !call print_tree(Tree)
      do iCP=1,4
         CPs(:,1) = CPs_test(:,icp)
         Uind2=0.0_ReKi; Uind1=0.0_ReKi
         call ui_tree_part(Tree, nCPs, CPs, BranchFactor, BranchSmall, Uind2, ErrStat, ErrMsg)
         call ui_part_nograd(nCPs, CPs, nPart, PartPoints, PartAlpha, RegFunction, RegParam, Uind1)
         !print*,'Uind',Uind1, Uind2
         ! Test
         call test_almost_equal(testname,'Uind tree 81 part', Uind1, Uind2, 1e-2_ReKi, .true.,.true.)
      enddo
      call cut_tree(Tree)
      ! --- Test that tree ui cannot be called after tree has been cut
      call ui_tree_part(Tree, nCPs, CPs, BranchFactor, BranchSmall, Uind2, ErrStat, ErrMsg)
      call test_equal(testname,'Err. stat tree cut',ErrStat,ErrID_Fatal)
      call dealloc()

   contains
      subroutine alloc(nPart, nCPs)
         integer(IntKi) :: nPart, nCPs
         allocate(PartPoints(3,nPart), PartAlpha(3,nPart), RegParam(nPart))
         allocate(CPs(3,nCPs), Uind1(3,nCPs), Uind2(3,nCPs))
         RegParam(:)=0.01
         PartAlpha(1,:)  = 0.0
         PartAlpha(2,:)  = 0.0
         PartAlpha(3,:)  = 1.0
         Uind1 =0.0_ReKi
         Uind2 =0.0_ReKi
      end subroutine 
      subroutine dealloc()
         deallocate(PartPoints, PartAlpha, RegParam)
         deallocate(CPs, Uind1, Uind2)
      end subroutine 
   end subroutine Test_BiotSavart_PartTree

   !> Compares the velocity field obtained from a segment and its convert to particle version
   subroutine Test_SegmentsToPart(testname, ErrStat, ErrMsg)
      character(len=*),     intent(in)  :: testname
      integer(IntKi)      , intent(out) :: ErrStat !< Error status of the operation
      character(ErrMsgLen), intent(out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
      real(ReKi),     dimension(:,:), allocatable :: PartPoints   !< Particle points
      real(ReKi),     dimension(:,:), allocatable :: PartAlpha    !< Particle circulation
      real(ReKi),     dimension(:)  , allocatable :: PartEpsilon  !< Regularization parameter
      integer(IntKi), parameter :: nSegTot  = 2
      integer(IntKi), parameter :: nSegPTot = 3
      integer(IntKi), parameter :: nCPsTot  = 10
      real(ReKi),     dimension(3,nSegPTot) :: SegPoints     !< Segment points
      integer(IntKi), dimension(2,nSegTot) :: SegConnct      !< Connectivity, indices of segments points iSeg1, iSeg2
      real(ReKi),     dimension(nSegTot)   :: SegGamma       !< Segment circulation
      real(ReKi),     dimension(nSegTot)   :: SegEpsilon     !< Regularization parameter
      real(ReKi),     dimension(3,nCPsTot) :: CPs            !< Control points
      real(ReKi),     dimension(3,nCPsTot) :: Uind1          !< Induced velocity vector - Side effects!!!
      real(ReKi),     dimension(3,nCPsTot) :: Uind2          !< Induced velocity vector - Side effects!!!
      real(ReKi) :: RegParam1 !< 
      integer(IntKi) :: nPartPerSeg, nPart, iHeadP
      integer(IntKi) :: RegFunctionPart, RegFunctionSeg
      ErrStat = ErrID_None
      ErrMsg  = ""
      RegParam1=1.0
      ! Creating two aligned segments
      SegConnct(:,1)=(/1,2/)
      SegConnct(:,2)=(/2,3/)
      SegPoints(:,1) = (/0.  ,0.,-1./)
      SegPoints(:,2) = (/0.  ,0., 0./)
      SegPoints(:,3) = (/0.  ,0., 1./)
      SegGamma(:)    =4
      SegEpsilon = RegParam1
      ! Points where velocity will be evaluated
      CPs(:,1) = SegPoints(:,1)
      CPs(:,2) = SegPoints(:,2)
      CPs(:,3) = SegPoints(:,3)
      CPs(:,4) = (/ 0.2, 0.2,   0.0/)
      CPs(:,6) = (/ 0.5, 0.5,   0. /)
      CPs(:,8) = (/ 1.0, 1.0,   0./)
      CPs(:,9) = (/ 10.0, 10.0, 0./)
      CPs(:,5) = (/ 0.2, 0.2,   0.5/)
      CPs(:,7) = (/ 0.5, 0.5,   0.5/)
      CPs(:,10) = (/ 1.0, 1.0,   1./)

      ! --- Test 1 - 10 particles, no regularization
      RegFunctionSeg  = idRegNone
      RegFunctionPart = idRegNone
      nPartPerSeg = 10

      nPart = nPartPerSeg * nSegTot
      call alloc(nPart)
      iHeadP=1
      call SegmentsToPart(SegPoints, SegConnct, SegGamma, SegEpsilon, 1, nSegTot, nPartPerSeg, PartPoints, PartAlpha, PartEpsilon, iHeadP)

      Uind1 =0.0_ReKi; Uind2 =0.0_ReKi;
      call ui_seg(1, nCPsTot, CPs, 1, nSegTot, SegPoints, SegConnct, SegGamma, RegFunctionSeg, SegEpsilon, Uind1)
      call ui_part_nograd(nCPSTot, CPs, nPart, PartPoints, PartAlpha, RegFunctionPart, PartEpsilon, Uind2)
      call test_almost_equal(testname,'Uind 10 part/sgmt no reg', Uind1, Uind2, 1e-3_ReKi, .true.,.true.)
      call dealloc()

      ! --- Test 1 - 2 particles, no regularization
      RegFunctionSeg  = idRegNone
      RegFunctionPart = idRegNone
      nPartPerSeg = 2

      nPart = nPartPerSeg * nSegTot
      call alloc(nPart)
      iHeadP=1
      call SegmentsToPart(SegPoints, SegConnct, SegGamma, SegEpsilon, 1, nSegTot, nPartPerSeg, PartPoints, PartAlpha, PartEpsilon, iHeadP)

      Uind1 =0.0_ReKi; Uind2 =0.0_ReKi;
      call ui_seg(1, nCPsTot, CPs, 1, nSegTot, SegPoints, SegConnct, SegGamma, RegFunctionSeg, SegEpsilon, Uind1)
      call ui_part_nograd(nCPsTot, CPs, nPart, PartPoints, PartAlpha, RegFunctionPart, PartEpsilon, Uind2)
      call test_almost_equal(testname,'Uind 2 part/sgmt noreg', Uind1, Uind2, 3e-1_ReKi, .true.,.true.)
      call dealloc()


      ! --- Test 3 - 10 particles, regularization 
      ! NOTE: more work needed to match the regularization functions and parameters optimally
      RegFunctionSeg = idRegLambOseen
      RegFunctionPart = idRegExp
      nPartPerSeg = 10

      nPart = nPartPerSeg * nSegTot
      call alloc(nPart)
      iHeadP=1
      call SegmentsToPart(SegPoints, SegConnct, SegGamma, SegEpsilon, 1, nSegTot, nPartPerSeg, PartPoints, PartAlpha, PartEpsilon, iHeadP)

      Uind1 =0.0_ReKi; Uind2 =0.0_ReKi;
      call ui_seg(1, nCPsTot, CPs, 1, nSegTot, SegPoints, SegConnct, SegGamma, RegFunctionSeg, SegEpsilon, Uind1)
      call ui_part_nograd(nCPSTot, CPs, nPart, PartPoints, PartAlpha, RegFunctionPart, PartEpsilon, Uind2)
      !print'(A,10F7.3)','Uind1',Uind1(1,:)
      !print'(A,10F7.3)','Uind2',Uind2(1,:)
      !print'(A,10F7.3)','Uind1',Uind1(2,:)
      !print'(A,10F7.3)','Uind2',Uind2(2,:)
      !print'(A,10F7.3)','Uind1',Uind1(3,:)
      !print'(A,10F7.3)','Uind2',Uind2(3,:)
      call test_almost_equal(testname,'Uind 10 part/sgmt w.reg', Uind1, Uind2, 5e-2_ReKi, .true.,.true.)
      call dealloc()

   contains 
      subroutine alloc(n)
         integer(IntKi) :: n
         allocate(PartPoints(3,n), PartAlpha(3,n), PartEpsilon(n))
         PartAlpha(:,:)  = -99999.99_ReKi
         PartPoints(:,:) = -99999.99_ReKi
         PartEpsilon(:)  = -99999.99_ReKi
      end subroutine 
      subroutine dealloc()
         deallocate(PartPoints, PartAlpha, PartEpsilon)
      end subroutine 
   end subroutine Test_SegmentsToPart

   !>
   subroutine Test_LatticeToSegment(mvtk,iStat)
      type(VTK_Misc),intent(inout) :: mvtk       !< miscvars for VTK output
      integer(IntKi), intent(  out)  :: iStat !< Status for test
      ! Local
      integer(IntKi),dimension(:,:), allocatable :: SegConnct !< Segment connectivity
      real(ReKi),    dimension(:,:), allocatable :: SegPoints !< Segment Points
      real(ReKi),    dimension(:)  , allocatable :: SegGamma  !< Segment Circulation
      real(ReKi),    dimension(:),   allocatable :: SegEpsilon !< 
      !
      real(ReKi),    dimension(:,:,:), allocatable :: LatticePoints1 !< Lattice Points
      real(ReKi),    dimension(:,:,:), allocatable :: LatticePoints2 !< Lattice Points
      real(ReKi),    dimension(:,:),   allocatable :: LatticeGamma1  !< Lattice Circulation
      real(ReKi),    dimension(:,:),   allocatable :: LatticeGamma2  !< Lattice Circulation
      real(ReKi),    dimension(:,:,:), allocatable :: LatticeEps1 !< Lattice Reg Param
      real(ReKi),    dimension(:,:,:), allocatable :: LatticeEps2 !< Lattice Reg Param
      real(ReKi),    dimension(:,:),   allocatable :: CPs   !< ControlPoints
      real(ReKi),    dimension(:,:),   allocatable :: Uind  !< Induced velocity
      integer(IntKi) :: iHeadC
      integer(IntKi) :: iHeadP
      integer(IntKi) :: i,j
      integer(IntKi) :: nP
      integer(IntKi) :: nC
      integer(IntKi) :: nP1, nP2
      integer(IntKi) :: nC1, nC2
      integer(IntKi) :: nDepth, nSpan
      integer(IntKi) :: SmoothModel
      logical :: bladeFrame   !< Output in blade frame instead of global coordinate frame
      iStat=0
      bladeFrame=.FALSE.

      ! --- Creating two lattice
      allocate(LatticePoints1(3,2,2)) 
      allocate(LatticePoints2(3,3,4)) 
      allocate(LatticeEps1(3,1,1)) 
      allocate(LatticeEps2(3,2,3)) 
      allocate(LatticeGamma1(1,1)) ; 
      allocate(LatticeGamma2(2,3)) ; 
      LatticeGamma1=1
      LatticeEps1(1,:,:) = 1
      LatticeEps1(2,:,:) = 2
      LatticeEps1(3,:,:) = 3
      ! Test shed vorticity
      LatticeGamma2(:,1)=1
      LatticeGamma2(:,2)=2
      LatticeGamma2(:,3)=3

      LatticeEps2(:,:,1) = 1
      LatticeEps2(:,:,2) = 2
      LatticeEps2(:,:,3) = 3
      ! Test trailed vorticity
!       LatticeGamma2(1,:)=1
!       LatticeGamma2(2,:)=2
      CALL MeshMe(LatticePoints1,(/0.0_ReKi,0.0_ReKi,0.0_ReKi/))
      CALL MeshMe(LatticePoints2,(/0.0_ReKi,0.0_ReKi,1.0_ReKi/))

      CALL WrVTK_Lattice('Points1.vtk',mvtk,LatticePoints1, LatticeGamma1, bladeframe=bladeframe)
      CALL WrVTK_Lattice('Points2.vtk',mvtk,LatticePoints2, LatticeGamma2, bladeframe=bladeframe)

      ! --- Convert lattice 1 to segments
      nSpan  = size(LatticePoints1,2)
      nDepth = size(LatticePoints1,3)
      nP1 = nSpan*nDepth
      nC1 = 2*(nSpan*nDepth)-nSpan-nDepth
      allocate(SegConnct(1:2,1:nC1)); SegConnct=-1
      allocate(SegPoints(1:3,1:nP1)); SegPoints=-1
      allocate(SegGamma (1:nC1)    ); SegGamma=-999
      allocate(SegEpsilon(1:nC1)    ); SegEpsilon=0.0_ReKi

      iHeadP=1
      iHeadC=1
      CALL LatticeToSegments(LatticePoints1, LatticeGamma1, LatticeEps1, 1, SegPoints, SegConnct, SegGamma, SegEpsilon, iHeadP, iHeadC, .true., .true., .false. )
      CALL printall()
      CALL WrVTK_Segments('Points1_seg.vtk', mvtk, SegPoints, SegConnct, SegGamma, SegEpsilon, bladeFrame) 

      allocate(Uind(1:3,1) ); Uind=0.0_ReKi
      allocate(CPs (1:3,1) ); 
      CPs(1:3,1)=(/1.5,1.5,0./)
      SegEpsilon=100.0_ReKi
      SmoothModel=0 ! No smooth
      CALL ui_seg(1, 1, CPs, &
      1, nC1, SegPoints, SegConnct, SegGamma,   &
      SmoothModel, SegEpsilon, Uind)
      !print*,'Uind',Uind

      ! --- Convert lattice 2 to segments
      nSpan  = size(LatticePoints2,2)
      nDepth = size(LatticePoints2,3)
      nP2 = nSpan*nDepth
      nC2 = 2*(nSpan*nDepth)-nSpan-nDepth
      deallocate(SegConnct)
      deallocate(SegPoints)
      deallocate(SegGamma)
      allocate(SegConnct(1:2,1:nC2)); SegConnct=-1
      allocate(SegPoints(1:3,1:nP2)); SegPoints=-1
      allocate(SegGamma (1:nC2)    ); SegGamma=-9999
      iHeadP=1
      iHeadC=1
      CALL LatticeToSegments(LatticePoints2, LatticeGamma2,  LatticeEps2, 1, SegPoints, SegConnct, SegGamma, SegEpsilon, iHeadP, iHeadC , .true., .true., .false.)
      CALL printall()
      CALL WrVTK_Segments('Points2_seg.vtk', mvtk, SegPoints, SegConnct, SegGamma, SegEpsilon, bladeFrame) 

      ! --- Concatenate both
      nP = nP1 + nP2
      nC = nC1 + nC2
      iHeadP=1
      iHeadC=1
      deallocate(SegConnct)
      deallocate(SegPoints)
      deallocate(SegGamma)
      allocate(SegConnct(1:2,1:nC)); SegConnct=-1
      allocate(SegPoints(1:3,1:nP)); SegPoints=-1
      allocate(SegGamma (1:nC)    ); SegGamma=-9999
      CALL LatticeToSegments(LatticePoints1, LatticeGamma1, LatticeEps1, 1, SegPoints, SegConnct, SegGamma, SegEpsilon, iHeadP, iHeadC, .true. , .true., .false.)
      CALL LatticeToSegments(LatticePoints2, LatticeGamma2, LatticeEps2, 1, SegPoints, SegConnct, SegGamma, SegEpsilon, iHeadP, iHeadC, .true. , .true., .false.)
      CALL printall()
      CALL WrVTK_Segments('PointsBoth_seg.vtk', mvtk, SegPoints, SegConnct, SegGamma, SegEpsilon, bladeFrame) 


   contains
      subroutine printall()
         print*,'Points'
         do i=1,size(SegPoints,2)
            print*,'i',i,'Coords:', SegPoints(1:3,i)
         enddo
         print*,'Connectivity'
         do i=1,size(SegConnct,2)
            print*,'i',i,'Conn:', SegConnct(1:2,i),'Gam:', SegGamma(i)
         enddo
         print*,'-----------------------------'
      endsubroutine

      subroutine MeshMe(M,offset)
         real(ReKi), dimension(:,:,:), intent(inout) :: M
         real(ReKi), dimension(3)    , intent(in   ):: offset
         do j=1,size(M,3)
            do i=1,size(M,2)
               M(1,i,j)=i + offset(1)
               M(2,i,j)=j + offset(2)
               M(3,i,j)=0 + offset(3)
            enddo
         enddo 
      endsubroutine 
   endsubroutine Test_LatticeToSegment

   !> Test Wake Induced velocity calcualtion when using nNWMax or nNWFree
   !! A dummy helical wake is created. The induced velocity is computed on 
   !! either the full wake, or just the "free" wake (which should be way faster)
   subroutine FVW_Test_WakeInducedVelocities(ErrStat, ErrMsg)
      integer(IntKi)      , intent(out) :: ErrStat !< Error status of the operation
      character(ErrMsgLen), intent(out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
      type(FVW_ParameterType)       :: p !< Parameters
      type(FVW_ContinuousStateType) :: x !< States
      type(FVW_MiscVarType)         :: m !< Initial misc/optimization variables
      !type(FVW_VTK_Misc)   :: mvtk
      integer :: iW, j, k, nSpan
      integer(IntKi)       :: ErrStat2
      character(ErrMsgLen) :: ErrMsg2
      character(*), parameter  :: RoutineName = 'FVW_Test_CPUTime'
      integer(ReKi), parameter :: nR          = 20
      real(ReKi), parameter    :: R           = 100
      real(ReKi), parameter    :: G           = 100
      real(ReKi), allocatable, dimension(:,:) :: V1
      real(ReKi), allocatable, dimension(:,:) :: V2
      real(ReKi) :: t1,t2 
      ErrStat = ErrID_None
      ErrMsg  = ""

      ! --- Create a helical wake TODO, put me into FVW_*
      p%nWings  = 3
      p%nNWMax           = 1600
      nSpan              = 50
      m%nNW              = p%nNWMax
      m%nFW              = 0
      p%nFWMax           = 0
      p%nFWFree          = 0
      p%ShearModel       = idShearNone
      p%RegFunction      = idRegVatistas
      p%VelocityMethod   = idVelocityTreePart
      p%FWShedVorticity  = .false.
      p%TreeBranchFactor = 1.5_ReKi
      p%PartPerSegment   = 1
      allocate(p%W(p%nWings))
      p%W(:)%nSpan       = nSpan
      call FVW_InitStates( x, p, ErrStat, ErrMsg )
      do iW=1,size(x%W); 
         do j=1,size(x%W(iW)%r_NW,2); 
            do k=1,size(x%W(iW)%r_NW,3); 
               x%W(iW)%r_NW(1,j,k) = real(k, ReKi)/p%nNWMax*(nR*R)
               x%W(iW)%r_NW(2,j,k) = real(j, ReKi)/nSpan*R*cos(iW*TwoPi/p%nWings + x%W(iW)%r_NW(1,j,k)/R*0.5)
               x%W(iW)%r_NW(3,j,k) = real(j, ReKi)/nSpan*R*sin(iW*TwoPi/p%nWings + x%W(iW)%r_NW(1,j,k)/R*0.5) + 1.5*R
            enddo
         enddo
         do j=1,size(x%W(iW)%r_NW,2)-1 
            do k=1,size(x%W(iW)%r_NW,3)-1
               x%W(iW)%Gamma_NW(j,k) = G*4.0_ReKi*(real((j-1),ReKi)/nSpan -0.5)**2
               x%W(iW)%Eps_NW(:,j,k) = 0.03*R*(real(k,ReKi)/p%nNWMax)
            enddo
         enddo
      enddo
      allocate(m%W(p%nWings))
      do iW = 1,p%nWings
         call AllocAry( m%W(iW)%Vind_NW , 3   ,  nSpan+1  ,p%nNWMax+1, 'Vind on NW ', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName); m%W(iW)%Vind_NW= -999_ReKi;
         call AllocAry( m%W(iW)%Vind_FW , 3   ,  FWnSpan+1,p%nFWMax+1, 'Vind on FW ', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName); m%W(iW)%Vind_FW= -999_ReKi;
      enddo
      call FVW_InitMiscVarsPostParam( p, m, ErrStat2, ErrMsg2) ! Alloc Sgmt, CPs, Uind

      ! --- Compute induced velocity on full wake
      allocate(V1(3,nSpan+1))
      p%nNWFree=p%nNWMax
      call cpu_time(t1)
      call WakeInducedVelocities(p, x, m, ErrStat2, ErrMsg2); 
      call cpu_time(t2)
      print*,'Ellapsed time',t2-t1
      V1 = m%W(1)%Vind_NW(:,:,1)

      ! --- Compute induced velocity on free wake only
      allocate(V2(3,nSpan+1))
      p%nNWFree=int(p%nNWMax/5)
      call cpu_time(t1)
      call WakeInducedVelocities(p, x, m, ErrStat2, ErrMsg2); 
      call cpu_time(t2)
      print*,'Ellapsed time',t2-t1
      V2 = m%W(1)%Vind_NW(:,:,1)
      !print*,'>>>Vx mean ',sum(abs((V1(1,:))))/nSpan
      !print*,'>>>Vx mean ',sum(abs((V2(1,:))))/nSpan
      !print*,'>>> Vx error',sum(abs(V1(1,:)-V2(1,:)))/nSpan
      !print*,'>>> Vy error',sum(abs(V1(2,:)-V2(2,:)))/nSpan
      !print*,'>>> Vz error',sum(abs(V1(3,:)-V2(3,:)))/nSpan
      !call WrVTK_Segments('_TEST.vtk', mvtk, m%Sgmt%Points(:,:), m%Sgmt%Connct(:,:), m%Sgmt%Gamma(:), m%Sgmt%Epsilon(:), .false.) 

      call test_almost_equal(RoutineName,'Uind nNW/nNWFree', V1, V2, 1e-6_ReKi, .true.,.true.)

      deallocate(V1)
      deallocate(V2)
      call FVW_DestroyParam(p, ErrStat2, ErrMsg2)
      call FVW_DestroyContState(x, ErrStat2, ErrMsg2)
      call FVW_DestroyMisc(m, ErrStat2, ErrMsg2)

   end subroutine FVW_Test_WakeInducedVelocities

   !> Main test function 
   subroutine FVW_RunTests(ErrStat,ErrMsg)
      integer(IntKi)      , intent(out) :: ErrStat !< Error status of the operation
      character(ErrMsgLen), intent(out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
      integer(IntKi)       :: ErrStat2
      character(ErrMsgLen) :: ErrMsg2
      character(len=255)   :: testname
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""
      testname='FVW'
      call Test_BiotSavart_Sgmt(testname, ErrStat2, ErrMsg2)
      call Test_BiotSavart_Part(testname, ErrStat2, ErrMsg2)
      call Test_BiotSavart_PartTree(testname, ErrStat2, ErrMsg2)
      call Test_SegmentsToPart(testname, ErrStat2, ErrMsg2)
      call FVW_Test_WakeInducedVelocities(ErrStat2, ErrMsg2)
   end subroutine FVW_RunTests

end module FVW_Tests
