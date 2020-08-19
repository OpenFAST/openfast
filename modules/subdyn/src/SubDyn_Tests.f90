module SubDyn_Tests
   use NWTC_Library
   use SubDyn_Types
   use SD_FEM
   use IntegerList
   
   implicit none

   public :: SD_Tests
   private

   character(len=255),save :: testname
   interface test_equal; module procedure &
         test_equal_i1, &
         test_equal_i0
   end interface
    interface test_almost_equal; module procedure &
          test_almost_equal_0, &
          test_almost_equal_1, &
          test_almost_equal_1d, &
          test_almost_equal_2,  &
          test_almost_equal_2d  
    end interface
contains

   ! --------------------------------------------------------------------------------
   ! --- Helper functions (should be part of NWTC library)
   ! --------------------------------------------------------------------------------
    subroutine test_success(info,bPrint_in)
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

    subroutine test_fail(info,bPrint_in,bStop_in)
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

    subroutine  test_equal_i0(Var,iTry,iRef)
        ! Arguments
        character(len=*), intent(in) :: Var
        integer, intent(in) :: iTry         !< 
        integer, intent(in) :: iRef         !< 
        ! Variables
        character(len=255) :: InfoAbs
        if(iRef/=iTry) then
            write(InfoAbs,'(A,I0,A,I0)') trim(Var),iRef,'/',iTry
            call test_fail(InfoAbs)
            STOP 
        else
            write(InfoAbs,'(A,A,I0)') trim(Var),' ok ',iRef
            call test_success(InfoAbs)
        endif
    end subroutine

    subroutine  test_equal_i1(Var,VecTry,VecRef,bTest,bPrintOnly,bPassed)
        ! Arguments
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
                call test_fail(InfoAbs)
                STOP 
            else
                call test_success(InfoAbs)
            endif
        endif
    end subroutine

    subroutine  test_almost_equal_0(Var,Ref,Try,MINNORM,bStop,bPrint,bPassed)
        ! Arguments
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
            call test_fail(InfoAbs,bPrint,bStop)
        else
            write(InfoAbs,'(A,ES8.1E2,A,ES8.1E2)') trim(Var)//' tol: ',MINNORM,', mean: ',delta
            call test_success(InfoAbs,bPrint)
        endif
        if(present(bPassed)) then
            bPassed=delta>MINNORM
        endif
    end subroutine
    subroutine  test_almost_equal_1(Var,VecRef,VecTry,MINNORM,bStop,bPrint,bPassed)
        ! Arguments
        character(len=*), intent(in) :: Var
        real(SiKi), dimension(:), intent(in) :: VecRef         !< 
        real(SiKi), dimension(:), intent(in) :: VecTry         !< 
        real(SiKi), intent(in) :: MINNORM
        logical, intent(in) :: bStop
        logical, intent(in) :: bPrint
        logical, intent(out),optional :: bPassed
        ! Variables
        character(len=255) :: InfoAbs
        integer :: i,cpt
        real(SiKi) :: delta
        real(SiKi) :: delta_cum
        ! 
        cpt=0
        delta_cum=0.0_SiKi
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
            call test_fail(InfoAbs,bPrint,bStop)
        else
            write(InfoAbs,'(A,ES8.1E2,A,ES8.1E2)') trim(Var)//' tol: ',MINNORM,', mean: ',delta_cum
            call test_success(InfoAbs,bPrint)
        endif
        if(present(bPassed)) then
            bPassed=(cpt==0)
        endif
    end subroutine
    subroutine  test_almost_equal_1d(Var,VecRef,VecTry,MINNORM,bStop,bPrint,bPassed)
        ! Arguments
        character(len=*), intent(in) :: Var
        real(R8Ki), dimension(:), intent(in) :: VecRef         !< 
        real(R8Ki), dimension(:), intent(in) :: VecTry         !< 
        real(R8Ki), intent(in) :: MINNORM
        logical, intent(in) :: bStop
        logical, intent(in) :: bPrint
        logical, intent(out),optional :: bPassed
        ! Variables
        character(len=255) :: InfoAbs
        integer :: i,cpt
        real(R8Ki) :: delta
        real(R8Ki) :: delta_cum
        ! 
        cpt=0
        delta_cum=0.0_R8Ki
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
            call test_fail(InfoAbs,bPrint,bStop)
        else
            write(InfoAbs,'(A,ES8.1E2,A,ES8.1E2)') trim(Var)//' tol: ',MINNORM,', mean: ',delta_cum
            call test_success(InfoAbs,bPrint)
        endif
        if(present(bPassed)) then
            bPassed=(cpt==0)
        endif
    end subroutine
    subroutine  test_almost_equal_2(Var,VecRef,VecTry,MINNORM,bStop,bPrint,bPassed)
        ! Arguments
        character(len=*), intent(in) :: Var
        real(SiKi), dimension(:,:), intent(in) :: VecRef         !< 
        real(SiKi), dimension(:,:), intent(in) :: VecTry         !< 
        real(SiKi), intent(in) :: MINNORM
        logical, intent(in) :: bStop
        logical, intent(in) :: bPrint
        logical, intent(out),optional :: bPassed
        ! Variables
        real(SiKi), dimension(:),allocatable :: VecRef2    !< 
        real(SiKi), dimension(:),allocatable :: VecTry2   !<
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
        call  test_almost_equal(Var,VecRef2,VecTry2,MINNORM,bStop,bPrint,bPassed)
    end subroutine
    subroutine  test_almost_equal_2d(Var,VecRef,VecTry,MINNORM,bStop,bPrint,bPassed)
        ! Arguments
        character(len=*), intent(in) :: Var
        real(R8Ki), dimension(:,:), intent(in) :: VecRef         !< 
        real(R8Ki), dimension(:,:), intent(in) :: VecTry         !< 
        real(R8Ki), intent(in) :: MINNORM
        logical, intent(in) :: bStop
        logical, intent(in) :: bPrint
        logical, intent(out),optional :: bPassed
        ! Variables
        real(R8Ki), dimension(:),allocatable :: VecRef2    !< 
        real(R8Ki), dimension(:),allocatable :: VecTry2   !<
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
        call  test_almost_equal(Var,VecRef2,VecTry2,MINNORM,bStop,bPrint,bPassed)
    end subroutine



   ! --------------------------------------------------------------------------------}
   ! --- Specific SubDyn tests 
   ! --------------------------------------------------------------------------------{
   subroutine Test_CB_Results(MBBt, MBMt, KBBt, OmegaM, DOFTP, DOFM, ErrStat, ErrMsg)
      INTEGER(IntKi)                                     :: DOFTP, DOFM
      REAL(ReKi)                                         :: MBBt(DOFTP, DOFTP)
      REAL(ReKi)                                         :: MBmt(DOFTP, DOFM)
      REAL(ReKi)                                         :: KBBt(DOFTP, DOFTP)
      REAL(ReKi)                                         :: OmegaM(DOFM)
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
      ! local variables
      INTEGER(IntKi) :: DOFT, NM, i
      REAL(ReKi), Allocatable     :: OmegaCB(:), PhiCB(:, :)
      REAL(ReKi), Allocatable     :: K(:, :)
      REAL(ReKi), Allocatable     :: M(:, :)
      Character(1024)             :: rootname
      ErrStat = ErrID_None
      ErrMsg  = ''
      print*,'This test is not a unit test'
      
      DOFT = DOFTP + DOFM
      NM = DOFT - 3
      Allocate( OmegaCB(NM), K(DOFT, DOFT), M(DOFT, DOFT), PhiCB(DOFT, NM) )
      K = 0.0
      M = 0.0
      OmegaCB = 0.0
      PhiCB = 0.0
      
      M(1:DOFTP, 1:DOFTP) = MBBt
      M(1:DOFTP, (DOFTP+1):DOFT ) = MBMt
      M((DOFTP+1):DOFT, 1:DOFTP ) = transpose(mbmt)

      DO i = 1, DOFM
         K(DOFTP+i, DOFTP+i) = OmegaM(i)*OmegaM(i)
         M(DOFTP+i, DOFTP+i) = 1.0
      ENDDO
         
      K(1:DOFTP, 1:DOFTP) = KBBt

      ! temporary rootname
      rootname = './test_assemble_C-B_out'
      
      ! NOTE: Eigensolve is in SubDyn
      !CALL EigenSolve(K, M, DOFT, NM,.False.,Init,p, PhiCB, OmegaCB,  ErrStat, ErrMsg)
      IF ( ErrStat /= 0 ) RETURN  
   end subroutine Test_CB_Results

   !> Transformation matrices tests
   subroutine Test_Transformations(ErrStat,ErrMsg)
      integer(IntKi)      , intent(out) :: ErrStat
      character(ErrMsgLen), intent(out) :: ErrMsg

      real(ReKi), dimension(3) :: P1, P2, e1, e2, e3
      real(FEKi), dimension(3,3) :: DirCos, Ref
      real(ReKi), dimension(6,6) :: T, Tref
      real(ReKi) :: L
      integer(IntKi) :: I
      testname='Transf'

      ! --- DirCos
      P1=(/0,0,0/)
      P2=(/2,0,0/)
      call GetDirCos(P1, P2, DirCos, L, ErrStat, ErrMsg)
      Ref = reshape( (/0_FEKi,-1_FEKi,0_FEKi, 0_FEKi, 0_FEKi, -1_FEKi, 1_FEKi, 0_FEKi, 0_FEKi/) , (/3,3/))
      call  test_almost_equal('DirCos',Ref,DirCos,1e-8_FEKi,.true.,.true.)

      ! --- Rigid Transo
      P1=(/1,2,-1/)
      P2=(/2,5,5/)
      call GetRigidTransformation(P1, P2, T, ErrStat, ErrMsg)
      Tref = 0; do I=1,6; Tref(I,I) =1_ReKi; enddo
      Tref(1,5) = 6._ReKi; Tref(1,6) =-3._ReKi;
      Tref(2,4) =-6._ReKi; Tref(2,6) = 1._ReKi;
      Tref(3,4) = 3._ReKi; Tref(3,5) =-1._ReKi;
      call  test_almost_equal('TRigid',Tref,T,1e-8_ReKi,.true.,.true.)
      

      ! --- Orthogonal vectors
      e1 = (/10,0,0/)
      call GetOrthVectors(e1,e2,e3,ErrStat, ErrMsg)
      call  test_almost_equal('orth',e2,(/0._ReKi,0._ReKi,-1._ReKi/),1e-8_ReKi,.true.,.true.)
      call  test_almost_equal('orth',e3,(/0._ReKi,1._ReKi, 0._ReKi/),1e-8_ReKi,.true.,.true.)
      e1 = (/0,10,0/)
      call GetOrthVectors(e1,e2,e3,ErrStat, ErrMsg)
      call  test_almost_equal('orth',e2,(/0._ReKi,0._ReKi, 1._ReKi/),1e-8_ReKi,.true.,.true.)
      call  test_almost_equal('orth',e3,(/1._ReKi,0._ReKi, 0._ReKi/),1e-8_ReKi,.true.,.true.)
      e1 = (/1,2,4/)
      call GetOrthVectors(e1,e2,e3,ErrStat, ErrMsg)
      call test_almost_equal('dot', 0._ReKi, dot_product(e1,e2),  1e-8_ReKi, .true., .true.)
      call test_almost_equal('dot', 0._ReKi, dot_product(e1,e3),  1e-8_ReKi, .true., .true.)
   end subroutine  Test_Transformations


   !> Linear algebra tests
   subroutine Test_Linalg(ErrStat,ErrMsg)
      integer(IntKi)      , intent(out) :: ErrStat
      character(ErrMsgLen), intent(out) :: ErrMsg
      real(FEKi), dimension(:,:), allocatable :: A, Ainv, Aref
      real(DbKi) :: det
      testname='Linalg'

      ! --- Determinant of a singular matrix
      ! Commented since might lead to floating invalid
      !allocate(A(3,3));
      !A(1,1) = 0 ; A(1,2) = 0 ; A(1,3) =  1 ;
      !A(2,1) = 0 ; A(2,2) = 0 ; A(2,3) = -1 ;
      !A(3,1) =-3 ; A(3,2) = 4 ; A(3,3) = -2 ;
      !det = Determinant(A,ErrStat, ErrMsg)
      !call test_almost_equal('Det of singular 3x3 matrix', real(det,ReKi), 0.0_ReKi, 1e-8_ReKi, .true. , .true.)
      !deallocate(A   )

      ! --- Inverse and determinant of a 3x3 matrix
      allocate(A(3,3)); allocate(Aref(3,3))
      A(1,1) = 7 ; A(1,2) = 2 ; A(1,3) =  1 ;
      A(2,1) = 0 ; A(2,2) = 3 ; A(2,3) = -1 ;
      A(3,1) =-3 ; A(3,2) = 4 ; A(3,3) = -2 ;
      Aref(1,1) =-2 ; Aref(1,2) = 8 ; Aref(1,3) = -5 ;
      Aref(2,1) = 3 ; Aref(2,2) =-11; Aref(2,3) =  7 ;
      Aref(3,1) = 9 ; Aref(3,2) =-34; Aref(3,3) =  21;
      call PseudoInverse(A, Ainv, ErrStat, ErrMsg)
      ! Determinant test
      det = Determinant(A,ErrStat, ErrMsg)
      call test_almost_equal('Det of 3x3 matrix', real(det,ReKi), 1.0_ReKi, 1e-8_ReKi, .true. , .true.)
      det = Determinant(Ainv,ErrStat, ErrMsg)
      call test_almost_equal('Det of 3x3 matrix', real(det,ReKi), 1.0_ReKi, 1e-8_ReKi, .true. , .true.)
      ! Inverse test
      call test_almost_equal('Inverse of 3x3 matrix', real(Aref,ReKi), real(Ainv,ReKi), 1e-8_ReKi, .true., .true.)
      deallocate(A   )
      deallocate(Ainv)
      deallocate(Aref)

      ! --- Inverse of a 3x6 matrix
      allocate(A(3,6))
      allocate(Aref(6,3))
      A(1,:) =  (/  0,   1,   2,   0,   1,   2 /)
      A(2,:) =  (/ -1,   1,   2,  -0,   0,   0 /)
      A(3,:) =  (/ -0,   0,   0,  -1,   1,   2 /)
      Aref(:,:) = transpose(reshape( (/ 0.500000,  -0.583333,  -0.416667, 0.100000,   0.083333,  -0.083333 , 0.200000,   0.166667,  -0.166667 , 0.500000,  -0.416667,  -0.583333 , 0.100000,  -0.083333,   0.083333 , 0.200000,  -0.166667,   0.166667 /), (/ 3, 6 /)))
      call PseudoInverse(A, Ainv, ErrStat, ErrMsg)
      call test_almost_equal('Inverse of 3x6 matrix', real(Aref,ReKi), real(Ainv,ReKi), 1e-6_ReKi, .true., .true.)
      deallocate(A   )
      deallocate(Ainv)
      deallocate(Aref)

      ! --- Inverse of a 6x3 matrix
      allocate(A(6,3))
      allocate(Aref(3,6))
      A(:,1) =  (/  0,   1,   2,   0,   1,   2 /)
      A(:,2) =  (/ -1,   1,   2,  -0,   0,   0 /)
      A(:,3) =  (/ -0,   0,   0,  -1,   1,   2 /)
      Aref(:,:) = reshape( (/ 0.500000,  -0.583333,  -0.416667, 0.100000,   0.083333,  -0.083333 , 0.200000,   0.166667,  -0.166667 , 0.500000,  -0.416667,  -0.583333 , 0.100000,  -0.083333,   0.083333 , 0.200000,  -0.166667,   0.166667 /), (/ 3, 6 /))
      call PseudoInverse(A, Ainv, ErrStat, ErrMsg)
      call test_almost_equal('Inverse of 6x3 matrix', real(Aref,ReKi), real(Ainv,ReKi), 1e-6_ReKi, .true., .true.)
   end subroutine  Test_Linalg

   !> Series of tests for integer lists
   subroutine Test_lists(ErrStat,ErrMsg)
      integer(IntKi)      , intent(out) :: ErrStat
      character(ErrMsgLen), intent(out) :: ErrMsg
      type(IList) :: L1
      type(IList) :: L2
      integer(IntKi) :: e
      ErrStat = ErrID_None
      ErrMsg  = ""
      testname='Lists'

      call init_list(L1, 0, 0 ,ErrStat, ErrMsg)
      call init_list(L2, 10, 12,ErrStat, ErrMsg)

      ! test len
      call test_equal('length',0 ,len(L1))
      call test_equal('length',10,len(L2))

      ! test append
      call append(L1, 5, ErrStat, ErrMsg)
      call append(L1, 3, ErrStat, ErrMsg)
      call append(L1, 1, ErrStat, ErrMsg)
      call test_equal('appnd',L1%List, (/5,3,1/) , .true. , .false.)

      ! test get
      call test_equal('get  ',get(L1,2, ErrStat, ErrMsg), 3)
      e = get(L1,0, ErrStat, ErrMsg)
      call test_equal('get <0 ', ErrStat, ErrID_Fatal)
      e = get(L1,7, ErrStat, ErrMsg)
      call test_equal('get >n ', ErrStat, ErrID_Fatal)

      ! test sort
      call sort(L1, ErrStat, ErrMsg)
      call test_equal('sort ',L1%List, (/1,3,5/) , .true. , .false.)

      ! test reverse
      call reverse(L1, ErrStat, ErrMsg)
      call test_equal('rev  ',L1%List, (/5,3,1/) , .true. , .false.)

      ! test pop
      e = pop(L1, ErrStat, ErrMsg)
      call test_equal('pop  ',e , 1)
      e = pop(L1, ErrStat, ErrMsg)
      call test_equal('pop  ',e , 3)
      e = pop(L1, ErrStat, ErrMsg)
      call test_equal('pop  ',e , 5)
      call destroy_list(L1, ErrStat, ErrMsg)

      ! test unique
      call init_list(L1,(/5,3,2,3,8/),ErrStat, ErrMsg)
      call unique(L1, ErrStat, ErrMsg)
      call test_equal('uniq ',L1%List, (/5,3,2,8/) , .true. , .false.)

      call destroy_list(L1, ErrStat, ErrMsg)
      call destroy_list(L2, ErrStat, ErrMsg)
   end subroutine Test_lists

   !> Test CheckBoard (from FEM), useful for joint stiffness
   subroutine Test_ChessBoard(ErrStat,ErrMsg)
      integer(IntKi)      , intent(out) :: ErrStat
      character(ErrMsgLen), intent(out) :: ErrMsg
      real(ReKi), dimension(:,:), allocatable :: M, Mref
      ErrStat = ErrID_None
      ErrMsg  = ""
      testname='ChessBoard'
      allocate(M(1:6,1:6), Mref(1:6,1:6))
      ! Typical example for pin joint Stiffness add
      Mref(1 , :)= (/4  , -1 , -1 , -1 , -1, -1/)
      Mref(2 , :)= (/-1 ,  4 , -1 , -1 , -1, -1/)
      Mref(3 , :)= (/-1 , -1 ,  4 , -1 , -1, -1/)
      Mref(4 , :)= (/-1 , -1 , -1 ,  4 , -1, -1/)
      Mref(5 , :)= (/-1 , -1 , -1 , -1 ,  4, -1/)
      Mref(6 , :)= (/-1 , -1 , -1 , -1 , -1,  4/)
      call ChessBoard(M, -1._ReKi, -10._ReKi, nSpace=0, diagVal=4._ReKi)
      call test_almost_equal('ChessBoardPin', Mref, M, 1e-6_ReKi, .true., .true.)

      ! Typical example for universal joint Stiffness add
      Mref=0.0_ReKi
      Mref(1 , :)= (/2  ,  0 , -1 ,  0 , -1,  0 /)
      Mref(2 , :)= (/0  ,  2 ,  0 , -1 ,  0, -1 /)
      Mref(3 , :)= (/-1 ,  0 ,  2 ,  0 , -1,  0 /)
      Mref(4 , :)= (/ 0 , -1 ,  0 ,  2 ,  0, -1 /)
      Mref(5 , :)= (/-1 ,  0 , -1 ,  0 ,  2,  0 /)
      Mref(6 , :)= (/ 0 , -1 ,  0 , -1 ,  0,  2 /)
      call ChessBoard(M, -1._ReKi, 0._ReKi, nSpace=1, diagVal=2._ReKi)
      call test_almost_equal('ChessBoardUnv', Mref, M, 1e-6_ReKi, .true., .true.)

      ! Typical example for ball joint Stiffness add
      Mref(1 , :)= (/ 1 ,  0 ,  0 , -1 ,  0,  0 /)
      Mref(2 , :)= (/ 0 ,  1 ,  0 ,  0 , -1,  0 /)
      Mref(3 , :)= (/ 0 ,  0 ,  1 ,  0 ,  0, -1 /)
      Mref(4 , :)= (/-1 ,  0 ,  0 ,  1 ,  0,  0 /)
      Mref(5 , :)= (/ 0 , -1 ,  0 ,  0 ,  1,  0 /)
      Mref(6 , :)= (/ 0 ,  0 , -1 ,  0 ,  0,  1 /)
      call ChessBoard(M, -1._ReKi, 0._ReKi, nSpace=2, diagVal=1._ReKi)
      call test_almost_equal('ChessBoardBll', Mref, M, 1e-6_ReKi, .true., .true.)

      deallocate(M,Mref)
   end subroutine Test_ChessBoard

   subroutine SD_Tests(ErrStat,ErrMsg)
      integer(IntKi)      , intent(out) :: ErrStat !< Error status of the operation
      character(ErrMsgLen), intent(out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
      integer(IntKi)       :: ErrStat2
      character(ErrMsgLen) :: ErrMsg2
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""

      call Test_lists(ErrStat2, ErrMsg2); if(Failed()) return
      call Test_Transformations(ErrStat2, ErrMsg2); if(Failed()) return
      call Test_Linalg(ErrStat2, ErrMsg2); if(Failed()) return
      call Test_ChessBoard(ErrStat2, ErrMsg2); if(Failed()) return
      contains
         logical function Failed()
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_Tests') 
            Failed =  ErrStat >= AbortErrLev
         end function failed
   end subroutine SD_Tests


end module SubDyn_Tests
