module SubDyn_Tests
   use NWTC_Library
   USE SubDyn_Types
   USE SD_FEM
   use IntegerList
   
   public :: SD_Tests
   private

   character(len=255),save :: testname
   interface test_equal; module procedure &
         test_equal_i1, &
         test_equal_i0
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
            STOP -1 !OTHER-COMPILER
            STOP ! COMPAQ-COMPILER
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
                STOP -1 !OTHER-COMPILER
	         	STOP ! COMPAQ-COMPILER
            else
                call test_success(InfoAbs)
            endif
        endif
    end subroutine




    ! --------------------------------------------------------------------------------}
    ! --- Specific SubDyn tests 
    ! --------------------------------------------------------------------------------{
   subroutine Test_CB_Results(MBBt, MBMt, KBBt, OmegaM, DOFTP, DOFM, ErrStat, ErrMsg,Init,p)
      TYPE(SD_InitType),      INTENT(  in)                :: Init         ! Input data for initialization routine
      TYPE(SD_ParameterType), INTENT(inout)                :: p           ! Parameters
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


   subroutine Test_lists(ErrStat,ErrMsg)
      integer(IntKi)      , intent(out) :: ErrStat
      character(ErrMsgLen), intent(out) :: ErrMsg
      type(IList) :: L1
      type(IList) :: L2
      integer(IntKi) :: e
      ErrStat = ErrID_None
      ErrMsg  = ""
      testname='Lists'

      call init(L1, 0, 0 ,ErrStat, ErrMsg)
      call init(L2, 10, 12,ErrStat, ErrMsg)

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

   end subroutine

   subroutine SD_Tests(ErrStat,ErrMsg)
      integer(IntKi)      , intent(out) :: ErrStat !< Error status of the operation
      character(ErrMsgLen), intent(out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
      integer(IntKi)       :: ErrStat2
      character(ErrMsgLen) :: ErrMsg2
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""

      call Test_lists(ErrStat2, ErrMsg2)
   end subroutine SD_Tests


end module SubDyn_Tests
