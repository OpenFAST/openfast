@test
subroutine test_AD_FVW()
    ! test branches
    ! - known valid checks for various FVW routines (contained in own module)
    ! - known invalid rotation matrix: halve the angle of the diagonal elements

    use pFUnit_mod
    use NWTC_Num
    use FVW_Tests 
    
    implicit none
    
    integer(IntKi)       :: ErrStat
    character(ErrMsgLen) :: ErrMsg
    character(1024)      :: testname

    ! initialize NWTC_Num constants
    call SetConstants()

!This is a single routine that contains the test cases below.   
   ! --------------------------------------------------------------------------    
   testname = "Set of FVW tests"
   call FVW_RunTests( ErrStat, ErrMsg )
   @assertEqual(0, ErrStat, testname)


! test routines from FVW_RunTests to be run individually -- except these are all private
!   ! --------------------------------------------------------------------------    
!   testname = "known valid Biot-Savart segment"
!   call Test_BiotSavart_Sgmt(testname, ErrStat, ErrMsg)
!   @assertEqual(0, ErrStat, testname)
!
!   ! --------------------------------------------------------------------------    
!   testname = "known valid Biot-Savart part"
!   call Test_BiotSavart_Part(testname, ErrStat, ErrMsg)
!   @assertEqual(0, ErrStat, testname)
!
!   ! --------------------------------------------------------------------------    
!   testname = "known valid Biot-Savart to part-tree"
!   call Test_BiotSavart_PartTree(testname, ErrStat, ErrMsg)
!   @assertEqual(0, ErrStat, testname)
!
!   ! --------------------------------------------------------------------------    
!   testname = "known valid segment split to parts"
!   call Test_SegmentsToPart(testname, ErrStat, ErrMsg)
!   @assertEqual(0, ErrStat, testname)
 
end subroutine test_AD_FVW
