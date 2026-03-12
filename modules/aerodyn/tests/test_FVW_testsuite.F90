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

   ! --- Run all tests at once
   call FVW_RunTests                  (errStat,errMsg); 
   @assertEqual(0, errStat, 'All FVW tests                 ')

   ! --- Run individual tests
   !call Test_LinSolve                 (errStat,errMsg); 
   !@assertEqual(0, errStat, 'Test_LinSolve                 ')
   !call Test_SrcPnl_Sphere            (errStat,errMsg); 
   !@assertEqual(0, errStat, 'Test_SrcPnl_Sphere            ')
   !call Test_BiotSavart_SrcPnl        (errStat,errMsg); 
   !@assertEqual(0, errStat, 'Test_BiotSavart_SrcPnl        ')
   !call Test_BiotSavart_Sgmt          (errStat,errMsg); 
   !@assertEqual(0, errStat, 'Test_BiotSavart_Sgmt          ')
   !call Test_BiotSavart_Part          (errStat,errMsg); 
   !@assertEqual(0, errStat, 'Test_BiotSavart_Part          ')
   !call Test_BiotSavart_PartTree      (errStat,errMsg); 
   !@assertEqual(0, errStat, 'Test_BiotSavart_PartTree      ')
   !call Test_SegmentsToPart           (errStat,errMsg); 
   !@assertEqual(0, errStat, 'Test_SegmentsToPart           ')
   !call FVW_Test_WakeInducedVelocities(errStat,errMsg); 
   !@assertEqual(0, errStat, 'FVW_Test_WakeInducedVelocities')

end subroutine test_AD_FVW
