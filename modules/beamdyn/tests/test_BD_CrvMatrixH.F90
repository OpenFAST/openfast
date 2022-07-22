@test
subroutine test_BD_CrvMatrixH()
    ! test branches
    ! - simple rotation with known parameters: Pi on xaxis
    ! - 0 rotation
    ! - small rotation with baseline WM parameters calculated
    
    ! TODO
    ! invalid wm parameters (if thats a thing)
    ! does the implemented WM formulation have any boundaries?
    
    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools
    
    implicit none
        
    real(BDKi), dimension(3,3) :: testH, baselineH
    real(BDKi), dimension(3)   :: wmparams
    real(BDKi)                 :: angle, n(3)
    character(1024)            :: testname
    real(BDKi)                 :: tolerance
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    ! set the rotation axis for all tests
    n = (/ 1., 0., 0. /) ! x axis
    
    
    ! --------------------------------------------------------------------------
    testname = "simple rotation with known parameters: Pi on xaxis:"
    angle = Pi_D

    ! Wiener-Milenkovic parameters are <4.0, 0.0, 0.0>
    wmparams = (/ 4.0, 0.0, 0.0 /)
        
    baselineH = H(wmparams)
    
    call BD_CrvMatrixH(wmparams, testH)
    
    @assertEqual(baselineH, testH, tolerance, testname)
    
    
    ! --------------------------------------------------------------------------
    testname = "0 rotation:"
    angle = 0
    
    ! Wiener-Milenkovic parameters are <0.0, 0.0, 0.0>
    wmparams = (/ 0.0, 0.0, 0.0 /)
        
    baselineH = H(wmparams)
    
    call BD_CrvMatrixH(wmparams, testH)
    
    @assertEqual(baselineH, testH, tolerance, testname)
    
    
    ! --------------------------------------------------------------------------
    testname = "small rotation with baseline WM parameters calculated:"
    angle = 0.1*Pi_D
    
    ! Wiener-Milenkovic parameters are calculated; note tangent is asymptotic at +/- pi/2
    wmparams = 4.*tan(angle/4.)*n
        
    baselineH = H(wmparams)
    
    call BD_CrvMatrixH(wmparams, testH)
    
    @assertEqual(baselineH, testH, tolerance, testname)
    
    contains
        function H(c)
            real(BDKi)          :: c0, c(3)
            real(BDKi)          :: H(3,3)
            
            c0 = 2.0 - dot_product(c,c) / 8.0
            
            H(1,:) = (/ c0 + c(1)*c(1)/4., c(1)*c(2)/4. - c(3), c(1)*c(3)/4. + c(2) /)
            H(2,:) = (/ c(1)*c(2)/4. + c(3), c0 + c(2)*c(2)/4., c(2)*c(3)/4. - c(1) /)
            H(3,:) = (/ c(1)*c(3)/4. - c(2), c(2)*c(3)/4. + c(1), c0 + c(3)*c(3)/4. /)
            H = 2.*H/(4.-c0)**2
            
        end function
end subroutine
