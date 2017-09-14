@test
subroutine test_ExtractRelativeRotation()
    ! this is actually an integration test not a unit test...

    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools
    
    implicit none
    
    real(BDKi), dimension(3,3) :: rotationmatrix1, rotationmatrix2, composedrotationmatrix
    real(BDKi), dimension(3)   :: wmparams, result_wmparams, rotationaxis, baseline_wmparams
    real(BDKi)                 :: rotationangle
    character(1024)            :: testname
    real(BDKi)                 :: tolerance
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg
    
    type(BD_ParameterType)     :: parametertype

    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
        
    testname = "extract xaxis by pi/2, zaxis by pi"
        
    ! input rotation matrix is rotation of pi/2 around x-axis
    ! 1.0  0.0  0.0
    ! 0.0  0.0 -1.0
    ! 0.0  1.0  0.0
    call calcRotationMatrix(rotationmatrix1, PiBy2, 1)
    
    ! input parameterized rotation is of pi around z-axis
    ! <0.0, 0.0, 4.0>
    call calcWMParameters(wmparams, Pi, (/ 0.0, 0.0, 1.0 /))
    
    ! call the test function
    call ExtractRelativeRotation(rotationmatrix1, parametertype, result_wmparams, ErrStat, ErrMsg)
    
    call calcRotationMatrix(rotationmatrix2, Pi, 3)
    composedrotationmatrix = matmul(rotationmatrix1, rotationmatrix2)
    
    ! get the axis and angle of the composed rotation
    ! https://en.wikipedia.org/wiki/Rotation_matrix#Determining_the_axis  
    rotationaxis = (/                                                          &
                      composedrotationmatrix(2,3)-composedrotationmatrix(3,2), &
                      composedrotationmatrix(3,1)-composedrotationmatrix(1,3), &
                      composedrotationmatrix(1,2)-composedrotationmatrix(2,1)  &
                   /)
    rotationaxis = rotationaxis/maxval(rotationaxis)
    rotationangle = acos((composedrotationmatrix(1,1) + composedrotationmatrix(2,2) + composedrotationmatrix(3,3) - 1)/2)
    call calcWMParameters(baseline_wmparams, rotationangle, rotationaxis)
    
    @assertEqual(baseline_wmparams, result_wmparams, tolerance, testname)

end subroutine
