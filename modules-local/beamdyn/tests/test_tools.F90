module test_tools
  
implicit none
  
contains  
    
    subroutine calcWMParameters(params, angle, n)
        use BeamDyn_Subs
        implicit none
        
        real(BDKi), intent(  out), dimension(3) :: params
        real(BDKi), intent(in   )               :: angle
        real(BDKi), intent(in   ), dimension(3) :: n
        
        params = 4.0 * tan(angle/4.0) * n
    end subroutine
    
    function calcRotationMatrix(angle, axis)
        use BeamDyn_Subs
        implicit none
        
        real(BDKi),             dimension(3,3) :: calcRotationMatrix
        real(BDKi), intent(in)                 :: angle
        real(BDKi), intent(in), dimension(3)   :: axis
        real(BDKi),             dimension(3,3) :: r

        r(1,:) = (/              cos(angle) + (1-cos(angle))*axis(1)**2, axis(1)*axis(2)*(1-cos(angle)) - axis(3)*sin(angle), axis(1)*axis(3)*(1-cos(angle)) + axis(2)*sin(angle) /)
        r(2,:) = (/ axis(2)*axis(1)*(1-cos(angle)) + axis(3)*sin(angle),              cos(angle) + (1-cos(angle))*axis(2)**2, axis(2)*axis(3)*(1-cos(angle)) - axis(1)*sin(angle) /)
        r(3,:) = (/ axis(3)*axis(1)*(1-cos(angle)) - axis(2)*sin(angle), axis(3)*axis(2)*(1-cos(angle)) + axis(1)*sin(angle),              cos(angle) + (1-cos(angle))*axis(3)**2 /)
        
        calcRotationMatrix = r
        
    end function

    function RonXAxis(angle)
        use BeamDyn_Subs
        implicit none
        
        real(BDKi) :: angle, r(3,3), RonXAxis(3,3)
        
        r(1,:) = (/  1.0,         0.0,         0.0 /)
        r(2,:) = (/  0.0,  cos(angle), -sin(angle) /)
        r(3,:) = (/  0.0,  sin(angle),  cos(angle) /)
        RonXAxis = r
    end function  
        
    function getMassMatrix()
        use BeamDyn_Subs
        implicit none
        
        real(BDKi), dimension(6,6) :: getMassMatrix
        getMassMatrix(1,:) =  (/   1.E0,    0.0,    0.0,    0.0,    0.0,   -0.5 /)
        getMassMatrix(2,:) =  (/    0.0,   1.E0,    0.0,    0.0,    0.0,    0.5 /)
        getMassMatrix(3,:) =  (/    0.0,    0.0,   1.E0,    0.5,   -0.5,    0.0 /)
        getMassMatrix(4,:) =  (/    0.0,    0.0,    0.5,   1.E0,   -1.0,    0.0 /)
        getMassMatrix(5,:) =  (/    0.0,    0.0,   -0.5,  -1.E0,   1.E0,    0.0 /)
        getMassMatrix(6,:) =  (/   -0.5,    0.5,    0.0,    0.0,    0.0,   2.E0 /)
    end function
    
    function getStiffnessMatrix()
        use BeamDyn_Subs
        implicit none
        
        real(BDKi), dimension(6,6) :: getStiffnessMatrix
        getStiffnessMatrix(1,:) = (/   1.E4,    0.0,    0.0,    0.0,    0.0,    0.0 /)
        getStiffnessMatrix(2,:) = (/    0.0,   1.E4,    0.0,    0.0,    0.0,    0.0 /)
        getStiffnessMatrix(3,:) = (/    0.0,    0.0,   1.E4,    0.0,    0.0,    0.0 /)
        getStiffnessMatrix(4,:) = (/    0.0,    0.0,    0.0,   1.E2,    0.0,    0.0 /)
        getStiffnessMatrix(5,:) = (/    0.0,    0.0,    0.0,    0.0,   1.E2,    0.0 /)
        getStiffnessMatrix(6,:) = (/    0.0,    0.0,    0.0,    0.0,    0.0, 200.E0 /)
    end function
    
end module