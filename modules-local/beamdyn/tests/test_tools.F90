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

    subroutine calcRotationMatrix(r, angle, axis)
        use BeamDyn_Subs
        implicit none
        
        real(BDKi), intent(  out), dimension(3,3) :: r
        real(BDKi), intent(in   )                 :: angle
        integer,    intent(in   )                 :: axis
        
        if (axis == 1) then
            r(:,1) = (/ 1.0,        0.0,         0.0 /)
            r(:,2) = (/ 0.0, cos(angle), -sin(angle) /)
            r(:,3) = (/ 0.0, sin(angle),  cos(angle) /)
        else if (axis == 2) then
            r(:,1) = (/ cos(angle), 0.0, -sin(angle) /)
            r(:,2) = (/       0.0,  1.0,         0.0 /)
            r(:,3) = (/ sin(angle), 0.0,  cos(angle) /)
        else if (axis ==3) then
            r(:,1) = (/ cos(angle), -sin(angle), 0.0 /)
            r(:,2) = (/ sin(angle),  cos(angle), 0.0 /)
            r(:,3) = (/        0.0,         0.0, 1.0 /)
        end if
    end subroutine
end module