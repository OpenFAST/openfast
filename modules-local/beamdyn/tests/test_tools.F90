module test_tools
  
use BeamDyn_Types

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
    
    function identity()
        use BeamDyn_Subs
        implicit none
        
        real(BDKi) :: identity(3,3)
        
        identity(1,:) = (/ 1.0, 0.0, 0.0 /)
        identity(2,:) = (/ 0.0, 1.0, 0.0 /)
        identity(3,:) = (/ 0.0, 0.0, 1.0 /)
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
        getMassMatrix(4,:) =  (/    0.0,    0.0,    0.5,   1.E0,    0.0,    0.0 /)
        getMassMatrix(5,:) =  (/    0.0,    0.0,   -0.5,    0.0,   1.E0,    0.0 /)
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

    function getGravityInZ()
        use BeamDyn_Subs
        implicit none
        
        real(BDKi), dimension(3) :: getGravityInZ
        getGravityInZ = (/ 0.0, 0.0, -9.806 /)
    end function
    
    type(BD_ParameterType) function simpleParameterType()
        
        type(BD_ParameterType) :: p
        integer                :: i, j
        integer                :: ErrStat
        character(1024)        :: ErrMsg
        
        ! scalars
        p%elem_total = 1
        p%nodes_per_elem = 16
        p%nqp = 16
        p%qp_indx_offset = 0
        
        ! fixed size arrays
        p%Glb_crv = (/ 0.0, 0.0, 0.0 /)
        p%GlbRot = identity()
        
        ! allocate arrays
        call AllocAry(p%qp%mmm, p%nqp, p%elem_total, 'qp_mmm', ErrStat, ErrMsg)
        call AllocAry(p%qp%mEta, 3, p%nqp, p%elem_total, 'qp_RR0mEta', ErrStat, ErrMsg)
        call AllocAry(p%Mass0_QP, 6, 6, p%nqp*p%elem_total, 'Mass0_QP', ErrStat, ErrMsg)
        
        ! construct arrays
        p%qp%mmm = getMassMatrix()
        
        do j=1, p%elem_total
            do i=1, p%nqp
                p%qp%mEta(:,i,j) = (/ 0.0, 0.0, 0.0 /)
                p%Mass0_QP(:,:,(i-1)*p%elem_total+j) = getMassMatrix()
            end do
        end do
        
        ! set the return value
        simpleParameterType = p
        
    end function
    
    type(BD_MiscVarType) function simpleMiscVarType(nqp, nelem)
        
        type(BD_MiscVarType) :: m
        integer, intent(in)  :: nqp, nelem
        integer              :: i, j
        integer              :: ErrStat
        character(1024)      :: ErrMsg
        
        ! scalars
        
        ! fixed size arrays
        
        ! allocate arrays
        call AllocAry(m%qp%Fg, 6, nqp, nelem, 'qp_Fg', ErrStat, ErrMsg)
        call AllocAry(m%qp%RR0, 3, 3, nqp, nelem, 'qp_RR0', ErrStat, ErrMsg)
        call AllocAry(m%qp%RR0mEta, 3, nqp, nelem, 'qp_RR0mEta', ErrStat, ErrMsg)
        call AllocAry(m%DistrLoad_QP, 6, nqp, nelem, 'DistrLoad_QP', ErrStat, ErrMsg)
        call AllocAry(m%qp%rho, 3, 3, nqp, nelem, 'qp_rho', ErrStat, ErrMsg)
        
        ! construct arrays
        do j=1, nelem
            do i=1, nqp
                m%qp%RR0(:,:,i,j) = identity()
                m%qp%RR0mEta(:,i,j) = (/ 0.0, 0.0, 0.0 /)
            end do
        end do
        
        ! set the return value
        simpleMiscVarType = m
        
    end function

    type(BD_InputType) function simpleInputType(nqp, nelem)
        
        type(BD_InputType)   :: i
        integer, intent(in)  :: nqp, nelem
        integer              :: j
        integer              :: ErrStat
        character(1024)      :: ErrMsg
        
        ! scalars
        
        ! fixed size arrays
        
        ! allocate arrays
        call AllocAry(i%DistrLoad%Force, 3, nqp*nelem, 'DistrLoadForce', ErrStat, ErrMsg)
        call AllocAry(i%DistrLoad%Moment, 3, nqp*nelem, 'DistrLoadMoment', ErrStat, ErrMsg)
        
        ! construct arrays
        do j = 1, nqp*nelem
            i%DistrLoad%Force(:,j)  = (/  3*(j-1)+1,  3*(j-1)+2,  3*(j-1)+3 /)
            i%DistrLoad%Moment(:,j) = (/ -3*(j-1)-1, -3*(j-1)-2, -3*(j-1)-3 /)
        end do
        
        ! set the return value
        simpleInputType = i
        
    end function
    
end module
