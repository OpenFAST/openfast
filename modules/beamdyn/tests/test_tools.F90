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
        
        r(1,:) = (/  1.0_BDKi,    0.0_BDKi,    0.0_BDKi /)
        r(2,:) = (/  0.0_BDKi,  cos(angle), -sin(angle) /)
        r(3,:) = (/  0.0_BDKi,  sin(angle),  cos(angle) /)
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
    
    type(BD_ParameterType) function simpleParameterType(elem_total, nodes_per_elem, nqp, qp_indx_offset, refine) RESULT(p)
  
        integer, intent(in   ) :: elem_total
        integer, intent(in   ) :: nodes_per_elem
        integer, intent(in   ) :: nqp
        integer, intent(in   ) :: qp_indx_offset
        integer, intent(in   ) :: refine

        integer                :: i, j
        integer                :: ErrStat
        character(1024)        :: ErrMsg

        p%elem_total = elem_total
        p%nodes_per_elem = nodes_per_elem
        p%nqp = nqp
        p%qp_indx_offset = qp_indx_offset
        p%refine = refine

        p%dof_node = 6

        ! fixed size arrays
        p%Glb_crv = (/ 0.0, 0.0, 0.0 /)
        p%GlbRot = identity()

        ! allocate arrays
        call AllocAry(p%qp%mmm, p%nqp, p%elem_total, 'qp_mmm', ErrStat, ErrMsg)
        call AllocAry(p%qp%mEta, 3, p%nqp, p%elem_total, 'qp_RR0mEta', ErrStat, ErrMsg)
        call AllocAry(p%Mass0_QP, 6, 6, p%nqp*p%elem_total, 'Mass0_QP', ErrStat, ErrMsg)  ! if called, this allocated in InitializeMassStiffnessMatrices
        call AllocAry(p%Stif0_QP, 6, 6, p%nqp*p%elem_total, 'Stif0_QP', ErrStat, ErrMsg)  ! if called, this allocated in InitializeMassStiffnessMatrices
        call AllocAry(p%QPtw_Shp_Shp_Jac, p%nqp, p%nodes_per_elem, p%nodes_per_elem, p%elem_total, 'QPtw_Shp_Shp_Jac', ErrStat, ErrMsg)
        call AllocAry(p%QPtw_ShpDer_ShpDer_Jac, p%nqp, p%nodes_per_elem, p%nodes_per_elem, p%elem_total, 'QPtw_ShpDer_ShpDer_Jac', ErrStat, ErrMsg)
        call AllocAry(p%QPtw_Shp_ShpDer, p%nqp, p%nodes_per_elem, p%nodes_per_elem, 'QPtw_Shp_ShpDer', ErrStat, ErrMsg)
        call AllocAry(p%QPtw_Shp_Jac, p%nqp, p%nodes_per_elem, p%elem_total, 'QPtw_Shp_Jac', ErrStat, ErrMsg)
        call AllocAry(p%Shp, p%nodes_per_elem, p%nqp, 'Shp', ErrStat, ErrMsg)
        call AllocAry(p%ShpDer, p%nodes_per_elem, p%nqp, 'ShpDer', ErrStat, ErrMsg)
        call AllocAry(p%QPtN, p%nqp, 'QPtN', ErrStat, ErrMsg)
        call AllocAry(p%QPtWeight, p%nqp, 'QPtWeight', ErrStat, ErrMsg)
        call AllocAry(p%QPtw_ShpDer, p%nqp, p%nodes_per_elem, 'QPtw_ShpDer', ErrStat, ErrMsg)
        call AllocAry(p%Jacobian, p%nqp, p%elem_total, 'Jacobian', ErrStat, ErrMsg)
        call AllocAry(p%uuN0, p%dof_node, p%nodes_per_elem, p%elem_total,'uuN0', ErrStat, ErrMsg)

        call AllocAry(p%uu0, p%dof_node, p%nqp,   p%elem_total,'uu0', ErrStat, ErrMsg)
        call AllocAry(p%E10, p%dof_node/2, p%nqp, p%elem_total,'E10', ErrStat, ErrMsg)
        call AllocAry(p%rrN0, p%dof_node/2, p%nodes_per_elem, p%elem_total,'rrN0', ErrStat, ErrMsg)

        CALL AllocAry(p%node_elem_idx,p%elem_total,2,'start and end node numbers of elements in p%node_total sized arrays',ErrStat,ErrMsg)

        ! construct arrays
        p%qp%mmm = 1.0

        DO i=1,p%elem_total
           p%node_elem_idx(i,1) =  (i-1)*(p%nodes_per_elem-1) + 1           ! First node in element
           p%node_elem_idx(i,2) =   i   *(p%nodes_per_elem-1) + 1           ! Last node in element
        ENDDO

        do j=1, p%elem_total
            do i=1, p%nqp
                p%qp%mEta(:,i,j) = (/ 0.0, 0.0, 0.0 /)
                p%Mass0_QP(:,:,(i-1)*p%elem_total+j) = getMassMatrix()
            end do
        end do

    end function
   
    type(BD_MiscVarType) function simpleMiscVarType(nqp, dof_node, elem_total, nodes_per_elem) RESULT(m)
        
        integer, intent(in)  :: nqp, elem_total, dof_node, nodes_per_elem
        integer              :: i, j
        integer              :: ErrStat
        character(1024)      :: ErrMsg
        
        ! scalars
        
        ! fixed size arrays
        
        ! allocate arrays
        call AllocAry(m%qp%Fg, 6, nqp, elem_total, 'qp_Fg', ErrStat, ErrMsg)
        call AllocAry(m%qp%RR0, 3, 3, nqp, elem_total, 'qp_RR0', ErrStat, ErrMsg)
        call AllocAry(m%qp%RR0mEta, 3, nqp, elem_total, 'qp_RR0mEta', ErrStat, ErrMsg)
        call AllocAry(m%DistrLoad_QP, 6, nqp, elem_total, 'DistrLoad_QP', ErrStat, ErrMsg)

        CALL AllocAry(m%qp%uuu, dof_node  ,nqp,elem_total, 'm%qp%uuu displacement at quadrature point',ErrStat,ErrMsg)
        CALL AllocAry(m%qp%uup, dof_node/2,nqp,elem_total, 'm%qp%uup displacement prime at quadrature point',ErrStat,ErrMsg)

         ! E1, kappa -- used in force calculations
        CALL AllocAry(m%qp%E1,    dof_node/2,nqp,elem_total, 'm%qp%E1    at quadrature point',ErrStat,ErrMsg)
        CALL AllocAry(m%qp%kappa, dof_node/2,nqp,elem_total, 'm%qp%kappa at quadrature point',ErrStat,ErrMsg)
        CALL AllocAry(m%qp%RR0,   3,3,       nqp,elem_total, 'm%qp%RR0 at quadrature point',ErrStat,ErrMsg)
        CALL AllocAry(m%qp%Stif,  6,6,       nqp,elem_total, 'm%qp%Stif at quadrature point',ErrStat,ErrMsg)

        CALL AllocAry(m%qp%RR0mEta, dof_node/2, nqp, elem_total, 'm%qp%RRo times p%qp%mEta at quadrature point',ErrStat,ErrMsg)
        call AllocAry(m%qp%rho, 3, 3, nqp, elem_total, 'qp_rho', ErrStat, ErrMsg)
        CALL AllocAry(m%qp%betaC, 6,6, nqp,elem_total, 'm%qp%betaC at quadrature point',ErrStat,ErrMsg)
       
        CALL AllocAry(m%qp%Fc, dof_node, nqp, elem_total, 'm%qp%Fc at quadrature point',ErrStat,ErrMsg)
        CALL AllocAry(m%qp%Fd, dof_node, nqp, elem_total, 'm%qp%Fd at quadrature point',ErrStat,ErrMsg)

        CALL AllocAry(m%Nrrr, dof_node/2, nodes_per_elem, elem_total,'Nrrr: rotation parameters relative to root',  ErrStat,ErrMsg)

        ! construct arrays
        do j=1, elem_total
            do i=1, nqp
                m%qp%RR0(:,:,i,j) = identity()
                m%qp%RR0mEta(:,i,j) = (/ 0.0, 0.0, 0.0 /)
            end do
        end do
        
    end function

    type(BD_InputType) function simpleInputType(nqp, nelem) RESULT(i)
        
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
        
    end function
    
    type(BD_InputFile) function simpleInputFile() RESULT(i)
        
        integer              :: j
        integer              :: ErrStat
        character(1024)      :: ErrMsg
        
        ! scalars
        i%QuasiStaticInit = .false.   ! -  - - "QuasiStaticInit" -
        i%member_total = 1    ! -  - - "Total number of members" -
        i%kp_total = 3        ! -  - - "Total number of key point" -
        i%order_elem = 15     ! -  - - "Order of interpolation (basis) function" -
        i%NRMax = 10          ! -  - - "Max number of iterations in Newton Ralphson algorithm" -
        i%quadrature = 1      ! -  - - "Quadrature: 1: Gauss; 2: Trapezoidal" -
        i%n_fact = 5          ! -  - - "Factorization frequency" -
        i%refine = 1          ! -  - - "FE mesh refinement factor for trapezoidal quadrature" -
        i%rhoinf = 0.0        ! -  - - "Numerical damping parameter for generalized-alpha integrator" -
        i%DTBeam = 2E-03      ! -  - - "Time interval for BeamDyn  calculations {or default} (s)" -
        i%UsePitchAct = .FALSE. ! -  - - "Whether to use a pitch actuator inside BeamDyn" (flag) 
        i%pitchJ = 0.0        ! - - -     "Pitch actuator inertia" (kg-m^2)
        i%pitchK = 0.0        ! - - -     "Pitch actuator stiffness" (kg-m^2/s^2) 
        i%pitchC = 0.0        ! - - -     "Pitch actuator damping" - (kg-m^2/s)  
        i%Echo = .TRUE.       ! -  - - "Echo"
        i%NNodeOuts = 1       ! -  - - "Number of node outputs [0 - 9]" -
        i%OutNd = 1           ! {9} - - "Nodes whose values will be output" -
        i%SumPrint = .TRUE.   ! -  - - "Print summary data to file? (.sum)" -
        i%OutFmt = "ES16.8E2" ! -  - - "Format specifier" -
        
        ! fixed size arrays
        i%kp_member = (/ 3 /) !{:} - - "Number of key points in each member" -
        i%OutList = (/ "TipTDxr, TipTDyr, TipTDzr", "TipRDxr, TipRDyr, TipRDzr" /)  ! {:} - - "List of user-requested output channels" -
        
        ! allocate arrays
        call AllocAry(i%kp_coordinate, 3, 4, 'kp_coordinate', ErrStat, ErrMsg)
        
        ! construct arrays
        i%kp_coordinate(1,:) = (/ 0.000000, 0.000000,  0.0000, 0.00000 /)  !  {:}{:} - - "Key point coordinates array" -
        i%kp_coordinate(2,:) = (/ 0.000000, 0.000000,  5.0000, 0.00000 /)
        i%kp_coordinate(3,:) = (/ 0.000000, 0.000000, 10.0000, 0.00000 /)
        
    end function

    type(BD_ContinuousStateType) function simpleContinuousStateType(node_total, nodes_per_elem, elem_total) RESULT(x)
        
        integer, intent(in)  :: node_total,nodes_per_elem, elem_total
        integer              :: j
        integer              :: ErrStat
        character(1024)      :: ErrMsg
        
        ! scalars
        
        ! fixed size arrays
        
        ! allocate arrays
        call AllocAry(x%q,    6, node_total, 'Displacement/Rotation Nodal DOF', ErrStat, ErrMsg)
        call AllocAry(x%dqdt, 6, node_total, 'Velocity Nodal DOF', ErrStat, ErrMsg)
        
    end function
    
end module
