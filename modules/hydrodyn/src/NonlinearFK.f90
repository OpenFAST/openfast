module NonlinearFK

    USE NonlinearFK_Types
    USE SeaSt_WaveField

    implicit none

    ! Quadrature coordinates (barycentric) and weights for 4-point rule
    real(ReKi), parameter :: Qdrt_L(3,4) = reshape([ &
    1.0_ReKi/3.0_ReKi, 1.0_ReKi/3.0_ReKi, 1.0_ReKi/3.0_ReKi, &
    3.0_ReKi/5.0_ReKi, 1.0_ReKi/5.0_ReKi, 1.0_ReKi/5.0_ReKi, &
    1.0_ReKi/5.0_ReKi, 3.0_ReKi/5.0_ReKi, 1.0_ReKi/5.0_ReKi, &
    1.0_ReKi/5.0_ReKi, 1.0_ReKi/5.0_ReKi, 3.0_ReKi/5.0_ReKi  &
    ], [3, 4])

    real(ReKi), parameter :: Qdrt_w(4) = [ &
    -27.0_ReKi/48.0_ReKi, &
     25.0_ReKi/48.0_ReKi, &
     25.0_ReKi/48.0_ReKi, &
     25.0_ReKi/48.0_ReKi ]

    ! Number of Quadrature Points
    integer(IntKi), parameter :: nQdrt = 4

    PUBLIC :: NonlinearFK_Init
    PUBLIC :: NonlinearFK_CalcOutput

contains

subroutine NonlinearFK_Init(InitInp, p, m, InitOut, ErrStat, ErrMsg)
    type(NonlinearFK_InitInputType), intent(in   ) :: InitInp
    type(NonlinearFK_ParameterType), intent(inout) :: p
    type(NonlinearFK_MiscVarType),   intent(inout) :: m
    type(NonlinearFK_InitOutputType),intent(  out) :: InitOut
    integer(IntKi),                  intent(  out) :: ErrStat     !< Error status of the operation
    character(*),                    intent(  out) :: ErrMsg      !< Error message if ErrStat /= ErrID_None

    INTEGER(IntKi)           :: i
    type(STLGeomType)        :: STLGeom
    INTEGER(IntKi)           :: ErrStat2
    CHARACTER(ErrMsgLen)     :: ErrMsg2
    Character(*), parameter  :: RoutineName = 'NonlinearFK_Init'

    ErrStat = ErrID_None
    ErrMsg  = ''

    p%WaveField => InitInp%WaveField
    p%nBody = InitInp%nBody

    CALL AllocAry( p%FKMod, p%nBody, 'FKMod', ErrStat2, ErrMsg2); if (Failed()) return
    p%FKMod = InitInp%FKMod

    allocate( p%Bodies(p%nBody), stat=ErrStat2)
    if (ErrStat2 /= 0) then
        ErrStat = ErrID_Fatal
        ErrMsg  = trim(RoutineName)//": Failed to allocate memory for p%Bodies."
        return
    end if
    allocate( m%Bodies(p%nBody), stat=ErrStat2)
    if (ErrStat2 /= 0) then
        ErrStat = ErrID_Fatal
        ErrMsg  = trim(RoutineName)//": Failed to allocate memory for m%Bodies."
        return
    end if
    allocate( InitOut%Buoyancy(6,p%nBody), stat=ErrStat2)
    if (ErrStat2 /= 0) then
        ErrStat = ErrID_Fatal
        ErrMsg  = trim(RoutineName)//": Failed to allocate memory for InitOut%Buoyancy."
        return
    end if
    InitOut%Buoyancy = 0.0_ReKi

    do i=1,p%nBody
        if (p%FKMod(i) /= FKMod_none) then
            call read_ascii_stl(InitInp%GeoFile(i), STLGeom, ErrStat2, ErrMsg2); if (Failed()) return
            call Body_Init(STLGeom, [InitInp%PtfmRefxt(i),InitInp%PtfmRefyt(i),InitInp%PtfmRefzt(i)], InitInp%PtfmRefztRot(i), &
                           p%Bodies(i), m%Bodies(i), p%WaveField%RhoXg, p%WaveField%MSL2SWL, InitOut%Buoyancy(:,i), ErrStat2, ErrMsg2)
                if (Failed()) return
        end if
    end do

    call cleanup()

contains
   subroutine cleanup()
       if (allocated(STLGeom%tris)) deallocate(STLGeom%tris)
   end subroutine cleanup

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
      Failed =  ErrStat >= AbortErrLev
      if (Failed) call cleanup()
   end function Failed
end subroutine NonlinearFK_Init

subroutine read_ascii_stl(filename, STLGeom, ErrStat, ErrMsg)
    character(len=*),  intent(in   ) :: filename
    type(STLGeomType), intent(inout) :: STLGeom
    integer(IntKi),    intent(  out) :: ErrStat     !< Error status of the operation
    character(*),      intent(  out) :: ErrMsg      !< Error message if ErrStat /= ErrID_None
    
    integer(IntKi)          :: n_tris
    integer(IntKi)          :: io_status, iu, i
    character(len=256)      :: line
    character(len=32)       :: dummy1, dummy2
    character(*), parameter :: RoutineName = "read_ascii_stl"

    ! Initialize variables
    ErrStat = ErrID_None
    ErrMsg  = ""

    if (allocated(STLGeom%tris)) deallocate(STLGeom%tris)

    ! 1. Open the file
    open(newunit=iu, file=filename, status='old', action='read', iostat=io_status)
    if (io_status /= 0) then
        ErrStat = ErrID_Fatal
        ErrMsg  = trim(RoutineName)//": Could not open STL file: "//trim(filename)//". "
        return
    end if

    ! 2. First Pass: Count the number of triangles
    n_tris = 0
    do
        read(iu, '(A)', iostat=io_status) line
        if (io_status /= 0) exit ! End of file
        
        line = adjustl(line)
        if (index(line, 'facet normal') == 1) then
            n_tris = n_tris + 1
        end if
    end do
    
    if (n_tris == 0_IntKi) then
        ErrStat = ErrID_Fatal
        ErrMsg  = trim(RoutineName)//": stl file "//trim(filename)//" is invalid or contains no triangles. "
        close(iu)
        return
    end if

    ! 3. Allocate the triangle array
    allocate(STLGeom%tris(n_tris), stat=io_status)
    if (io_status /= 0) then
        ErrStat = ErrID_Fatal
        ErrMsg  = trim(RoutineName)//": Failed to allocate memory for mesh."
        close(iu)
        return
    end if

    ! 4. Rewind for the Second Pass
    rewind(iu)

    ! 5. Second Pass: Read the data
    i = 1
    do
        read(iu, '(A)', iostat=io_status) line
        if (io_status /= 0) exit
        
        line = adjustl(line)
        if (index(line, 'facet normal') == 1) then
            ! Read the normal vector
            read(line, *, iostat=io_status) dummy1, dummy2, STLGeom%tris(i)%n(1), STLGeom%tris(i)%n(2), STLGeom%tris(i)%n(3)
            if (io_status /= 0) exit
            
            ! Skip the 'outer loop' line
            read(iu, *, iostat=io_status) dummy1, dummy2
            if (io_status /= 0) exit
            
            ! Read the 3 vertices (Format: 'vertex X Y Z')
            ! v(1:3, 1) means [x,y,z] of vertex 1
            read(iu, *, iostat=io_status) dummy1, STLGeom%tris(i)%v(1,1), STLGeom%tris(i)%v(2,1), STLGeom%tris(i)%v(3,1); if (io_status /= 0) exit
            read(iu, *, iostat=io_status) dummy1, STLGeom%tris(i)%v(1,2), STLGeom%tris(i)%v(2,2), STLGeom%tris(i)%v(3,2); if (io_status /= 0) exit
            read(iu, *, iostat=io_status) dummy1, STLGeom%tris(i)%v(1,3), STLGeom%tris(i)%v(2,3), STLGeom%tris(i)%v(3,3); if (io_status /= 0) exit
            
            ! Skip 'endloop' and 'endfacet' (we just let the outer loop read past them)
            i = i + 1
        end if
    end do

    if (i<=n_tris) then
        ErrStat = ErrID_Fatal
        ErrMsg  = trim(RoutineName)//": failed to read patch "//trim(num2lstr(i))//" from stl file "//trim(filename)//". "
        close(iu)
        return
    end if

    STLGeom%n_tris = n_tris

    close(iu)

end subroutine read_ascii_stl

subroutine Body_Init(STLGeom,PtfmRefPt,PtfmRefztRot,body,m_body,RhoXg,MSL2SWL,Buoyancy,ErrStat,ErrMsg)
    type(STLGeomType),  intent(in   ) :: STLGeom
    real(ReKi),         intent(in   ) :: PtfmRefPt(3)
    real(ReKi),         intent(in   ) :: PtfmRefztRot
    type(BodyType),     intent(inout) :: body
    type(BodyMiscType), intent(inout) :: m_body
    real(SiKi),         intent(in   ) :: RhoXg
    real(ReKi),         intent(in   ) :: MSL2SWL
    real(ReKi),         intent(  out) :: buoyancy(6)
    integer(IntKi),     intent(  out) :: ErrStat     !< Error status of the operation
    character(*),       intent(  out) :: ErrMsg      !< Error message if ErrStat /= ErrID_None

    ! Local variables
    integer(IntKi)          :: iTri, iCorner, iUnique
    integer(IntKi)          :: current_idx
    real(ReKi)              :: v_raw(3)
    real(ReKi)              :: dist2
    logical                 :: found_match
    integer(IntKi)          :: i, j, iQdrt, n_sub
    real(ReKi)              :: vol(3)
    real(ReKi)              :: q_pos(3,nQdrt)
    real(ReKi)              :: nds(3)
    real(R8Ki)              :: R(3,3)
    real(ReKi)              :: d(3)
    real(ReKi)              :: RhoXgLocal
    real(R8Ki)              :: PtfmRefztRotLocal
    real(ReKi)              :: dF(3), rXnds(3), Force(3), Moment(3)
    type(triangle3D)        :: sub_tris(2)
    integer(IntKi)          :: ErrStat2
    character(ErrMsgLen)    :: ErrMsg2
    character(*), parameter :: RoutineName = "Body_Init"

    ! Tolerance for merging vertices (squared to avoid sqrt() calls)
    real(ReKi), parameter   :: tol = 1.0E-4_ReKi
    real(ReKi), parameter   :: tol2 = tol * tol
    
    ! Initialize variables
    ErrStat = ErrID_None
    ErrMsg  = ""

    RhoXgLocal        = real(RhoXg,ReKi)
    PtfmRefztRotLocal = real(PtfmRefztRot,R8Ki)

    body%PtfmRefPt = PtfmRefPt

    ! 1. Pre-allocate maximum possible sizes
    ! Worst case: every triangle is completely disconnected
    allocate(body%nodes(3, STLGeom%n_tris * 3), stat=ErrStat2)
    if (ErrStat2 /= 0) then
        ErrStat = ErrID_Fatal
        ErrMsg  = trim(RoutineName)//": Failed to allocate body%nodes."
        return
    end if
    allocate(body%tris(3, STLGeom%n_tris), stat=ErrStat2)
    if (ErrStat2 /= 0) then
        ErrStat = ErrID_Fatal
        ErrMsg  = trim(RoutineName)//": Failed to allocate body%tris."
        return
    end if
    
    body%n_nodes = 0_IntKi
    body%n_tris  = STLGeom%n_tris

    ! 2. The Deduplication Loop
    do iTri = 1, STLGeom%n_tris

        ! Loop over the 3 corners of the raw triangle
        do iCorner = 1, 3
            v_raw(:) = STLGeom%tris(iTri)%v(:, iCorner)
            
            found_match = .false.
            
            ! Check against already accepted unique nodes
            do iUnique = 1, body%n_nodes
                ! Calculate squared distance
                dist2 = (body%nodes(1, iUnique) - v_raw(1))**2 + &
                        (body%nodes(2, iUnique) - v_raw(2))**2 + &
                        (body%nodes(3, iUnique) - v_raw(3))**2
                
                if (dist2 < tol2) then
                    ! Match found! Use this existing index.
                    current_idx = iUnique
                    found_match = .true.
                    exit
                end if
            end do
            
            ! If it is a completely new vertex, add it to the list
            if (.not. found_match) then
                body%n_nodes = body%n_nodes + 1
                body%nodes(:, body%n_nodes) = v_raw(:)
                current_idx = body%n_nodes
            end if
            
            ! Assign the index to the triangle connectivity array
            body%tris(iCorner, iTri) = current_idx
            
        end do
    end do

    ! 3. Trim the nodes array to the actual number of unique nodes
    ! (Fortran 2003+ allows reallocation on assignment, but doing it explicitly is safer)
    if (body%n_nodes < STLGeom%n_tris * 3) then
        block
            real(ReKi), allocatable :: temp_nodes(:,:)
            allocate(temp_nodes(3, body%n_nodes),stat=ErrStat2)
            if (ErrStat2 /= 0) then
               ErrStat = ErrID_Fatal
               ErrMsg  = trim(RoutineName)//": Failed to allocate body%tris."
               return
            end if
            temp_nodes(:, :) = body%nodes(:, 1:body%n_nodes)
            deallocate(body%nodes)
            call move_alloc(from=temp_nodes, to=body%nodes)
        end block
    end if

    ! 4. Rotate the nodes based on PtfmRefztRot
    R = reshape([cos(PtfmRefztRotLocal),sin(PtfmRefztRotLocal),0.0_R8Ki,-sin(PtfmRefztRotLocal),cos(PtfmRefztRotLocal),0.0_R8Ki,0.0_R8Ki,0.0_R8Ki,1.0_R8Ki],[3,3])
    body%Nodes = matmul(R,body%nodes)

    ! Allocate MiscVars for this body
    allocate(m_body%Nodes(3,body%n_nodes),stat=ErrStat2)
    if (ErrStat2 /= 0) then
       ErrStat = ErrID_Fatal
       ErrMsg  = trim(RoutineName)//": Failed to allocate m%body%Nodes."
       return
    end if
    allocate(m_body%WaveElev(body%n_nodes),stat=ErrStat2)
    if (ErrStat2 /= 0) then
       ErrStat = ErrID_Fatal
       ErrMsg  = trim(RoutineName)//": Failed to allocate m%body%WaveElev."
       return
    end if
    m_body%Nodes    = 0.0_ReKi
    m_body%WaveElev = 0.0_ReKi

    ! Compute and check total volume
    vol = 0.0_ReKi
    do i = 1,body%n_tris
        associate( v1 => body%Nodes(:,body%tris(1,i)), &
                   v2 => body%Nodes(:,body%tris(2,i)), &
                   v3 => body%Nodes(:,body%tris(3,i)) )
            nds = 0.5_ReKi * cross_product(v2-v1, v3-v2)
            q_pos = matmul( reshape([v1,v2,v3],[3,3]) , Qdrt_L )
            vol = vol + matmul( q_pos, Qdrt_W) * nds
        end associate
    end do
    body%volume = sum(vol)/3.0_ReKi
    do i = 1,3
       if (abs(vol(i)-body%volume)>body%volume*1.0E-6_ReKi) then
          ErrStat = ErrID_Fatal
          ErrMsg  = " Inconsistent volumes computed for nonlinear F-K body. Check mesh validity and gaps. "
          return
       end if
    end do
    if (body%volume<=0.0_ReKi) then
        ErrStat = ErrID_Fatal
        ErrMsg  = " Nonlinear F-K body has negative volume. Check normal direction of stl file. "
        return
    end if

    ! Compute buoyancy on undisplaced structure (moment is about global origin)
    d    = PtfmRefPt
    d(3) = d(3) - MSL2SWL
    m_body%Nodes(1,:) = body%Nodes(1,:) + d(1)
    m_body%Nodes(2,:) = body%Nodes(2,:) + d(2)
    m_body%Nodes(3,:) = body%Nodes(3,:) + d(3)
    Force  = 0.0_ReKi
    Moment = 0.0_ReKi
    do i = 1,body%n_tris
        associate( v1 => m_body%Nodes(:,body%tris(1,i)), &
                   v2 => m_body%Nodes(:,body%tris(2,i)), &
                   v3 => m_body%Nodes(:,body%tris(3,i))  )
            nds = 0.5_ReKi * cross_product(v2-v1, v3-v2)
            call Clip_Triangle(reshape([v1,v2,v3],[3,3]), [0.0_ReKi,0.0_ReKi,0.0_ReKi], nds, sub_tris, n_sub)
            do j=1,n_sub
                q_pos = matmul( sub_tris(j)%v , Qdrt_L )
                do iQdrt = 1,nQdrt
                    dF = -RhoXgLocal * q_pos(3,iQdrt) * Qdrt_W(iQdrt) * (-sub_tris(j)%nds)
                    Force  = Force  + dF
                    Moment = Moment + cross_product(q_pos(:,iQdrt), dF)
                end do
            end do
        end associate
    end do
    buoyancy(1:3) = Force
    buoyancy(4:6) = Moment

end subroutine Body_Init

subroutine computeBodyFK(bodyIdx,Time,p,m,Position,Orientation,Force,Moment,ErrStat,ErrMsg)

    integer(IntKi),                  intent(in   ) :: bodyIdx
    real(R8Ki),                      intent(in   ) :: Time
    type(NonlinearFK_ParameterType), intent(in   ) :: p
    type(NonlinearFK_MiscVarType),   intent(inout) :: m
    real(ReKi),                      intent(in   ) :: Position(3)
    real(R8Ki),                      intent(in   ) :: Orientation(3,3)
    real(ReKi),                      intent(  out) :: Force(3)
    real(ReKi),                      intent(  out) :: Moment(3)
    integer(IntKi),                  intent(  out) :: ErrStat     !< Error status of the operation
    character(*),                    intent(  out) :: ErrMsg      !< Error message if ErrStat /= ErrID_None

    integer(IntKi)          :: i, j, iQdrt, n_sub
    real(ReKi)              :: d(3)
    real(R8Ki)              :: R(3,3)
    real(ReKi)              :: q_pos(3,nQdrt)
    real(ReKi)              :: nds(3)
    real(ReKi)              :: dF(3), rXnds(3)
    type(triangle3D)        :: sub_tris(2)
    integer(IntKi)          :: nodeInWater
    real(SiKi)              :: FDynP
    real(ReKi)              :: RhoXgLocal
    integer(IntKi)          :: ErrStat2     !< Error status of the operation
    character(ErrMsgLen)    :: ErrMsg2      !< Error message if ErrStat /= ErrID_None
    character(*), parameter :: RoutineName = "computeBodyFK"
    
    ErrStat = ErrID_None
    ErrMsg  = ""

    Force  = 0.0_ReKi
    Moment = 0.0_ReKi

    if (p%FKMod(bodyIdx) /= FKMod_full) return

    RhoXgLocal = real(p%WaveField%RhoXg,ReKi)
    associate( body=>p%Bodies(bodyIdx), m_body=>m%Bodies(bodyIdx) )

        d    = Position
        d(3) = d(3) - p%WaveField%MSL2SWL
        R    = transpose(Orientation)

        ! Compute displaced node positions and the wave elevation at these nodes
        m_body%Nodes = matmul(R,body%Nodes)
        m_body%Nodes(1,:) = m_body%Nodes(1,:) + d(1)
        m_body%Nodes(2,:) = m_body%Nodes(2,:) + d(2)
        m_body%Nodes(3,:) = m_body%Nodes(3,:) + d(3)

        if (p%WaveField%WaveStMod /= 0_IntKi) then
            do i = 1,body%n_nodes
                m_body%WaveElev(i) = real( WaveField_GetNodeTotalWaveElev( p%WaveField, m%WaveField_m, Time, m_body%Nodes(:,i), ErrStat2, ErrMsg2 ), ReKi)
            end do
            if (Failed()) return
        ! else
            ! m_body%WaveElev already initialized to zero at initialization
        end if

        do i = 1,body%n_tris
            associate( v1    => m_body%Nodes(:,body%tris(1,i)),  &
                       v2    => m_body%Nodes(:,body%tris(2,i)),  &
                       v3    => m_body%Nodes(:,body%tris(3,i)),  &
                       zeta1 => m_body%WaveElev(body%tris(1,i)), &
                       zeta2 => m_body%WaveElev(body%tris(2,i)), &
                       zeta3 => m_body%WaveElev(body%tris(3,i))  )
                nds = 0.5_ReKi * cross_product(v2-v1, v3-v2)
                call Clip_Triangle(reshape([v1,v2,v3],[3,3]), [zeta1,zeta2,zeta3], nds, sub_tris, n_sub)
                do j=1,n_sub
                    q_pos = matmul( sub_tris(j)%v , Qdrt_L )
                    do iQdrt = 1,nQdrt
                        call WaveField_GetDynP( p%WaveField, m%WaveField_m, Time, q_pos(:,iQdrt), .false., nodeInWater, FDynP, ErrStat2, ErrMsg2 )
                        dF = ( real(FDynP,ReKi) - RhoXgLocal * q_pos(3,iQdrt) ) * Qdrt_W(iQdrt) * (-sub_tris(j)%nds)
                        Force  = Force  + dF
                        Moment = Moment + cross_product(q_pos(:,iQdrt)-d, dF)
                    end do
                end do
            end associate
        end do
        if (Failed()) return

    end associate

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
      Failed =  ErrStat >= AbortErrLev
   end function Failed

end subroutine computeBodyFK

subroutine NonlinearFK_CalcOutput(Time, u_Mesh, p, m, Force, Moment, ErrStat, ErrMsg)
    real(R8Ki),                      intent(in   ) :: Time
    type(MeshType),                  intent(in   ) :: u_Mesh
    type(NonlinearFK_ParameterType), intent(in   ) :: p
    type(NonlinearFK_MiscVarType),   intent(inout) :: m
    real(ReKi),                      intent(  out) :: Force(:,:)
    real(ReKi),                      intent(  out) :: Moment(:,:)
    integer(IntKi),                  intent(  out) :: ErrStat     !< Error status of the operation
    character(*),                    intent(  out) :: ErrMsg      !< Error message if ErrStat /= ErrID_None

    integer(IntKi)          :: i
    real(ReKi)              :: Position(3)
    real(R8Ki)              :: Orientation(3,3)
    integer(IntKi)          :: ErrStat2     !< Error status of the operation
    character(ErrMsgLen)    :: ErrMsg2      !< Error message if ErrStat /= ErrID_None
    character(*), parameter :: RoutineName = "NonlinearFK_CalcOutput"

    ErrStat = ErrID_None
    ErrMsg  = ""

    do i = 1,p%nBody
        Position = u_Mesh%position(:,i) + u_Mesh%TranslationDisp(:,i)
        Orientation = u_Mesh%Orientation(:,:,i)
        call computeBodyFK(i,Time,p,m,Position,Orientation,Force(:,i),Moment(:,i),ErrStat2,ErrMsg2)
        if (Failed()) return
    end do

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
      Failed =  ErrStat >= AbortErrLev
   end function Failed

end subroutine NonlinearFK_CalcOutput

subroutine Clip_Triangle(v_in, zeta_in, nds_orig, sub_tris, n_sub)
    real(ReKi),       intent(in   ) :: v_in(3, 3)  ! Original vertices (Global)
    real(ReKi),       intent(in   ) :: zeta_in(3)  ! Wave elevation at vertices
    real(ReKi),       intent(in   ) :: nds_orig(3) ! Original n*dS vector (for flat triangles, normal is constant)
    type(Triangle3D), intent(  out) :: sub_tris(2) ! Max 2 sub-triangles
    integer(IntKi),   intent(  out) :: n_sub       ! Number of wetted sub-triangles (0, 1, or 2)

    real(ReKi)     :: H(3)
    integer(IntKi) :: n_wet, i, j
    integer(IntKi) :: wet_idx(3), dry_idx(3)
    integer(IntKi) :: n_w, n_d
    real(ReKi)     :: t_interp
    real(ReKi)     :: p1(3), p2(3), z1, z2

    ! 1. Calculate Submergence metric (H <= 0 means submerged)
    H(:) = v_in(3, :) - zeta_in(:)

    ! 2. Categorize vertices into Wet and Dry lists
    n_w = 0
    n_d = 0
    do i = 1, 3
        if (H(i) <= 0.0_ReKi) then
            n_w = n_w + 1
            wet_idx(n_w) = i
        else
            n_d = n_d + 1
            dry_idx(n_d) = i
        end if
    end do

    ! 3. The 4 Clipping Cases
    n_sub = 0

    select case(n_w)
    
        case(0)
            ! Completely Dry
            return

        case(3)
            ! Completely Wet: Return the original triangle
            n_sub = 1
            sub_tris(1)%v(:, :)  = v_in(:, :)
            sub_tris(1)%nds(:)   = nds_orig(:)
            return

        case(1)
            ! 1 Vertex Wet -> Results in 1 smaller Triangle
            n_sub = 1
            
            ! The one wet vertex is V_w
            ! The two dry vertices are V_d1, V_d2
            
            ! Intersection 1: Between V_w and V_d1
            t_interp = H(wet_idx(1)) / (H(wet_idx(1)) - H(dry_idx(1)))
            p1(:) = v_in(:, wet_idx(1)) + t_interp * (v_in(:, dry_idx(1)) - v_in(:, wet_idx(1)))

            ! Intersection 2: Between V_w and V_d2
            t_interp = H(wet_idx(1)) / (H(wet_idx(1)) - H(dry_idx(2)))
            p2(:) = v_in(:, wet_idx(1)) + t_interp * (v_in(:, dry_idx(2)) - v_in(:, wet_idx(1)))

            ! Build the new wetted triangle
            sub_tris(1)%v(:, 1) = v_in(:, wet_idx(1))
            sub_tris(1)%v(:, 2) = p1(:)
            sub_tris(1)%v(:, 3) = p2(:)
                        
            ! Calculate new area vector (nds)
            sub_tris(1)%nds = 0.5_ReKi * cross_product(p1 - v_in(:, wet_idx(1)), p2 - v_in(:, wet_idx(1)))
            
            ! Ensure normal points in the same direction as the original
            if (dot_product(sub_tris(1)%nds, nds_orig) < 0.0_ReKi) then
                 sub_tris(1)%nds = -sub_tris(1)%nds
            end if

        case(2)
            ! 2 Vertices Wet -> Results in a Quadrilateral -> Split into 2 Triangles
            n_sub = 2
            
            ! The two wet vertices are V_w1, V_w2
            ! The one dry vertex is V_d
            
            ! Intersection 1: Between V_w1 and V_d
            t_interp = H(wet_idx(1)) / (H(wet_idx(1)) - H(dry_idx(1)))
            p1(:) = v_in(:, wet_idx(1)) + t_interp * (v_in(:, dry_idx(1)) - v_in(:, wet_idx(1)))

            ! Intersection 2: Between V_w2 and V_d
            t_interp = H(wet_idx(2)) / (H(wet_idx(2)) - H(dry_idx(1)))
            p2(:) = v_in(:, wet_idx(2)) + t_interp * (v_in(:, dry_idx(1)) - v_in(:, wet_idx(2)))

            ! Sub-Triangle 1: V_w1, p1, p2
            sub_tris(1)%v(:, 1) = v_in(:, wet_idx(1))
            sub_tris(1)%v(:, 2) = p1(:)
            sub_tris(1)%v(:, 3) = p2(:)
            
            ! Sub-Triangle 2: V_w1, p2, V_w2
            sub_tris(2)%v(:, 1) = v_in(:, wet_idx(1))
            sub_tris(2)%v(:, 2) = p2(:)
            sub_tris(2)%v(:, 3) = v_in(:, wet_idx(2))

            ! Calculate new area vectors
            sub_tris(1)%nds = 0.5_ReKi * cross_product(p1 - v_in(:, wet_idx(1)), p2 - v_in(:, wet_idx(1)))
            sub_tris(2)%nds = 0.5_ReKi * cross_product(p2 - v_in(:, wet_idx(1)), v_in(:, wet_idx(2)) - v_in(:, wet_idx(1)))

            ! Ensure normals point outwards
            if (dot_product(sub_tris(1)%nds, nds_orig) < 0.0_ReKi) sub_tris(1)%nds = -sub_tris(1)%nds
            if (dot_product(sub_tris(2)%nds, nds_orig) < 0.0_ReKi) sub_tris(2)%nds = -sub_tris(2)%nds

    end select

end subroutine Clip_Triangle

end module NonlinearFK