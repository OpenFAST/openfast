module TestMeshMapping_Mod
   
   USE NWTC_Library
   
   implicit none

   TYPE(meshtype) :: mesh1_I, mesh1_O
   TYPE(meshtype) :: mesh2_I, mesh2_O 

   TYPE(meshtype) :: mesh1_I_op, mesh1_O_op
   TYPE(meshtype) :: mesh2_I_op, mesh2_O_op 
   
   
   TYPE(MeshMapType)       :: Map_Mod1_Mod2        ! Data for mapping meshes from mod1 to mod2
   TYPE(MeshMapType)       :: Map_Mod2_Mod1        ! Data for mapping meshes from mod1 to mod2
   
   REAL(R8Ki)              :: Orientation(3,3), angles(3)
   REAL(ReKi)              :: position(3)
   REAL(ReKi)              :: Angle, e, mf, mm
   
   CHARACTER(*), PARAMETER :: Fmt = '(ES10.3E2)'
   CHARACTER(*), PARAMETER :: ErrFmt = '(A28,i2,1x,F18.8,1x,F18.4,13(1x,F18.8))'
   CHARACTER(*), PARAMETER :: ErrTxtFmt = '(A28,A2,15(1x,A18))'
   
   !
   INTEGER :: NNodes, I,J, n1, n2
      
   INTEGER :: Mesh1Type
   INTEGER :: Mesh2Type
   
   INTEGER(IntKi) :: ErrStat
   CHARACTER(1024) :: ErrMsg   
   
   logical :: debug_print = .false.
   
      
        
contains   
      
   ! ..............................................   
   subroutine getRotationPerturb(Vec)
      real(reki), intent(inout) :: vec(:)
      
      real(reki)                :: mx,mn
      integer                   :: nq
      
      !call getRandomVector(Vec, 0.4_ReKi)
      
      call random_number( Vec )
      mx = maxval(Vec)
      mn = minval(Vec)
            
      Vec = 0.4_ReKi*(2.*(Vec-mn)/(mx-mn)-1.)/real(n1**2)  !note n1 from the loop this is called in
          
      do nq=2,size(Vec),3
         Vec(nq-1)=0
         Vec(nq+1)=0
      end do
                 
      
      !vec = 0.0_ReKi
   
   end subroutine getRotationPerturb
   ! ..............................................   
   subroutine getRandomVector(Vec, InLimits)
      real(reki), intent(inout)           :: vec(:)
      real(reki), intent(in   ), optional :: InLimits
      real(reki)                          :: Limits

      real(reki)                          :: mx,mn
      
      if (present(InLimits)) then
         Limits = abs(InLimits)
         if (EqualRealNos( Limits, 0.0_ReKi ) ) Limits = 1.0_ReKi
      else
         Limits = 1.0_ReKi
      end if
      
      call random_number( Vec )
      mx = maxval(Vec)
      mn = minval(Vec)

      Vec = Limits*(2.*(Vec-mn)/(mx-mn)-1.)/real(n1**3)  !note n1 from the loop this is called in
      !Vec = (Vec)/real(n1**3)  !note n1 from the loop this is called in
      
   end subroutine getRandomVector
   ! ..............................................   
   subroutine getLinearOrient(theta, Orientation)
      real(reki), intent(in   )           :: theta(3)
      real(R8ki), intent(inout)           :: Orientation(3,3)
   
      !call eye(Orientation,ErrStat,ErrMsg)
      !Orientation = Orientation - SkewSymMat(theta)
            
      call SmllRotTrans( 'orientation', theta(1) &
                                      , theta(2) &
                                      , theta(3) & 
                                      , Orientation, ErrStat=ErrStat, ErrMsg=ErrMsg)  
      
      
   end subroutine getLinearOrient
   ! ..............................................   
   subroutine WrErrorLine(Desc, Actual, Approx, DestField, DeltaS1, DeltaS2, DeltaS3, DeltaS4)
   character(*), intent(in)            :: Desc
   real(reki)  , intent(in)            :: Actual(:)
   real(reki)  , intent(in)            :: Approx(:)
   real(reki)  , intent(in)            :: DestField(:,:)
   real(reki)  , intent(in)            :: DeltaS1(:)
   real(reki)  , intent(in), optional  :: DeltaS2(:)
   real(reki)  , intent(in), optional  :: DeltaS3(:)
   real(reki)  , intent(in), optional  :: DeltaS4(:)
   
   real(reki)                          :: MaxAbsErr, MaxRelErr
   real(reki)                          :: AbsErr,    RelErr, Denom
   real(reki)                          :: v(3)
   
   integer(intki)                      :: i
   REAL(ReKi), parameter               :: mmin = 1D-7; !sqrt(epsilon(mmin)) !1D-7 !
   
   
   MaxAbsErr = 0.0_ReKi
   MaxRelErr = 0.0_ReKi

   do i=1,size(Actual)-1,3
      v = Actual(i:i+2) - Approx(i:i+2)
      AbsErr = TwoNorm( v )

      v = Actual(i:i+2)
      Denom = TwoNorm(v)
      RelErr = AbsErr / max(mmin, Denom)
      
      if (AbsErr > MaxAbsErr) then
         MaxAbsErr = max(MaxAbsErr,AbsErr) 
         
         MaxRelErr = RelErr !jmj wants the relative error associated with the node with greatest absolute error
         !MaxRelErr = max(MaxAbsErr,RelErr) 
      end if
      
      if (debug_print) then
         call WrNumAryFileNR ( 58, (/ AbsErr /), 'ES15.8', ErrStat, ErrMsg  )   
         call WrNumAryFileNR ( 59, (/ RelErr /), 'ES15.8', ErrStat, ErrMsg  )   
         call WrNumAryFileNR ( 60, (/ Denom /),  'ES15.8', ErrStat, ErrMsg  )   
      end if
      
   end do
   
   if (debug_print) then
      write( 58, *)   
      write( 59, *)   
      write( 60, *)   
   end if
      
   
   !do i=1,size(Actual)
   !   
   !   AbsErr = abs(Actual(i) - Approx(i))
   !   if (AbsErr > MaxAbsErr) then
   !      MaxAbsErr = max(MaxAbsErr,AbsErr) 
   !         
   !      RelErr = AbsErr / max(mmin, Actual(i))
   !      
   !      MaxRelErr = RelErr !jmj wants the relative error associated with the node with greatest absolute error
   !      !MaxRelErr = max(MaxAbsErr,RelErr) 
   !   end if
   !                     
   !end do

   MaxRelErr = MaxRelErr*100.
   !note: n1 comes from loops where this is called (global variable)!
   
   if (present(DeltaS4)) then
      write(*,ErrFmt) Desc//': ', n1, MaxAbsErr, MaxRelErr, maxval(abs(Actual)), maxval(abs(DestField)), maxval(abs(DeltaS1)), maxval(abs(DeltaS2)), maxval(abs(DeltaS3)), maxval(abs(DeltaS4))
   elseif (present(DeltaS3)) then
      write(*,ErrFmt) Desc//': ', n1, MaxAbsErr, MaxRelErr, maxval(abs(Actual)), maxval(abs(DestField)), maxval(abs(DeltaS1)), maxval(abs(DeltaS2)), maxval(abs(DeltaS3))
   elseif (present(DeltaS2)) then
      write(*,ErrFmt) Desc//': ', n1, MaxAbsErr, MaxRelErr, maxval(abs(Actual)), maxval(abs(DestField)), maxval(abs(DeltaS1)), maxval(abs(DeltaS2))
   else
      write(*,ErrFmt) Desc//': ', n1, MaxAbsErr, MaxRelErr, maxval(abs(Actual)), maxval(abs(DestField)), maxval(abs(DeltaS1))
   end if
   
      
   end subroutine
   ! .............................................. 
   subroutine CreateOutputMeshes_Test1()   
      ! this is a point-to-point mapping, with one point going to many.
      ! it is a figure in the AIAA paper.
   
      real(ReKi)   :: dz
         
      Mesh1Type = ELEMENT_POINT
      Mesh2Type = ELEMENT_POINT
      
      !.........................
      ! Mesh1 (Output: Motions)
      !.........................
      
      Nnodes = 1
            
      CALL MeshCreate( BlankMesh          = mesh1_O           &
                       , IOS              = COMPONENT_OUTPUT  &
                       , NNodes           = NNodes            &
                       , Orientation      = .TRUE.            &
                       , TranslationDisp  = .TRUE.            &
                       , TranslationVel   = .TRUE.            &
                       , RotationVel      = .TRUE.            &
                       , TranslationAcc   = .TRUE.            &
                       , RotationAcc      = .TRUE.            &                                
                       , ErrStat          = ErrStat           &
                       , ErrMess          = ErrMsg            )
      
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

                                    
      do j=1,NNodes 
            ! place nodes in a line
         position = (/0.0_ReKi, 0.0_ReKi, 1.0_ReKi*(j-1) /)
         CALL MeshPositionNode ( mesh1_O, j, position, ErrStat, ErrMsg )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))             

         CALL MeshConstructElement ( mesh1_O, ELEMENT_POINT, ErrStat, ErrMsg, j )
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
         
      END DO   
         
      CALL MeshCommit ( mesh1_O, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh1 output for test 1.")

      !..............
      ! initialize output fields:
      !..............      
      
      do j=1,Mesh1_O%NNodes
      
   !      Angle = 0      
   !      Angle = (20*j)*D2R      
         Angle = (20)*D2R      
         !note this "looks" like the transpose, but isn't
         Mesh1_O%Orientation(:,:,j) = GetDCM(Angle, 1)
            
         Mesh1_O%TranslationDisp(:,j) = (/ 2.00, 0.,  0. /)
         Mesh1_O%TranslationVel(:,j)  = (/ 0.5,  0.0, 0.0 /)
         Mesh1_O%RotationVel(:,j)     = (/ 2.0,  0.0, 0.0 /)
         Mesh1_O%TranslationAcc(:,j)  = 0.0_ReKi ! (/ 1., 1., 0. /)*.115
         Mesh1_O%RotationAcc(:,j)     = (/ 2.0, 0.0, 0.0 /)
      
      end do
      
               
      !.........................
      ! Mesh2 (Output: Loads)
      !.........................
            
        
      NNodes = 5
      dz = 1.0/(Nnodes-1)
   
      CALL MeshCreate(  BlankMesh       = mesh2_O           &
                        , IOS           = COMPONENT_OUTPUT  &
                        , NNodes        = NNodes            &
                        , Force         = .TRUE.            &
                        , Moment        = .TRUE.            &
                        , ErrStat       = ErrStat           &
                        , ErrMess       = ErrMsg            )   
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

      do j=1,NNodes
                  
         Angle = (-25. + j*j)*D2R  !note this "looks" like the transpose, but isn't
         Orientation(:,1) = (/ COS(Angle), -1.*SIN(Angle), 0.0_ReKi /)
         Orientation(:,2) = (/ SIN(Angle),     COS(Angle), 0.0_ReKi /)
         Orientation(:,3) = (/      0.,        0.0,        1.0 /)
      
            ! place nodes in a line
         position = (/0.0_ReKi, 0.0_ReKi, dz*(j-1) /)
         CALL MeshPositionNode ( mesh2_O, j, position, ErrStat, ErrMsg, &
               Orient= Orientation )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   
         CALL MeshConstructElement ( mesh2_O, ELEMENT_POINT, ErrStat, ErrMsg, j )
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))            
         

      END DO
                  
         ! that's our entire mesh:
      CALL MeshCommit ( mesh2_O, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))   
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh2 output for test 1.")
         
      !..............
      ! initialize output fields:
      !..............
      do j=1,Mesh2_O%NNodes
         Mesh2_O%Force( :,j) = (/  1.0, 0.,  0.   /)  !*(j*0.5)
         Mesh2_O%Moment(:,j) = (/  0.0, 0.5, 0.5  /)*(-j*0.0)
      end do

                  
   end subroutine CreateOutputMeshes_Test1
   ! ..............................................
   subroutine CreateOutputMeshes_Test2(ThisCase)   
      ! this is a line2-to-line2 mapping, with the same meshes in an upside-down "T" shape, though
      ! one is discretized more than the other.
      character(1), intent(in) :: ThisCase
      real(reki) :: dx, dz
      integer NumHorizNodes, NumVertNodes, ConnectionNode
      REAL(ReKi), PARAMETER :: Len_Horiz = 4
      REAL(ReKi), PARAMETER :: Len_Vert  = 2
      type(meshtype) :: rotateMesh
      type(meshmaptype) :: rotateMesh_map

   
      Mesh1Type = ELEMENT_LINE2
      Mesh2Type = ELEMENT_LINE2
      
      !.........................
      ! Mesh1 (Output: Motions)
      !.........................
      SELECT CASE (ThisCase)
      CASE ('A')
         NumHorizNodes  = 9  
         ConnectionNode = 5
         NumVertNodes   = 8
      CASE ('B','C','D')
         NumHorizNodes  =  5
         ConnectionNode =  3
         NumVertNodes   =  2
      END SELECT
      
            
      Nnodes = NumHorizNodes + NumVertNodes !7 
      dx = Len_Horiz/(NumHorizNodes-1) !1.0
      dz = Len_Vert/(NumVertNodes) !1.0
   
      CALL MeshCreate( BlankMesh          = mesh1_O           &
                       , IOS              = COMPONENT_OUTPUT  &
                       , NNodes           = NNodes            &
                       , Orientation      = .TRUE.            &
                       , TranslationDisp  = .TRUE.            &
                       , TranslationVel   = .TRUE.            &
                       , RotationVel      = .TRUE.            &
                       , TranslationAcc   = .TRUE.            &
                       , RotationAcc      = .TRUE.            &                                
                       , ErrStat          = ErrStat           &
                       , ErrMess          = ErrMsg            )
      
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      

      CALL PositionNodesElements(mesh1_O, NumHorizNodes, ConnectionNode, dx,dz,ErrStat,ErrMsg)                 
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh1 output for test 2.")

      !..............
      ! initialize Mesh1 output fields:
      !..............      
      
      !... temp mesh, for rotation only:
      CALL MeshCreate(BlankMesh  = RotateMesh                 &
                       , IOS              = COMPONENT_OUTPUT  &
                       , NNodes           = 1                 &
                       , Orientation      = .TRUE.            &
                       , TranslationDisp  = .TRUE.            &
                       , TranslationVel   = .TRUE.            &
                       , TranslationAcc   = .TRUE.            &
                       , RotationVel      = .TRUE.            &
                       , RotationAcc      = .TRUE.            &
                       , ErrStat          = ErrStat           &
                       , ErrMess          = ErrMsg            )
      
      CALL MeshPositionNode(RotateMesh, 1, mesh1_O%Position(:,ConnectionNode), ErrStat=ErrStat, ErrMess=ErrMsg, &
                            Orient=mesh1_O%RefOrientation(:,:,ConnectionNode) )      ; IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      CALL MeshConstructElement ( RotateMesh, ELEMENT_POINT, ErrStat, ErrMsg, 1 )    ; IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      CALL MeshCommit ( RotateMesh, ErrStat, ErrMsg )                                ; IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   
      Angle = (15.)*D2R      
      !note this "looks" like the transpose, but isn't
      RotateMesh%Orientation(:,1,1) = (/  COS(Angle) , -1.*SIN(Angle), 0.0_ReKi /)
      RotateMesh%Orientation(:,2,1) = (/  SIN(Angle),      COS(Angle), 0.0_ReKi /)
      RotateMesh%Orientation(:,3,1) = (/  0.0,         0.0,            1.0 /)
      RotateMesh%TranslationDisp(:,1) = (/ 1.00, 1.50,  0.00 /) 
      RotateMesh%TranslationVel( :,1) = (/   Pi,   Pi,  0.00_ReKi /) *0.0
      RotateMesh%RotationVel(    :,1) = (/ 0.00, 0.00,  0.5 /) 
       
      CALL MeshMapCreate( RotateMesh, Mesh1_O, rotateMesh_map, ErrStat, ErrMsg) ; IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      CALL Transfer_Point_to_Line2( RotateMesh, Mesh1_O, rotateMesh_map, ErrStat, ErrMsg); IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      
      call MeshDestroy(rotateMesh,ErrStat,ErrMsg)
      call MeshMapDestroy( rotateMesh_map,ErrStat,ErrMsg)
                                       
      !.........................
      ! Mesh2 (Output: Loads)
      !.........................
                 
      SELECT CASE (ThisCase)
      CASE ('B')
         NumHorizNodes  = 9  
         ConnectionNode = 5
         NumVertNodes   = 8
      CASE ('A','C')
         NumHorizNodes  =  5
         ConnectionNode =  3
         NumVertNodes   =  2
      CASE ('D')
         NumHorizNodes  =  3
         ConnectionNode =  2
         NumVertNodes   =  2
      END SELECT
      
      
      dx = Len_Horiz/(NumHorizNodes-1) !1/2
      dz = Len_Vert/(NumVertNodes) !1/4
      
      Nnodes = NumHorizNodes + NumVertNodes !17 
               
      CALL MeshCreate(  BlankMesh       = mesh2_O           &
                        , IOS           = COMPONENT_OUTPUT  &
                        , NNodes        = NNodes            &
                        , Force         = .TRUE.            &
                        , Moment        = .TRUE.            &
                        , ErrStat       = ErrStat           &
                        , ErrMess       = ErrMsg            )   
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
            
      CALL PositionNodesElements(mesh2_O, NumHorizNodes, ConnectionNode, dx,dz,ErrStat,ErrMsg)                 
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh2 output for test 2.")
                                           
            
      !..............
      ! initialize Mesh2 output fields:
      !..............
      
      Mesh2_O%Force    = 0.0
      Mesh2_O%Moment   = 0.0
      
      Mesh2_O%Force(3, 1: NumHorizNodes) = 1.         
      Mesh2_O%Force(1, ConnectionNode)   = 0.5
      do j=NumHorizNodes+1,Mesh2_O%NNodes
         Mesh2_O%Force(1, j) = (0.5*Mesh2_O%Position(3,j))**5 + Mesh2_O%Force(1, ConnectionNode)
      end do
                       
      
      
      return
   end subroutine CreateOutputMeshes_Test2   
   !-------------------------------------------------------------------------------------------------
   subroutine PositionNodesElements(ThisMesh, NumHorizNodes, ConnectionNode, dx,dz,ErrStat,ErrMsg)
   ! this positions nodes, creates elements, and committs the mesh for Test2
   
      TYPE(MeshType), INTENT(INOUT) :: ThisMesh
      real(reki),     INTENT(IN)    :: dx, dz
      integer,        INTENT(IN)    :: NumHorizNodes, ConnectionNode
      INTEGER(IntKi), INTENT(OUT)   :: ErrStat      
      character(*),   INTENT(OUT)   :: ErrMsg      
      
      integer :: j
      
            ! create a "T"
            
      ! position nodes:
            
            
      ! horizontal line:            
      DO j=1,NumHorizNodes 
         position = (/dx*(j-ConnectionNode), 0.0_ReKi, 0.0_ReKi  /)
         CALL MeshPositionNode ( ThisMesh, j, position, ErrStat, ErrMsg )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))             
      END DO   
      
      
      ! vertical line:
      DO j=(NumHorizNodes+1),ThisMesh%Nnodes
         position = (/0.0_ReKi, 0.0_ReKi, dz*(j-NumHorizNodes)  /)
         CALL MeshPositionNode ( ThisMesh, j, position, ErrStat, ErrMsg )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))             
      END DO                  
   
         ! construct elements:

        
      ! horizontal line:
      do j=2,NumHorizNodes
         CALL MeshConstructElement ( ThisMesh, ELEMENT_LINE2, ErrStat, ErrMsg, P1=j-1, P2=J )
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))            
      END DO
         
      ! vertical line:
      CALL MeshConstructElement ( ThisMesh, ELEMENT_LINE2, ErrStat, ErrMsg, P1=ConnectionNode, P2=NumHorizNodes+1 )
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))            
         
      DO j=(NumHorizNodes+2),ThisMesh%Nnodes      
         CALL MeshConstructElement ( ThisMesh, ELEMENT_LINE2, ErrStat, ErrMsg, P1=j-1, P2=j )
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))            
      END DO
                                    
         ! that's our entire mesh:
      CALL MeshCommit ( ThisMesh, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
      if (ErrStat >= AbortErrLev) CALL ProgAbort("")        
               
   end subroutine PositionNodesElements   
   !-------------------------------------------------------------------------------------------------
   subroutine CreateOutputMeshes_Test5
   
      real(reki) :: z
      
      
      Mesh1Type = ELEMENT_LINE2
      Mesh2Type = ELEMENT_LINE2
      
      !.........................
      ! Mesh1 (Output: Motions)
      !.........................
      
      Nnodes = 6            
      CALL MeshCreate( BlankMesh          = mesh1_O           &
                       , IOS              = COMPONENT_OUTPUT  &
                       , NNodes           = NNodes            &
                       , Orientation      = .TRUE.            &
                       , TranslationDisp  = .TRUE.            &
                       , TranslationVel   = .TRUE.            &
                       , RotationVel      = .TRUE.            &
                       , TranslationAcc   = .TRUE.            &
                       , RotationAcc      = .TRUE.            &                                
                       , ErrStat          = ErrStat           &
                       , ErrMess          = ErrMsg            )
      
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

                                    
      do j=1,NNodes 
            ! place nodes in a line
         position = (/0.0_ReKi, 0.0_ReKi, 1.0_ReKi*(j-1)/(Nnodes-1) /)
         CALL MeshPositionNode ( mesh1_O, j, position, ErrStat, ErrMsg )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                      
      END DO   


      do j=2,NNodes 
         CALL MeshConstructElement ( mesh1_O, ELEMENT_LINE2, ErrStat, ErrMsg, J-1, j )
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
      end do
      
      
      CALL MeshCommit ( mesh1_O, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh1 output for test 3.")

      !..............
      ! initialize output fields:
      !..............      
      
      do j=1,Mesh1_O%NNodes
      
         !Angle = 0      
         Angle = 0.5*Mesh1_O%Position(3,j)
         
         !note this "looks" like the transpose, but isn't
         Mesh1_O%Orientation(:,1,j) = (/  COS(Angle),   0.0_ReKi,  1.*SIN(Angle) /)
         Mesh1_O%Orientation(:,2,j) = (/  0.0,               1.0,            0.0 /)
         Mesh1_O%Orientation(:,3,j) = (/-1.*SIN(Angle), 0.0_ReKi,     COS(Angle) /)
                     
         Mesh1_O%TranslationDisp(:,j) = (/ 2.*(1.0-COS(Angle)), 0._ReKi, 2.0*SIN(Angle)-Mesh1_O%Position(3,j) /)

         Mesh1_O%TranslationVel( :,j) = Mesh1_O%TranslationDisp(:,j)*1.5
         Mesh1_O%RotationVel(:,j)     = (/ 0.0_ReKi, Angle, 0.0_ReKi /)*1.5
         Mesh1_O%TranslationAcc(:,j)  = 0.0_ReKi 
         Mesh1_O%RotationAcc(:,j)     = 0.0_ReKi
                           
      end do
      
               
      !.........................
      ! Mesh2 (Output: Loads)
      !.........................
                 
      NNodes = 8
   
      CALL MeshCreate(  BlankMesh       = mesh2_O           &
                        , IOS           = COMPONENT_OUTPUT  &
                        , NNodes        = NNodes            &
                        , Force         = .TRUE.            &
                        , Moment        = .TRUE.            &
                        , ErrStat       = ErrStat           &
                        , ErrMess       = ErrMsg            )   
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

      do j=1,NNodes
                  
         Angle = 30.*D2R  !note this "looks" like the transpose, but isn't
         !Angle = 0
         Orientation(:,1) = (/    COS(Angle), 0.0_ReKi,  1.*SIN(Angle) /)
         Orientation(:,2) = (/    0.0,             1.0,            0.0 /)
         Orientation(:,3) = (/-1.*SIN(Angle), 0.0_ReKi,     COS(Angle) /)
      
            ! place nodes in a line
         z = (j-1.0)/(NNodes-1.0)
         position = (/(sqrt(3.0)/3.0)*z-(sqrt(3.0)/6.0), 0.0_ReKi, z /)
         CALL MeshPositionNode ( mesh2_O, j, position, ErrStat, ErrMsg, &
               Orient= Orientation )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))        

      END DO
      
      do j=2,NNodes 
         CALL MeshConstructElement ( mesh2_O, ELEMENT_LINE2, ErrStat, ErrMsg, J-1, j )
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
      end do      
                  
         ! that's our entire mesh:
      CALL MeshCommit ( mesh2_O, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))   
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh2 output for test 1.")
         
      !..............
      ! initialize output fields:
      !..............
      ! This is done in another subroutine
      
      !do j=1,Mesh2_O%NNodes
      !   Mesh2_I%Orientation(:,:,j)
      !   
      !   Mesh2_O%Force( :,j) = Mesh2_I%Orientation(:,1,j)*5
      !   Mesh2_O%Moment(:,j) = 0.0 ! (/  0.0, 0.0, 0.5  /)*(-j*0.0)
      !end do
      

                  
   end subroutine CreateOutputMeshes_Test5
   subroutine InitTest5Loads()
      !..............
      ! initialize output fields:
      !..............
      do j=1,Mesh2_O%NNodes         
         Mesh2_O%Force( :,j) = Mesh2_I%Orientation(1,:,j)*0.75
         Mesh2_O%Moment(:,j) = 0.0 ! (/  0.0, 0.0, 0.5  /)*(-j*0.0)
      end do   
      
   end subroutine InitTest5Loads   
   subroutine CreateOutputMeshes_Test5Orient
   
      real(reki) :: z
      real(r8ki) :: DCM(3,3)
      
      
      Mesh1Type = ELEMENT_LINE2
      Mesh2Type = ELEMENT_LINE2
         
      !.........................
      ! Mesh1 (Output: Motions)
      !.........................
      
      Nnodes =5 !100 !100 !5            
      CALL MeshCreate( BlankMesh          = mesh1_O           &
                       , IOS              = COMPONENT_OUTPUT  &
                       , NNodes           = NNodes            &
                       , Orientation      = .TRUE.            &
                       , TranslationDisp  = .TRUE.            &
                       , TranslationVel   = .TRUE.            &
                       , RotationVel      = .TRUE.            &
                       , TranslationAcc   = .TRUE.            &
                       , RotationAcc      = .TRUE.            &                                
                       , ErrStat          = ErrStat           &
                       , ErrMess          = ErrMsg            )
      
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

                                    
      do j=1,NNodes 
            ! place nodes in a line
            
         position = (/0.0_ReKi, 0.0_ReKi, 1.0_ReKi*(j-1)/(Nnodes-1) /)          
         
         !CALL MeshPositionNode ( mesh1_O, j, position, ErrStat, ErrMsg, orient=orientation )     
         CALL MeshPositionNode ( mesh1_O, j, position, ErrStat, ErrMsg )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                      
      END DO   

      !do j=1,NNodes 
      !   CALL MeshConstructElement ( mesh1_O, Mesh1Type, ErrStat, ErrMsg, j )
      !   IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
      !end do      

      do j=2,NNodes 
         CALL MeshConstructElement ( mesh1_O, Mesh1Type, ErrStat, ErrMsg, J-1, j )
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
      end do
      
      
      CALL MeshCommit ( mesh1_O, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh1 output for test 3.")

      !..............
      ! initialize output fields:
      !..............      
      
      do j=1,Mesh1_O%NNodes
      
         !Angle = 0      
         !Angle = 0.5* Mesh1_O%Position(3,j) !28.6479
         Angle = 45.0_ReKi*Mesh1_O%Position(3,j)*D2R
         Mesh1_O%Orientation(:,:,j) = GetDCM(Angle,2)
                              
         Mesh1_O%TranslationDisp(:,j) = (/ 2.*(1.0-COS(Angle)), 0._ReKi, 2.0*SIN(Angle)-Mesh1_O%Position(3,j) /)

         Mesh1_O%TranslationVel( :,j) = Mesh1_O%TranslationDisp(:,j)*1.5
         Mesh1_O%RotationVel(:,j)     = (/ 0.0_ReKi, Angle, 0.0_ReKi /)*1.5
         Mesh1_O%TranslationAcc(:,j)  = 0.0_ReKi 
         Mesh1_O%RotationAcc(:,j)     = 0.0_ReKi
                           
      end do
      
               
      !.........................
      ! Mesh2 (Output: Loads)
      !.........................
                 
      NNodes = 8 !105 !
   
      CALL MeshCreate(  BlankMesh       = mesh2_O           &
                        , IOS           = COMPONENT_OUTPUT  &
                        , NNodes        = NNodes            &
                        , Force         = .TRUE.            &
                        , Moment        = .TRUE.            &
                        , ErrStat       = ErrStat           &
                        , ErrMess       = ErrMsg            )   
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

      do j=1,NNodes
                  
         Angle = 30.*D2R  !note this "looks" like the transpose, but isn't      
         DCM = GetDCM(Angle,2)
         
            ! place nodes in a line
         z = (j-1.0)/(NNodes-1.0)
         position = (/(sqrt(3.0)/3.0)*z-(sqrt(3.0)/6.0), 0.0_ReKi, z /)
         CALL MeshPositionNode ( mesh2_O, j, position, ErrStat, ErrMsg, Orient=DCM )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))        

      END DO
      
      !do j=1,NNodes 
      !   CALL MeshConstructElement ( mesh2_O, Mesh2Type, ErrStat, ErrMsg, J )
      !   IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
      !end do      
                  
      do j=2,NNodes 
         CALL MeshConstructElement ( mesh2_O, Mesh2Type, ErrStat, ErrMsg, J-1, j )
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
      end do      
      
      
         ! that's our entire mesh:
      CALL MeshCommit ( mesh2_O, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))   
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh2 output for test 1.")
         
      !..............
      ! initialize output fields:
      !..............
      ! This is done in another subroutine
      
      !do j=1,Mesh2_O%NNodes
      !   Mesh2_I%Orientation(:,:,j)
      !   
      !   Mesh2_O%Force( :,j) = Mesh2_I%Orientation(:,1,j)*5
      !   Mesh2_O%Moment(:,j) = 0.0 ! (/  0.0, 0.0, 0.5  /)*(-j*0.0)
      !end do
      

                  
   end subroutine CreateOutputMeshes_Test5Orient   
   ! ..............................................   
   subroutine CreateOutputMeshes_Test6()   
   
      REAL(reKi)  :: dx
      
      
      Mesh1Type = ELEMENT_POINT
      Mesh2Type = ELEMENT_LINE2
      
      !.........................
      ! Mesh1 (Output: Motions)
      !.........................
      
      Nnodes = 1
            
      CALL MeshCreate( BlankMesh          = mesh1_O           &
                       , IOS              = COMPONENT_OUTPUT  &
                       , NNodes           = NNodes            &
                       , Orientation      = .TRUE.            &
                       , TranslationDisp  = .TRUE.            &
                       , TranslationVel   = .TRUE.            &
                       , RotationVel      = .TRUE.            &
                       , TranslationAcc   = .TRUE.            &
                       , RotationAcc      = .TRUE.            &                                
                       , ErrStat          = ErrStat           &
                       , ErrMess          = ErrMsg            )
      
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

                                    
      do j=1,NNodes 
            ! place nodes in a line
         position = (/0.0_ReKi, 0.0_ReKi, 1.0_ReKi*(j-1) /)
         CALL MeshPositionNode ( mesh1_O, j, position, ErrStat, ErrMsg )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))             

         CALL MeshConstructElement ( mesh1_O, ELEMENT_POINT, ErrStat, ErrMsg, j )
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
         
      END DO   
         
      CALL MeshCommit ( mesh1_O, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh1 output for test 4.")

      !..............
      ! initialize output fields:
      !..............      
      
      do j=1,Mesh1_O%NNodes
      
         Angle = 0      
   !      Angle = (20*j)*D2R      
        Angle = (20)*D2R      
         !note this "looks" like the transpose, but isn't
         Mesh1_O%Orientation(:,1,j) = (/  1.0,            0.0 ,            0.0 /)
         Mesh1_O%Orientation(:,2,j) = (/  0.0_ReKi, COS(Angle), -1.*SIN(Angle) /)
         Mesh1_O%Orientation(:,3,j) = (/  0.0_ReKi, SIN(Angle),     COS(Angle) /)
                                 
         Mesh1_O%TranslationDisp(:,j) = (/ 2.00, 0.,  0. /)
         Mesh1_O%TranslationVel(:,j)  = (/ 0.5,  0.0, 0.0 /)
         Mesh1_O%RotationVel(:,j)     = (/ 2.0,  0.0, 0.0 /)
         Mesh1_O%TranslationAcc(:,j)  = 0.0_ReKi ! (/ 1., 1., 0. /)*.115
         Mesh1_O%RotationAcc(:,j)     = (/ 2.0, 0.0, 0.0 /)
                        
      end do
                           
      !.........................
      ! Mesh2 (Output: Loads)
      !.........................
            
        
      NNodes = 5
      dx = 0.25
      CALL MeshCreate(  BlankMesh       = mesh2_O           &
                        , IOS           = COMPONENT_OUTPUT  &
                        , NNodes        = NNodes            &
                        , Force         = .TRUE.            &
                        , Moment        = .TRUE.            &
                        , ErrStat       = ErrStat           &
                        , ErrMess       = ErrMsg            )   
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

      do j=1,NNodes
                  
         Angle = (-25. + j*j)*D2R  !note this "looks" like the transpose, but isn't
         Orientation(:,1) = (/ COS(Angle), -1.*SIN(Angle), 0.0_ReKi /)
         Orientation(:,2) = (/ SIN(Angle),     COS(Angle), 0.0_ReKi /)
         Orientation(:,3) = (/      0.,        0.0,        1.0      /)
      
            ! place nodes in a line
         position = (/0.0_ReKi, 0.0_ReKi, dx*(j-1) /)
         CALL MeshPositionNode ( mesh2_O, j, position, ErrStat, ErrMsg, &
               Orient= Orientation )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))         

      END DO
      
      do j=2,NNodes
         CALL MeshConstructElement ( mesh2_O, ELEMENT_LINE2, ErrStat, ErrMsg, j-1,j )
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))            
      END DO
      
                  
         ! that's our entire mesh:
      CALL MeshCommit ( mesh2_O, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))   
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh2 output for test 4.")
         
      !..............
      ! initialize output fields:
      !..............
      do j=1,Mesh2_O%NNodes
         Mesh2_O%Force( :,j) = (/ Nnodes / (dx*(Nnodes-1)), 0.0_ReKi, 0.0_ReKi /)
         Mesh2_O%Moment(:,j) = 0.0
      end do

                  
   end subroutine CreateOutputMeshes_Test6   
   ! ..............................................   
   subroutine CreateOutputMeshes_Test7(ThisCase)   
   
      character(1), intent(in) :: ThisCase
      integer, parameter :: Nnodes = 10
      
      SELECT CASE (ThisCase)
      CASE ('A')
         Mesh1Type = ELEMENT_POINT
         Mesh2Type = ELEMENT_LINE2
      CASE ('B')
         Mesh1Type = ELEMENT_LINE2
         Mesh2Type = ELEMENT_POINT
      END SELECT                     

      CALL CreateTest7_Motions(Mesh1Type,Nnodes)
      CALL CreateTest7_Loads(Mesh2Type)
      
                  
   end subroutine CreateOutputMeshes_Test7   
   ! ..............................................
   subroutine CreateTest7_Motions(MeshType,Nnodes)
   
      INTEGER, INTENT(IN) :: MeshType
      INTEGER, INTENT(IN) :: Nnodes

      !.........................
      ! Mesh1 (Output: Motions)
      !.........................
                  
      CALL MeshCreate( BlankMesh          = mesh1_O           &
                       , IOS              = COMPONENT_OUTPUT  &
                       , NNodes           = NNodes            &
                       , Orientation      = .TRUE.            &
                       , TranslationDisp  = .TRUE.            &
                       , TranslationVel   = .TRUE.            &
                       , RotationVel      = .TRUE.            &
                       , TranslationAcc   = .TRUE.            &
                       , RotationAcc      = .TRUE.            &                                
                       , ErrStat          = ErrStat           &
                       , ErrMess          = ErrMsg            )
      
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

                                    
      do j=1,NNodes 
            ! place nodes in a line
         position = (/0.0_ReKi, 0.0_ReKi, 1.0_ReKi*(j-1) /)
         CALL MeshPositionNode ( mesh1_O, j, position, ErrStat, ErrMsg )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))             
         
      END DO   

      IF (MeshType == ELEMENT_POINT) THEN
         do j=1,NNodes 
            CALL MeshConstructElement ( mesh1_O, ELEMENT_POINT, ErrStat, ErrMsg, j )
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
         END DO         
      ELSEIF (MeshType == ELEMENT_LINE2) THEN
         do j=2,NNodes 
            CALL MeshConstructElement ( mesh1_O, ELEMENT_LINE2, ErrStat, ErrMsg, j-1,j )
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
         END DO         
      END IF
      
      
      
      CALL MeshCommit ( mesh1_O, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh1 output for test 5.")

      !..............
      ! initialize output fields:
      !..............      
      
      do j=1,Mesh1_O%NNodes
      
         Angle = 0      
   !      Angle = (20*j)*D2R      
   !      Angle = (20)*D2R      
         !note this "looks" like the transpose, but isn't
         Mesh1_O%Orientation(:,1,j) = (/  1.0,       0.0 ,            0.0 /)
         Mesh1_O%Orientation(:,2,j) = (/  0.0_ReKi, COS(Angle), -1.*SIN(Angle) /)
         Mesh1_O%Orientation(:,3,j) = (/  0.0_ReKi, SIN(Angle),     COS(Angle) /)
            
         Mesh1_O%TranslationDisp(:,j) = (/ 2., 0.,  0. /)        
         Mesh1_O%TranslationVel(:,j)  = (/ 2., 0.0, 0.0 /)
         Mesh1_O%RotationVel(:,j)     = 0.0_ReKi 
         Mesh1_O%TranslationAcc(:,j)  = (/  -1., 0.0, 0.0 /)
         Mesh1_O%RotationAcc(:,j)     = 0.0_ReKi 
                         
      end do
                    
   
   
   end subroutine CreateTest7_Motions
   ! ..............................................
   subroutine CreateTest7_Loads(MeshType)
   
      INTEGER, INTENT(IN) :: MeshType
      integer :: nnodes
      
      !.........................
      ! Mesh2 (Output: Loads)
      !.........................
                    
      NNodes = mesh1_O%Nnodes
      CALL MeshCreate(  BlankMesh       = mesh2_O           &
                        , IOS           = COMPONENT_OUTPUT  &
                        , NNodes        = NNodes            &
                        , Force         = .TRUE.            &
                        , Moment        = .TRUE.            &
                        , ErrStat       = ErrStat           &
                        , ErrMess       = ErrMsg            )   
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

      do j=1,NNodes
                        
            ! place nodes in a line
         CALL MeshPositionNode ( mesh2_O, j, mesh1_O%Position(:,j), ErrStat, ErrMsg, &
               Orient= mesh1_O%RefOrientation(:,:,j) )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))         

      END DO
      
      
      
      IF (MeshType == ELEMENT_POINT) THEN
         do j=1,NNodes 
            CALL MeshConstructElement ( mesh2_O, ELEMENT_POINT, ErrStat, ErrMsg, j )
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
         END DO         
      ELSEIF (MeshType == ELEMENT_LINE2) THEN
         do j=2,NNodes
            CALL MeshConstructElement ( mesh2_O, ELEMENT_LINE2, ErrStat, ErrMsg, j-1,j )
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))            
         END DO
      END IF      
            
                  
         ! that's our entire mesh:
      CALL MeshCommit ( mesh2_O, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))   
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh2 output for test 5.")
         
      !..............
      ! initialize output fields:
      !..............
      do j=1,Mesh2_O%NNodes
         Mesh2_O%Force( :,j) = (/ 1.0, 0.0, 0.0 /)
         Mesh2_O%Moment(:,j) = 0.0
      end do   
   
   END subroutine CreateTest7_Loads
   ! ..............................................
   subroutine CreateOutputMeshes_Test9

      Mesh1Type = ELEMENT_Line2
      Mesh2Type = ELEMENT_POINT
      
      !.........................
      ! Mesh1 (Output: Motions)
      !.........................
      
      Nnodes = 2
            
      CALL MeshCreate( BlankMesh          = mesh1_O           &
                       , IOS              = COMPONENT_OUTPUT  &
                       , NNodes           = NNodes            &
                       , Orientation      = .TRUE.            &
                       , TranslationDisp  = .TRUE.            &
                       , TranslationVel   = .TRUE.            &
                       , RotationVel      = .TRUE.            &
                       , TranslationAcc   = .TRUE.            &
                       , RotationAcc      = .TRUE.            &                                
                       , ErrStat          = ErrStat           &
                       , ErrMess          = ErrMsg            )
      
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

                                    
      do j=1,NNodes 
            ! place nodes in a line
         position = (/0.0_ReKi, 0.0_ReKi, 1.0_ReKi*(j-1) /)
         CALL MeshPositionNode ( mesh1_O, j, position, ErrStat, ErrMsg )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                      
      END DO   

      do j=2,NNodes 
         CALL MeshConstructElement ( mesh1_O, ELEMENT_LINE2, ErrStat, ErrMsg, j-1, j )
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
      end do ! 
            
      CALL MeshCommit ( mesh1_O, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh1 output for test 9.")

      !..............
      ! initialize output fields:
      !..............      
      
      do j=1,Mesh1_O%NNodes
      
         Angle = 0      
   !      Angle = (20*j)*D2R      
    !     Angle = (20)*D2R      
         !note this "looks" like the transpose, but isn't
         Mesh1_O%Orientation(:,1,j) = (/  1.0,       0.0 ,            0.0 /)
         Mesh1_O%Orientation(:,2,j) = (/  0.0_ReKi, COS(Angle), -1.*SIN(Angle) /)
         Mesh1_O%Orientation(:,3,j) = (/  0.0_ReKi, SIN(Angle),     COS(Angle) /)
            
         Mesh1_O%TranslationVel(:,j)  = 0.0_ReKi ! (/ 1., 1.,  0. /)*.5
         Mesh1_O%RotationVel(:,j)     = 0.0_ReKi ! (/ 0., 0.5, 0.5 /)*.5
         Mesh1_O%TranslationAcc(:,j)  = 0.0_ReKi ! (/ 1., 1., 0. /)*.115
         Mesh1_O%RotationAcc(:,j)     = 0.0_ReKi ! (/ 1., 1., 1. /)*.115
         
          Mesh1_O%TranslationDisp(:,J) = 0.0_ReKi !(/ 2., 0.,  0. /)
      
      end do
                    
               
      !.........................
      ! Mesh2 (Output: Loads)
      !.........................
                    
      NNodes = 1
   
      CALL MeshCreate(  BlankMesh       = mesh2_O           &
                        , IOS           = COMPONENT_OUTPUT  &
                        , NNodes        = NNodes            &
                        , Force         = .TRUE.            &
                        , Moment        = .TRUE.            &
                        , ErrStat       = ErrStat           &
                        , ErrMess       = ErrMsg            )   
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

      do j=1,NNodes
           
         Angle = 0
         !Angle = (-25. + j*j)*D2R  !note this "looks" like the transpose, but isn't
         Orientation(:,1) = (/ COS(Angle), -1.*SIN(Angle), 0.0_ReKi /)
         Orientation(:,2) = (/ SIN(Angle),     COS(Angle), 0.0_ReKi /)
         Orientation(:,3) = (/      0.,        0.0,        1.0 /)
      
            ! place nodes in a line
         position = (/0.25_ReKi, 0.5_ReKi, 1.0_ReKi*(j-1) + 0.25 /)
         CALL MeshPositionNode ( mesh2_O, j, position, ErrStat, ErrMsg, &
               Orient= Orientation )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   
         CALL MeshConstructElement ( mesh2_O, ELEMENT_POINT, ErrStat, ErrMsg, j )
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))            
         

      END DO
                  
         ! that's our entire mesh:
      CALL MeshCommit ( mesh2_O, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))   
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh2 output for test 9.")
         
      !..............
      ! initialize output fields:
      !..............
      do j=1,Mesh2_O%NNodes
         Mesh2_O%Force( :,j) = (/  0.6, 1.,-0.5   /)  !*(j*0.5)
         Mesh2_O%Moment(:,j) = (/  4.0, 0.1, 0.3  /)
      end do   
   
   
   END subroutine CreateOutputMeshes_Test9

   ! ..............................................
   subroutine CreateOutputMeshes_Test10
   ! monopile case:
      real(ReKi),parameter :: z1(5)=(/ 10., -12.5000743865967, -5.00004959106445, 2.49997520446777, -20.0000991821289/)
      real(ReKi),parameter :: z2(6)=(/-20.0000991821289, 10., -20. ,0., -20., 0./)
   
      Mesh1Type = ELEMENT_POINT
      Mesh2Type = ELEMENT_POINT
      

      !.........................
      ! Mesh1 (Output: Motions)  i.e., subdyn's y2 mesh
      !.........................
      
      Nnodes = size(z1)
            
      CALL MeshCreate( BlankMesh          = mesh1_O           &
                       , IOS              = COMPONENT_OUTPUT  &
                       , NNodes           = NNodes            &
                       , Orientation      = .TRUE.            &
                       , TranslationDisp  = .TRUE.            &
                       , TranslationVel   = .TRUE.            &
                       , RotationVel      = .TRUE.            &
                       , TranslationAcc   = .TRUE.            &
                       , RotationAcc      = .TRUE.            &                                
                       , ErrStat          = ErrStat           &
                       , ErrMess          = ErrMsg            )
      
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))


                                    
      do j=1,NNodes 
            ! place nodes in a line
         position = (/0.0_ReKi, 0.0_ReKi, z1(j) /)
         CALL MeshPositionNode ( mesh1_O, j, position, ErrStat, ErrMsg )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                      
      END DO   

      do j=1,NNodes 
         CALL MeshConstructElement ( mesh1_O, ELEMENT_POINT, ErrStat, ErrMsg, j )
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
      end do ! 
            
      CALL MeshCommit ( mesh1_O, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh1 output for test 10.")

      !..............
      ! initialize output fields:
      !..............      
      
      ! node 1:
      Mesh1_O%Orientation(1,:,1)=(/ 1.000000000000000, -0.000000000011870, 0.000000005946981/)
      Mesh1_O%Orientation(2,:,1)=(/ 0.000000000011870,  0.999999880790710, 0.000000000735951/)
      Mesh1_O%Orientation(3,:,1)=(/-0.000000005946981, -0.000000000735951, 1.000000000000000/)
      
      ! node 2:
      Mesh1_O%Orientation(1,:,2)=(/ 1.000000000000000,  0.000000000001336, -0.000000006893778/)
      Mesh1_O%Orientation(2,:,2)=(/-0.000000000001336,  1.000000000000000,  0.000000000002210/)
      Mesh1_O%Orientation(3,:,2)=(/ 0.000000006893778, -0.000000000002210,  1.000000000000000/)                    
      
      ! node 3:
      Mesh1_O%Orientation(1,:,3)=(/ 1.000000000000000,  0.000000000001093,  0.000000012251162/)
      Mesh1_O%Orientation(2,:,3)=(/-0.000000000001093,  1.000000000000000,  0.000000000027043/)
      Mesh1_O%Orientation(3,:,3)=(/-0.000000012251162, -0.000000000027043,  1.000000000000000/)
      
      ! node 4:
      Mesh1_O%Orientation(1,:,4)=(/ 1.000000000000000,  0.000000000007391, -0.000000013338419/)
      Mesh1_O%Orientation(2,:,4)=(/-0.000000000007391,  1.000000000000000,  0.000000000005140/)
      Mesh1_O%Orientation(3,:,4)=(/ 0.000000013338419, -0.000000000005140,  1.000000000000000/)
      
      ! node 5:
      Mesh1_O%Orientation(1,:,5)=(/ 1.0, 0.0, 0.0 /)
      Mesh1_O%Orientation(2,:,5)=(/ 0.0, 1.0, 0.0 /)
      Mesh1_O%Orientation(3,:,5)=(/ 0.0, 0.0, 1.0 /)
   
      
      do i=1,4
         print *, 'DCM * DCM^T ', i
         call wrmatrix( MATMUL( Mesh1_O%Orientation(:,:,1) , TRANSPOSE( Mesh1_O%Orientation(:,:,1) ) ),cu,'f15.5' )
      END DO
                  
      
      Mesh1_O%TranslationDisp(:,1) = (/ 2.37597159724601e-08, 2.38855668577287e-09,-0.000899768609087914 /)
      Mesh1_O%TranslationDisp(:,2) = (/-4.48765788974015e-08,-2.53428389385135e-11,-0.000226117306738161 /)
      Mesh1_O%TranslationDisp(:,3) = (/-6.85023735513823e-08,-1.91699989215977e-11,-0.000451665051514283 /)
      Mesh1_O%TranslationDisp(:,4) = (/-7.33556220211540e-08,-1.62197810738007e-10,-0.000676169584039599 /)
      Mesh1_O%TranslationDisp(:,5) = (/ 0, 0, 0 /)
            
      Mesh1_O%TranslationVel(:,1) = (/ 9.30248279473744e-05, 9.57865177042550e-06, 0.000917342957109213 /)
      Mesh1_O%TranslationVel(:,2) = (/-0.000179907845449634,-6.50597939966247e-08,-0.00444589508697391 /)
      Mesh1_O%TranslationVel(:,3) = (/-0.000275392638286576, 1.26055965665728e-09,-0.00660990970209241/)
      Mesh1_O%TranslationVel(:,4) = (/-0.000295042205834761,-5.36239895154722e-07,-0.00463324552401900/)
      Mesh1_O%TranslationVel(:,5) = (/ 0, 0, 0 /)
                    									                
      Mesh1_O%RotationVel(:,1) = (/2.96345683636901e-06,-2.34051149163861e-05,-5.91553508400011e-08 /)
      Mesh1_O%RotationVel(:,2) = (/4.80349626741372e-09, 2.67697760136798e-05, 2.29704255616525e-09 /)
      Mesh1_O%RotationVel(:,3) = (/1.04646687759669e-07,-4.83136973343790e-05,-1.12983400413214e-09 /)
      Mesh1_O%RotationVel(:,4) = (/2.53845655606710e-08, 5.26454241480678e-05, 2.01115391007534e-08/)
      Mesh1_O%RotationVel(:,5) = (/ 0, 0, 0 /)      
            
      Mesh1_O%TranslationAcc(:,1) = (/ 0.172380834817886, 0.0190481301397085,    1.77754533290863 /)
      Mesh1_O%TranslationAcc(:,2) = (/-0.362261772155762,-5.16042709932663e-05, -8.71160507202148 /)
      Mesh1_O%TranslationAcc(:,3) = (/-0.560356795787811, 0.000119802658446133,-12.8353881835938 /)
      Mesh1_O%TranslationAcc(:,4) = (/-0.600724101066589,-0.000828713586088270, -8.91974353790283 /)
      Mesh1_O%TranslationAcc(:,5) = (/ 0, 0, 0 /)
      
      Mesh1_O%RotationAcc(:,1) = (/ 0.00601908657699823, -0.0441803149878979, -0.000328276102663949/)
      Mesh1_O%RotationAcc(:,2) = (/ 4.16704733652296e-06, 0.0473403967916966,  4.39332507085055e-06/)
      Mesh1_O%RotationAcc(:,3) = (/ 0.000204067997401580,-0.0911907851696014, -1.95350694411900e-05/)
      Mesh1_O%RotationAcc(:,4) = (/ 7.36967340344563e-05, 0.0997802764177322,  7.97090106061660e-05/)
      Mesh1_O%RotationAcc(:,5) = (/ 0, 0, 0 /)
      		
               
      
! test:
CALL EYE(Mesh1_O%Orientation,ErrStat,ErrMsg)
Mesh1_O%RotationAcc = 0
Mesh1_O%TranslationAcc = 0
Mesh1_O%RotationVel = 0      
Mesh1_O%TranslationVel = 0
Mesh1_O%TranslationDisp = 0

      !.........................
      ! Mesh2 (Output: Loads) i.e., HydroDyn's morision lumped mesh
      !.........................
                    
      NNodes = size(z2)
   
      CALL MeshCreate(  BlankMesh       = mesh2_O           &
                        , IOS           = COMPONENT_OUTPUT  &
                        , NNodes        = NNodes            &
                        , Force         = .TRUE.            &
                        , Moment        = .TRUE.            &
                        , ErrStat       = ErrStat           &
                        , ErrMess       = ErrMsg            )   
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

      do j=1,NNodes
                 
            ! place nodes in a line
         position =      (/0.0_ReKi, 0.0_ReKi, z2(j) /)       
         CALL MeshPositionNode ( mesh2_O, j, position, ErrStat, ErrMsg )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   
         CALL MeshConstructElement ( mesh2_O, ELEMENT_POINT, ErrStat, ErrMsg, j )
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     

      END DO
                  
         ! that's our entire mesh:
      CALL MeshCommit ( mesh2_O, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))   
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh2 output for test 10.")
         
      !..............
      ! initialize output fields:
      !..............
      Mesh2_O%Moment = 0.0
      Mesh2_O%Force  = 0.0
      Mesh2_O%Force(3,3:6) = (/5431999., 274118.281250000, -5431999., -274118.281250000/) 
      !Mesh2_O%Force(3,3:6) = (/5431999., 274118.281250000, -5431999., -274118.281250000/) / 1000.0_ReKi
         
   
   END subroutine CreateOutputMeshes_Test10
   ! ..............................................
   subroutine CreateOutputMeshes_Test11

      
   
      Mesh1Type = ELEMENT_Line2
      Mesh2Type = ELEMENT_Line2
            
      !.........................
      ! Mesh1 (Output: Motions)
      !.........................
      
      Nnodes = 10
            
      CALL MeshCreate( BlankMesh          = mesh1_O           &
                       , IOS              = COMPONENT_OUTPUT  &
                       , NNodes           = NNodes            &
                       , Orientation      = .TRUE.            &
                       , TranslationDisp  = .TRUE.            &
                       , TranslationVel   = .TRUE.            &
                       , RotationVel      = .TRUE.            &
                       , TranslationAcc   = .TRUE.            &
                       , RotationAcc      = .TRUE.            &                                
                       , ErrStat          = ErrStat           &
                       , ErrMess          = ErrMsg            )
      
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

                                    
      do j=1,NNodes 
            ! place nodes in a line
         position =  (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi /)           
         CALL MeshPositionNode ( mesh1_O, j, position, ErrStat, ErrMsg )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                      
      END DO   
      mesh1_O%Position(3,:) = (/0., 0.5, 1.5, 2., 4., 4.1, 5., 8., 8.5, 10. /)

      do j=2,NNodes 
         CALL MeshConstructElement ( mesh1_O, ELEMENT_LINE2, ErrStat, ErrMsg, j-1, j )
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
      end do ! 
            
      CALL MeshCommit ( mesh1_O, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh1 output for test 11.")

      !..............
      ! initialize output fields:
      !..............      
      
      do j=1,Mesh1_O%NNodes
      
         Angle = 0      
   !      Angle = (20*j)*D2R      
    !     Angle = (20)*D2R      
         !note this "looks" like the transpose, but isn't
         Mesh1_O%Orientation(:,1,j) = (/  1.0,       0.0 ,            0.0 /)
         Mesh1_O%Orientation(:,2,j) = (/  0.0_ReKi, COS(Angle), -1.*SIN(Angle) /)
         Mesh1_O%Orientation(:,3,j) = (/  0.0_ReKi, SIN(Angle),     COS(Angle) /)
            
         Mesh1_O%TranslationVel(:,j)  = (/ 1., 0.,  0. /)
         Mesh1_O%RotationVel(:,j)     = 0.0_ReKi ! (/ 0., 0.5, 0.5 /)*.5
         Mesh1_O%TranslationAcc(:,j)  = 0.0_ReKi ! (/ 1., 1., 0. /)*.115
         Mesh1_O%RotationAcc(:,j)     = 0.0_ReKi ! (/ 1., 1., 1. /)*.115
         
          Mesh1_O%TranslationDisp(:,J) = 0.0_ReKi !(/ 2., 0.,  0. /)
      
      end do
                    
               
      !.........................
      ! Mesh2 (Output: Loads)
      !.........................
                    
      NNodes = 7   
      CALL MeshCreate(  BlankMesh       = mesh2_O           &
                        , IOS           = COMPONENT_OUTPUT  &
                        , NNodes        = NNodes            &
                        , Force         = .TRUE.            &
                        , Moment        = .TRUE.            &
                        , ErrStat       = ErrStat           &
                        , ErrMess       = ErrMsg            )   
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

      do j=1,NNodes
           
         Angle = 0
         !Angle = (-25. + j*j)*D2R  !note this "looks" like the transpose, but isn't
         Orientation(:,1) = (/ COS(Angle), -1.*SIN(Angle), 0.0_ReKi /)
         Orientation(:,2) = (/ SIN(Angle),     COS(Angle), 0.0_ReKi /)
         Orientation(:,3) = (/      0.,        0.0,        1.0 /)
      
            ! place nodes in a line
         position = (/0._ReKi, 0._ReKi, (j-1)*10.0_ReKi/(Nnodes-1) /)
         CALL MeshPositionNode ( mesh2_O, j, position, ErrStat, ErrMsg, &
               Orient= Orientation )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))            

      END DO

      
      
      do j=2,NNodes 
         CALL MeshConstructElement ( mesh2_O, ELEMENT_LINE2, ErrStat, ErrMsg, j-1, j )
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
      end do !       
      
                  
         ! that's our entire mesh:
      CALL MeshCommit ( mesh2_O, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))   
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh2 output for test 11.")
         
      !..............
      ! initialize output fields:
      !..............
      do j=1,Mesh2_O%NNodes
         Mesh2_O%Force( :,j) = (/  0.0, 1.0, 0.0  /) 
         Mesh2_O%Moment(:,j) = (/  0.0, 0.0, 0.0  /)
      end do   
   
   
   END subroutine CreateOutputMeshes_Test11
   ! ..............................................
   subroutine CreateOutputMeshes_Test13
      real(reki),parameter  :: a = 1.
      real(reki),parameter  :: omega=2.
      real(reki)            :: L, LineLen
      real(reKi)            :: x, dx, yp, omega2, a2

      real(reKi),parameter  :: xAry1_disp(11) = (/ &
      0.00000000000000_DbKi, 0.51960358494079_DbKi, 1.05066474308880_DbKi, 1.58938188333378_DbKi, 2.12044304148178_DbKi, 2.64004662642258_DbKi, 3.15959741148684_DbKi, 3.69065856963484_DbKi, 4.22937570987983_DbKi, 4.76043686802783_DbKi, 5.27998765309209_DbKi /)
      
      real(reKi),parameter  :: xAry1(11) = (/ &
      0.00000000000000_DbKi, 0.24325351916746_DbKi, 0.55141234255808_DbKi, 1.01941540016335_DbKi, 1.32757422355397_DbKi, 1.57082774272143_DbKi, 1.81404984596235_DbKi, 2.12220866935298_DbKi, 2.59021172695825_DbKi, 2.89837055034887_DbKi, 3.14159265358979_DbKi  /)     
      real(reKi),parameter  :: xAry2(31) = (/ &
      0.00000000000000_DbKi, 0.07885397560510_DbKi, 0.15927874753700_DbKi, 0.24325351916746_DbKi, 0.33329156461934_DbKi, 0.43357120212193_DbKi, 0.55141234255808_DbKi, 0.69925569283602_DbKi, 0.87157204988542_DbKi, 1.01941540016335_DbKi, 1.13725654059951_DbKi, 1.23753617810209_DbKi, 1.32757422355397_DbKi, 1.41154899518443_DbKi, 1.49197376711633_DbKi, 1.57082774272143_DbKi, 1.64965030240000_DbKi, 1.73007507433190_DbKi, 1.81404984596235_DbKi, 1.90408789141424_DbKi, 2.00436752891682_DbKi, 2.12220866935298_DbKi, 2.27005201963091_DbKi, 2.44236837668031_DbKi, 2.59021172695825_DbKi, 2.70805286739440_DbKi, 2.80833250489699_DbKi, 2.89837055034887_DbKi, 2.98234532197933_DbKi, 3.06277009391123_DbKi, 3.14159265358979_DbKi  /)
      
      Mesh1Type = ELEMENT_Line2
      Mesh2Type = ELEMENT_Line2
                        
      
      !.........................
      ! Mesh1 (Output: Motions)
      !.........................
      
      Nnodes = size(xAry1) !11
            
      CALL MeshCreate( BlankMesh          = mesh1_O           &
                       , IOS              = COMPONENT_OUTPUT  &
                       , NNodes           = NNodes            &
                       , Orientation      = .TRUE.            &
                       , TranslationDisp  = .TRUE.            &
                       , TranslationVel   = .TRUE.            &
                       , RotationVel      = .TRUE.            &
                       , TranslationAcc   = .TRUE.            &
                       , RotationAcc      = .TRUE.            &                                
                       , ErrStat          = ErrStat           &
                       , ErrMess          = ErrMsg            )
      
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

                             
      L = TwoPi/omega
      dx = L/(NNodes-1.0_ReKi)   
      do j=1,NNodes 
            ! place nodes in a line
         x        = xAry1(j) !(j-1.0_ReKi) * dx
         yp       = a*omega*cos(omega*x)
         position = (/x, a*sin(omega*x), 0.0_ReKi /)                    
         Angle    = atan(yp)
         
         Orientation(:,1) = (/ COS(Angle), -1.*SIN(Angle), 0.0_ReKi /)
         Orientation(:,2) = (/ SIN(Angle),     COS(Angle), 0.0_ReKi /)
         Orientation(:,3) = (/      0.,        0.0,        1.0 /)
!call wrmatrix(orientation,cu,'f15.5','mesh1 '//trim(num2lstr(j))//' '//trim(num2lstr(angle)))         
                  
         CALL MeshPositionNode ( mesh1_O, j, position, ErrStat, ErrMsg, orient=Orientation )     
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                      
      END DO   

      do j=2,NNodes 
         CALL MeshConstructElement ( mesh1_O, ELEMENT_LINE2, ErrStat, ErrMsg, j-1, j )
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
      end do ! 
            
      CALL MeshCommit ( mesh1_O, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh1 output for test 13.")

!      LineLen = 0;
!      do j=1,mesh1_O%ElemTable(ELEMENT_Line2)%nelem
!         LineLen = LineLen + mesh1_O%ElemTable(ELEMENT_Line2)%Elements(j)%det_jac*2
!      end do
!      
!print *, 'length of line is =', LineLen      
LineLen = 5.27038      

      !..............
      ! initialize output fields:
      !..............      
      
      omega2 = 1.19 !0.5*omega
      a2     = 0.25*a !2*a !
      L = TwoPi/omega2
      dx = L/(NNodes-1.0_ReKi)   

      do j=1,Mesh1_O%NNodes
         
         x        = xAry1_disp(j) !(j-1.0_ReKi) * dx 
         yp       = a2*omega2*cos(omega2*x)
         position = (/x, a2*sin(omega2*x), 0.0_ReKi /)             
         Mesh1_O%TranslationDisp(:,J) = position - Mesh1_O%Position(:,j)
         
!print *,j,  Mesh1_O%TranslationDisp(:,J)

         Angle = atan(yp)
                  
         Mesh1_O%Orientation(:,1,j) = (/ COS(Angle), -1.*SIN(Angle), 0.0_ReKi /)
         Mesh1_O%Orientation(:,2,j) = (/ SIN(Angle),     COS(Angle), 0.0_ReKi /)
         Mesh1_O%Orientation(:,3,j) = (/      0.,        0.0,        1.0 /)
         
         Mesh1_O%TranslationVel(:,j)  = (/ 0., 0.,  0. /)
         Mesh1_O%RotationVel(:,j)     = 0.0_ReKi ! (/ 0., 0.5, 0.5 /)*.5
         Mesh1_O%TranslationAcc(:,j)  = 0.0_ReKi ! (/ 1., 1., 0. /)*.115
         Mesh1_O%RotationAcc(:,j)     = 0.0_ReKi ! (/ 1., 1., 1. /)*.115
         
      
      end do
       
!Mesh1_O%TranslationDisp = 0.0_ReKi         
!Mesh1_O%Orientation = Mesh1_O%RefOrientation         
         
      
               
      !.........................
      ! Mesh2 (Output: Loads)
      !.........................
                    
      NNodes = size(xAry2) !31 !19   
      L = TwoPi/omega      
      dx = L/(NNodes-1.0_ReKi)   
      CALL MeshCreate(  BlankMesh       = mesh2_O           &
                        , IOS           = COMPONENT_OUTPUT  &
                        , NNodes        = NNodes            &
                        , Force         = .TRUE.            &
                        , Moment        = .TRUE.            &
                        , ErrStat       = ErrStat           &
                        , ErrMess       = ErrMsg            )   
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

            
            
      do j=1,NNodes
            ! place nodes in a line
         x           = xAry2(j) !(j-1.0_ReKi) * dx
         yp          = a*omega*cos(omega*x)
         position    = (/x, a*sin(omega*x), 0.0_ReKi /)                    
         Angle       = atan(yp) 
         
         Orientation(:,1) = (/ COS(Angle), -1.*SIN(Angle), 0.0_ReKi /)
         Orientation(:,2) = (/ SIN(Angle),     COS(Angle), 0.0_ReKi /)
         Orientation(:,3) = (/      0.,        0.0,        1.0 /)
                 
         CALL MeshPositionNode ( mesh2_O, j, position, ErrStat, ErrMsg, Orient= Orientation )    
         
         !CALL MeshPositionNode ( mesh2_O, j, position, ErrStat, ErrMsg )     
         
         
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))            

      END DO
            
      do j=2,NNodes 
         CALL MeshConstructElement ( mesh2_O, ELEMENT_LINE2, ErrStat, ErrMsg, j-1, j )
         IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
      end do !       
      
                  
         ! that's our entire mesh:
      CALL MeshCommit ( mesh2_O, ErrStat, ErrMsg )   
      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))   
      if (ErrStat >= AbortErrLev) CALL ProgAbort("Error creating Mesh2 output for test 11.")
         
      !..............
      ! initialize output fields:
      !..............
      do j=1,Mesh2_O%NNodes
         Mesh2_O%Force( :,j) = (/  0.0, 0.5, 0.0  /) 
         Mesh2_O%Moment(:,j) = (/  0.0, 0.0, 0.0  /)
      end do   
   
   
   END subroutine CreateOutputMeshes_Test13  
   
   subroutine InitTest13Loads()
      !..............
      ! initialize output fields:
      !..............
      do j=1,Mesh2_O%NNodes         
         Mesh2_O%Force( :,j) = Mesh2_I%Orientation(2,:,j)*0.5 ! (/ 0.0_reKi, Mesh2_I%Orientation(2,2,j)*0.5, 0.0_ReKi /)
         Mesh2_O%Moment(:,j) = 0.0 ! (/  0.0, 0.0, 0.5  /)*(-j*0.0)
      end do   
   end subroutine InitTest13Loads      
   
   ! ..............................................
   SUBROUTINE TestOrientations()
   
      real(reKi)              :: DCM(3,3)
      real(reKi)              :: angles(3), da, astart
      !REAL(ReKi), parameter   :: da = 4.0   ! degrees
      integer                 :: k
   
      NNodes = 15
      astart = -180.0_ReKi
      
      da=abs(astart)*2.0_ReKi/(NNodes-1)
                   
      do i=1,NNodes 
         angles(1) = (astart + da*(i-1) )*D2R  
         do j = 1,NNodes 
            angles(2) = (astart + da*(j-1) )*D2R  
            do k = 1,NNodes 
               angles(3) = (astart + da*(k-1) )*D2R  
            
         
               CALL SmllRotTrans( 'Mesh Orientation Extrapolation', angles(1), angles(2),  angles(3), DCM, 'angle='//trim(num2lstr(angles(1))), ErrStat, ErrMsg )                 
write(75,'(I10,2(1x,I10),15(1x,F15.5))') i,j,k, angles, DCM

               DCM = DCM_exp ( angles )
write(76,'(I10,2(1x,I10),15(1x,F15.5))') i,j,k, angles, DCM
            end do
         end do         
      END DO   
      
      close(75)
      close(76)
   
   END SUBROUTINE TestOrientations
   
   function GetDCM( angle, dim )
   
      REAL(ReKi),     intent(in) :: angle
      REAL(ReKi)                 :: GetDCM(3,3)
      integer(IntKi), intent(in) :: dim
      
                  
      ! note that these look like the transpose of the matrix, but each line is a column      
      select case (dim)
      case( 1 )
         GetDCM(:,1) = (/  1.0,       0.0 ,            0.0 /)
         GetDCM(:,2) = (/  0.0_ReKi, COS(Angle), -1.*SIN(Angle) /)
         GetDCM(:,3) = (/  0.0_ReKi, SIN(Angle),     COS(Angle) /)
      case (2)
         GetDCM(:,1) = (/  COS(Angle),   0.0_ReKi,  1.*SIN(Angle) /)
         GetDCM(:,2) = (/  0.0,               1.0,            0.0 /)
         GetDCM(:,3) = (/-1.*SIN(Angle), 0.0_ReKi,     COS(Angle) /)
      case (3)      
         Orientation(:,1) = (/ COS(Angle), -1.*SIN(Angle), 0.0_ReKi /)
         Orientation(:,2) = (/ SIN(Angle),     COS(Angle), 0.0_ReKi /)
         Orientation(:,3) = (/      0.,        0.0,        1.0 /)
         
      end select
            
      
   end function GetDCM
   ! ..............................................

   subroutine TransferMotionData()
   
         IF (Mesh1Type == ELEMENT_POINT ) THEN
            IF ( Mesh2Type == ELEMENT_POINT ) THEN                        
               CALL Transfer_Point_to_Point( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     
            ELSEIF ( Mesh2Type == ELEMENT_LINE2) THEN                                                                
               CALL Transfer_Point_to_Line2( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );
            END IF                                                                                                   
         ELSEIF ( Mesh1Type == ELEMENT_LINE2 ) THEN                                                                  
            IF ( Mesh2Type == ELEMENT_LINE2 ) THEN                                                                   
               CALL Transfer_Line2_to_Line2( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );
            ELSEIF ( Mesh2Type == ELEMENT_POINT ) THEN        
               CALL Transfer_Line2_to_Point( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg )
            END IF            
         END IF         
         IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))
   
   end subroutine TransferMotionData
   
   subroutine TransferLoadData()
   
         IF (Mesh1Type == ELEMENT_POINT ) THEN
            IF ( Mesh2Type == ELEMENT_POINT ) THEN                        
               CALL Transfer_Point_to_Point( Mesh2_O, Mesh1_I, Map_Mod2_Mod1, ErrStat, ErrMsg, Mesh2_I, Mesh1_O )
            ELSEIF ( Mesh2Type == ELEMENT_LINE2) THEN                                                                
               CALL Transfer_Line2_to_Point( Mesh2_O, Mesh1_I, Map_Mod2_Mod1, ErrStat, ErrMsg, Mesh2_I, Mesh1_O )
            END IF                                                                                                   
         ELSEIF ( Mesh1Type == ELEMENT_LINE2 ) THEN                                                                  
            IF ( Mesh2Type == ELEMENT_LINE2 ) THEN                                                                   
               CALL Transfer_Line2_to_Line2( Mesh2_O, Mesh1_I, Map_Mod2_Mod1, ErrStat, ErrMsg, Mesh2_I, Mesh1_O )
            ELSEIF ( Mesh2Type == ELEMENT_POINT ) THEN        
               CALL Transfer_Point_to_Line2( Mesh2_O, Mesh1_I, Map_Mod2_Mod1, ErrStat, ErrMsg, Mesh2_I, Mesh1_O )
            END IF            
         END IF         
         IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))
   
   end subroutine TransferLoadData
   
end module TestMeshMapping_Mod

