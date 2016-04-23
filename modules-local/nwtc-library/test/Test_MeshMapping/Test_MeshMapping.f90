PROGRAM Test_TestMeshMapping

   USE NWTC_Library
   IMPLICIT NONE   

   TYPE(meshtype) :: mesh1_I, mesh1_O
   TYPE(meshtype) :: mesh2_I, mesh2_O 

   TYPE(meshtype) :: mesh1_I_op, mesh1_O_op
   TYPE(meshtype) :: mesh2_I_op, mesh2_O_op 
   
   
   TYPE(MeshMapType)       :: Map_Mod1_Mod2        ! Data for mapping meshes from mod1 to mod2
   TYPE(MeshMapType)       :: Map_Mod2_Mod1        ! Data for mapping meshes from mod1 to mod2
   
   REAL(R8Ki)              :: Orientation(3,3)
   REAL(ReKi)              :: position(3)
   REAL(ReKi)              :: Angle
   REAL(ReKi), allocatable :: LinVec_1(:), LinVec_1_theta(:), LinVec_1_a(:), LinVec_1_b(:), LinVec_1_c(:)
   REAL(ReKi), allocatable :: LinVec_2(:), LinVec_2_theta(:) ! data for linearization testing
   
   CHARACTER(*), PARAMETER :: Fmt = '(ES10.3E2)'
   
   !
   INTEGER :: NNodes, I,J, n1, n2
      
   INTEGER :: Mesh1Type
   INTEGER :: Mesh2Type
   
   INTEGER(IntKi) :: ErrStat
   CHARACTER(1024) :: ErrMsg   
   
   logical, parameter :: TestLinearization = .true.
   
   INTEGER :: TestNumber 
   CHARACTER(256) :: BinOutputName 
      

   call NWTC_Init()
   
   DO TestNumber=1,13 !1,9

      print *, '---------------------------------------------------------------'
      print *, '   Test ', TestNumber
      print *, '---------------------------------------------------------------'
      ! ..............................................................................................................................   
      ! Mesh1 fields:
      !   Mesh1_I (input) has loads
      !   Mesh1_O (output) has motions
      ! ..............................................................................................................................   

      ! ..............................................................................................................................   
      ! Mesh2 fields:
      !   Mesh2_I (input) has motions
      !   Mesh2_O (output) has loads
      ! ..............................................................................................................................   
      
   
   
      ! ..............................................................................................................................   
      ! Create output meshes: 
      !   Mesh1_O (output) has motions
      !   Mesh2_O (output) has loads
      ! ..............................................................................................................................   
      
      ! These subroutines will set the Mesh1Type and Mesh2Type variables, then create the output meshes, set reference  
      ! position/orientation, and set initial outputs for each module (to test mapping to inputs)
      SELECT CASE ( TestNumber )
      CASE(1) ! 1 point to 5 points
         CALL CreateOutputMeshes_Test1()
      CASE(2) ! 'T' with resolution gain
         CALL CreateOutputMeshes_Test2('A') ! was 'A'
      CASE(3) ! 'T' with loss of resolution
         CALL CreateOutputMeshes_Test2('B')
      CASE(4) ! 'T' with equal nodes
         CALL CreateOutputMeshes_Test2('C')
      CASE(5)
         CALL CreateOutputMeshes_Test5()
      CASE(6)
         CALL CreateOutputMeshes_Test6()
      CASE(7)
         CALL CreateOutputMeshes_Test7('A')
      CASE(8)
         CALL CreateOutputMeshes_Test7('B')
      CASE(9)
         CALL CreateOutputMeshes_Test9()
      CASE(10) !monopile (point-to-point)
         CALL CreateOutputMeshes_Test10()
      CASE(11) 
         CALL CreateOutputMeshes_Test11()
      CASE(12) 
         CALL CreateOutputMeshes_Test5Orient
      CASE(13) 
         CALL CreateOutputMeshes_Test13
      END SELECT
   
      WRITE(BinOutputName,'(A,A,A)') 'Test', TRIM(Num2LStr(TestNumber)),'Meshes.bin'
      
   
      ! ..............................................................................................................................   
      ! Create sibling input meshes: 
      !   Mesh1_I (input) has loads
      !   Mesh2_I (input) has motions
      ! ..............................................................................................................................   
   
   
      CALL MeshCopy (        SrcMesh      = mesh1_O                &
                           , DestMesh     = mesh1_I                &
                           , CtrlCode     = MESH_SIBLING           &
                           , IOS          = COMPONENT_INPUT        &
                           ,Force         = .TRUE.                 &
                           ,Moment        = .TRUE.                 &
                           ,ErrStat       = ErrStat                &
                           ,ErrMess       = ErrMsg                 )   
               IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   
      !..............................
   
      CALL MeshCopy (        SrcMesh          = mesh2_O          &
                           , DestMesh         = mesh2_I          &
                           , CtrlCode         = MESH_SIBLING     &
                           , IOS              = COMPONENT_INPUT  &
                           , Orientation      = .TRUE.           &
                           , TranslationDisp  = .TRUE.           &
                           , TranslationVel   = .TRUE.           &
                           , RotationVel      = .TRUE.           &
                           , TranslationAcc   = .TRUE.           &
                           , RotationAcc      = .TRUE.           &
                           , ErrStat          = ErrStat          &
                           , ErrMess          = ErrMsg           )
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
    
      ! ..............................................................................................................................   
      ! Initialize the mapping data: 
      ! ..............................................................................................................................   

      CALL MeshMapCreate( Mesh1_O, Mesh2_I,     Map_Mod1_Mod2, ErrStat, ErrMsg );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))   
      CALL MeshMapCreate( Mesh2_O, Mesh1_I,     Map_Mod2_Mod1, ErrStat, ErrMsg );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      
      ! ..............................................................................................................................   
      ! Map the outputs to inputs:
      ! ..............................................................................................................................   
      IF (Mesh1Type == ELEMENT_POINT ) THEN
         IF ( Mesh2Type == ELEMENT_POINT ) THEN            
            CALL Transfer_Point_to_Point( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))     
            CALL Transfer_Point_to_Point( Mesh2_O, Mesh1_I, Map_Mod2_Mod1, ErrStat, ErrMsg, Mesh2_I, Mesh1_O );   IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))               
         ELSEIF ( Mesh2Type == ELEMENT_LINE2) THEN                                                                
            CALL Transfer_Point_to_Line2( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))     
            CALL Transfer_Line2_to_Point( Mesh2_O, Mesh1_I, Map_Mod2_Mod1, ErrStat, ErrMsg, Mesh2_I, Mesh1_O );   IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))         
         END IF                                                                                                   
      ELSEIF ( Mesh1Type == ELEMENT_LINE2 ) THEN                                                                  
         IF ( Mesh2Type == ELEMENT_LINE2 ) THEN                                                                   
            CALL Transfer_Line2_to_Line2( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))   
            if (TestNumber == 5 .or. TestNumber == 12 ) call InitTest5Loads()
            if (TestNumber == 13 ) call InitTest13Loads()
            CALL Transfer_Line2_to_Line2( Mesh2_O, Mesh1_I, Map_Mod2_Mod1, ErrStat, ErrMsg, Mesh2_I, Mesh1_O );   IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))        
         ELSEIF ( Mesh2Type == ELEMENT_POINT ) THEN        
            CALL Transfer_Line2_to_Point( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))                  
            CALL Transfer_Point_to_Line2( Mesh2_O, Mesh1_I, Map_Mod2_Mod1, ErrStat, ErrMsg, Mesh2_I, Mesh1_O );   IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))     
         END IF            
      END IF

      
      ! ..............................................................................................................................   
      ! Write results to file(s)
      ! ..............................................................................................................................   
      call WriteMappingTransferToFile(Mesh1_I, Mesh1_O, Mesh2_I, Mesh2_O, Map_Mod1_Mod2, Map_Mod2_Mod1, BinOutputName)                       
      
      if (TestNumber == 13) then

         !write( 17, '(A)' )
         DO i=1,mesh2_o%nnodes
            n1 = mesh1_o%ElemTable(ELEMENT_LINE2)%Elements(Map_Mod1_Mod2%MapMotions(i)%OtherMesh_Element)%ElemNodes(1)
            n2 = mesh1_o%ElemTable(ELEMENT_LINE2)%Elements(Map_Mod1_Mod2%MapMotions(i)%OtherMesh_Element)%ElemNodes(2)
      
            WRITE (17, '(6(F15.5))' ) atan2( mesh2_i%refOrientation(1,2,i ) , mesh2_i%refOrientation(1,1,i ) ),&
                                      atan2( mesh2_i%Orientation(   1,2,i ) , mesh2_i%Orientation(   1,1,i ) ),&
                                      atan2( mesh1_o%refOrientation(1,2,n1) , mesh1_o%refOrientation(1,1,n1) ),&
                                      atan2( mesh1_o%Orientation(   1,2,n1) , mesh1_o%Orientation(   1,1,n1) ),&      
                                      atan2( mesh1_o%refOrientation(1,2,n2) , mesh1_o%refOrientation(1,1,n2) ),&
                                      atan2( mesh1_o%Orientation(   1,2,n2) , mesh1_o%Orientation(   1,1,n2) )
         end do                              
   


      end if
     
if ( TestLinearization ) then      
      ! ..............................................................................................................................   
      ! Linearize about the operating points:
      ! ..............................................................................................................................   
      
      call AllocAry( LinVec_1, Mesh1_O%nnodes*3, 'LinVec_1', ErrStat, ErrMsg);       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      call AllocAry( LinVec_2, Mesh2_O%nnodes*3, 'LinVec_2', ErrStat, ErrMsg);       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      call AllocAry( LinVec_1_theta, Mesh1_O%nnodes*3, 'LinVec_1_theta', ErrStat, ErrMsg);       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      call AllocAry( LinVec_2_theta, Mesh2_O%nnodes*3, 'LinVec_2_theta', ErrStat, ErrMsg);       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      call AllocAry( LinVec_1_a, Mesh1_O%nnodes*3, 'LinVec_1_a', ErrStat, ErrMsg);       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      call AllocAry( LinVec_1_b, Mesh1_O%nnodes*3, 'LinVec_1_b', ErrStat, ErrMsg);       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      call AllocAry( LinVec_1_c, Mesh2_O%nnodes*3, 'LinVec_1_c', ErrStat, ErrMsg);       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      
      
      CALL MeshCopy (        SrcMesh          = mesh1_O          &
                           , DestMesh         = mesh1_O_op       &
                           , CtrlCode         = MESH_NEWCOPY     &
                           , ErrStat          = ErrStat          &
                           , ErrMess          = ErrMsg           )
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      
      CALL MeshCopy (        SrcMesh          = mesh2_O          &
                           , DestMesh         = mesh2_O_op       &
                           , CtrlCode         = MESH_NEWCOPY     &
                           , ErrStat          = ErrStat          &
                           , ErrMess          = ErrMsg           )
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      
            
      CALL MeshCopy (        SrcMesh          = mesh1_I          &
                           , DestMesh         = mesh1_I_op       &
                           , CtrlCode         = MESH_NEWCOPY     &
                           , ErrStat          = ErrStat          &
                           , ErrMess          = ErrMsg           )
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      
      CALL MeshCopy (        SrcMesh          = mesh2_I          &
                           , DestMesh         = mesh2_I_op       &
                           , CtrlCode         = MESH_NEWCOPY     &
                           , ErrStat          = ErrStat          &
                           , ErrMess          = ErrMsg           )
            IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
            
      
write(99, *) '---------------------------------------------------------------'
write(99, *) '   Test ', TestNumber
write(99, *) '---------------------------------------------------------------'
      
      IF (Mesh1Type == ELEMENT_POINT ) THEN
         IF ( Mesh2Type == ELEMENT_POINT ) THEN            
      
            CALL Linearize_Point_to_Point( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))     
            CALL Linearize_Point_to_Point( Mesh2_O, Mesh1_I, Map_Mod2_Mod1, ErrStat, ErrMsg, Mesh2_I, Mesh1_O );   IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))   
                     
            call wrmatrix(map_mod1_mod2%dm%mi,  99,Fmt,'Mi')
            call wrmatrix(map_mod1_mod2%dm%fx_p,99,Fmt,'fx_p')
            call wrmatrix(map_mod1_mod2%dm%tv,  99,Fmt,'tv')
            call wrmatrix(map_mod1_mod2%dm%ta1, 99,Fmt,'ta1')
            call wrmatrix(map_mod1_mod2%dm%ta2, 99,Fmt,'ta2')
            

            call wrmatrix(map_mod2_mod1%dm%li,99,Fmt,'li')
            call wrmatrix(map_mod2_mod1%dm%M1,99,Fmt,'m1')
            call wrmatrix(map_mod2_mod1%dm%M2,99,Fmt,'m2')      
            
            
         ELSEIF ( Mesh2Type == ELEMENT_LINE2) THEN                                                                
            !CALL Linearize_Point_to_Line2( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))     
            !CALL Linearize_Line2_to_Point( Mesh2_O, Mesh1_I, Map_Mod2_Mod1, ErrStat, ErrMsg, Mesh2_I, Mesh1_O );   IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))         
         END IF                                                                                                   
      ELSEIF ( Mesh1Type == ELEMENT_LINE2 ) THEN                                                                  
         IF ( Mesh2Type == ELEMENT_LINE2 ) THEN                                                                   
            !CALL Linearize_Line2_to_Line2( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))   
            !if (TestNumber == 5 .or. TestNumber == 12 ) call InitTest5Loads()
            !if (TestNumber == 13 ) call InitTest13Loads()
            !CALL Linearize_Line2_to_Line2( Mesh2_O, Mesh1_I, Map_Mod2_Mod1, ErrStat, ErrMsg, Mesh2_I, Mesh1_O );   IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))        
         ELSEIF ( Mesh2Type == ELEMENT_POINT ) THEN        
            !CALL Linearize_Line2_to_Point( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))                  
            !CALL Linearize_Point_to_Line2( Mesh2_O, Mesh1_I, Map_Mod2_Mod1, ErrStat, ErrMsg, Mesh2_I, Mesh1_O );   IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))     
         END IF            
      END IF
      
      
      ! src=Mesh1_O, dest=Mesh2_I, Map_Mod1_Mod2
   if (Mesh1Type == ELEMENT_POINT  .and.  Mesh2Type == ELEMENT_POINT ) THEN            
      
      
      !........................ 
      ! get Mesh1_O_op orientation parameters, i.e., LinVec_1_theta = theta^S
      !........................
      j=1
      do i=1,Mesh1_O%nnodes
         LinVec_1_theta(j:j+2) = GetSmllRotAngs(Mesh1_O_op%Orientation(:,:,i),ErrStat,ErrMsg)
         j = j+3
         !IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))
      end do
      
      !call wrmatrix(LinVec_1_theta,  99,Fmt,'LinVec_1_theta')
         
      ! ..................
      ! Rotational displacement:
      ! ..................
      
      do n1=1,15
         
         ! perturb x^S:
         call random_number( LinVec_1 )
         LinVec_1 = LinVec_1/real(n1**3)               
         
         LinVec_2_theta = matmul( map_mod1_mod2%dm%mi, LinVec_1 ) ! approximate delta theta^D         
             
         !call wrmatrix(LinVec_1,  99,Fmt,'LinVec_1')
         !call wrmatrix(LinVec_2_theta,  99,Fmt,'LinVec_2_theta')
         
         j=1
         do i=1,Mesh1_O%nnodes

               ! Mesh1_O%Orientation = theta^S|_op + delta theta^S
            call SmllRotTrans( 'orientation', LinVec_1_theta(j  )+LinVec_1(j  ) &
                                            , LinVec_1_theta(j+1)+LinVec_1(j+1) &
                                            , LinVec_1_theta(j+2)+LinVec_1(j+2) & 
                                            , Mesh1_O%Orientation(:,:,i), ErrStat=ErrStat, ErrMsg=ErrMsg)            
                                    
            j = j+3            
            !IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))
         end do

         !M( theta^S|_op + delta theta^S )
         CALL Transfer_Point_to_Point( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))     

         
         j=1
         do i=1,Mesh2_O%nnodes
            LinVec_2(j:j+2) = GetSmllRotAngs(Mesh2_I%Orientation(:,:,i),ErrStat,ErrMsg) 
            !IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))
            
            LinVec_2(j:j+2) = LinVec_2(j:j+2) - GetSmllRotAngs(Mesh2_I_op%Orientation(:,:,i),ErrStat,ErrMsg)
            !IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))
            
            j = j+3
         end do
         !call wrmatrix(LinVec_2-LinVec_2_theta,  99,Fmt,'LinVec_2')
         
         
         print *, 'Rotational Displacement: ', n1, TwoNorm( LinVec_2 - LinVec_2_theta )
            
      end do
      
      ! ..................
      ! Translational displacement:
      ! ..................
      
      call meshcopy( Mesh1_O_op, Mesh1_O, MESH_UPDATECOPY, ErrStat, ErrMsg) 
      do n1=1,15
         
         ! perturb x^S:
         call random_number( LinVec_1 )
         LinVec_1 = LinVec_1/real(n1**3)               
         call random_number( LinVec_1_a )
         LinVec_1_a = LinVec_1_a/real(n1**3)               
         
         LinVec_2_theta = matmul( map_mod1_mod2%dm%mi, LinVec_1 ) + matmul( map_mod1_mod2%dm%fx_p, LinVec_1_a ) ! approximate delta theta^D         
                      
         j=1
         do i=1,Mesh1_O%nnodes

               ! Mesh1_O%TranslationDisp = theta^S|_op + delta theta^S
            Mesh1_O%TranslationDisp(:,i) = Mesh1_O_op%TranslationDisp(:,i) + LinVec_1(j:j+2)
            
               ! Mesh1_O%Orientation = theta^S|_op + delta theta^S
            call SmllRotTrans( 'orientation', LinVec_1_theta(j  )+LinVec_1_a(j  ) &
                                            , LinVec_1_theta(j+1)+LinVec_1_a(j+1) &
                                            , LinVec_1_theta(j+2)+LinVec_1_a(j+2) & 
                                            , Mesh1_O%Orientation(:,:,i), ErrStat=ErrStat, ErrMsg=ErrMsg)            
                                    
            j = j+3            
            !IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))
                        
         end do
         

         !M( theta^S|_op + delta theta^S )
         CALL Transfer_Point_to_Point( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))     

         ! delta theta^D = M( theta^S|_op + delta theta^S ) - x^D|_op
         j=1
         do i=1,Mesh2_O%nnodes
            LinVec_2(j:j+2) = Mesh2_I%TranslationDisp(:,i) - Mesh2_I_op%TranslationDisp(:,i)          
            j = j+3
         end do
                  
         print *, 'Translational Displacement: ', n1, TwoNorm( LinVec_2 - LinVec_2_theta )
                              
      end do      
      
      
      ! ..................
      ! Rotational velocity:
      ! ..................
      
      call meshcopy( Mesh1_O_op, Mesh1_O, MESH_UPDATECOPY, ErrStat, ErrMsg) 
      do n1=1,15
         
         ! perturb x^S:
         call random_number( LinVec_1 )
         LinVec_1 = LinVec_1/real(n1**3)               
         
         LinVec_2_theta = matmul( map_mod1_mod2%dm%mi, LinVec_1 )  ! approximate delta theta^D         
                      
         j=1
         do i=1,Mesh1_O%nnodes

               ! Mesh1_O%RotationVel = theta^S|_op + delta theta^S
            Mesh1_O%RotationVel(:,i) = Mesh1_O_op%RotationVel(:,i) + LinVec_1(j:j+2)
                                                
            j = j+3            
                        
         end do
         

         !M( theta^S|_op + delta theta^S )
         CALL Transfer_Point_to_Point( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))     

         ! delta theta^D = M( theta^S|_op + delta theta^S ) - x^D|_op
         j=1
         do i=1,Mesh2_O%nnodes
            LinVec_2(j:j+2) = Mesh2_I%RotationVel(:,i) - Mesh2_I_op%RotationVel(:,i)          
            j = j+3
         end do
                  
         print *, 'Rotational Velocity: ', n1, TwoNorm( LinVec_2 - LinVec_2_theta )                              
      end do            
      
      
      ! ..................
      ! Translational velocity:
      ! ..................
      
      call meshcopy( Mesh1_O_op, Mesh1_O, MESH_UPDATECOPY, ErrStat, ErrMsg) 
      do n1=1,15
         
         ! perturb x^S:
         call random_number( LinVec_1 )
         LinVec_1 = LinVec_1/real(n1**3)               
         call random_number( LinVec_1_a )
         LinVec_1_a = LinVec_1_a/real(n1**3)               
         call random_number( LinVec_1_b )
         LinVec_1_b = LinVec_1_b/real(n1**3)               
         
         LinVec_2_theta = matmul( map_mod1_mod2%dm%tv, LinVec_1 ) + matmul( map_mod1_mod2%dm%mi, LinVec_1_a ) &
                        + matmul( map_mod1_mod2%dm%fx_p, LinVec_1_b )! approximate delta theta^D         
                      
         j=1
         do i=1,Mesh1_O%nnodes
            
               ! Mesh1_O%Orientation = theta^S|_op + delta theta^S
            call SmllRotTrans( 'orientation', LinVec_1_theta(j  )+LinVec_1(j  ) &
                                            , LinVec_1_theta(j+1)+LinVec_1(j+1) &
                                            , LinVec_1_theta(j+2)+LinVec_1(j+2) & 
                                            , Mesh1_O%Orientation(:,:,i), ErrStat=ErrStat, ErrMsg=ErrMsg)            

                        
               ! Mesh1_O%TranslationVel = theta^S|_op + delta theta^S
            Mesh1_O%TranslationVel(:,i) = Mesh1_O_op%TranslationVel(:,i) + LinVec_1_a(j:j+2)
            
               ! Mesh1_O%RotationVel = theta^S|_op + delta theta^S
            Mesh1_O%RotationVel(:,i) = Mesh1_O_op%RotationVel(:,i) + LinVec_1_b(j:j+2)
            
            j = j+3            
            !IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))
                        
         end do
         

         !M( theta^S|_op + delta theta^S )
         CALL Transfer_Point_to_Point( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))     

         ! delta theta^D = M( theta^S|_op + delta theta^S ) - x^D|_op
         j=1
         do i=1,Mesh2_O%nnodes
            LinVec_2(j:j+2) = Mesh2_I%TranslationVel(:,i) - Mesh2_I_op%TranslationVel(:,i)          
            j = j+3
         end do
                  
         print *, 'Translational Velocity: ', n1, TwoNorm( LinVec_2 - LinVec_2_theta )
                              
      end do      
      
      ! ..................
      ! Rotational acceleration:
      ! ..................
      
      call meshcopy( Mesh1_O_op, Mesh1_O, MESH_UPDATECOPY, ErrStat, ErrMsg) 
      do n1=1,15
         
         ! perturb x^S:
         call random_number( LinVec_1 )
         LinVec_1 = LinVec_1/real(n1**3)               
         
         LinVec_2_theta = matmul( map_mod1_mod2%dm%mi, LinVec_1 )  ! approximate delta theta^D         
                      
         j=1
         do i=1,Mesh1_O%nnodes

               ! Mesh1_O%RotationAcc = theta^S|_op + delta theta^S
            Mesh1_O%RotationAcc(:,i) = Mesh1_O_op%RotationAcc(:,i) + LinVec_1(j:j+2)
                                                
            j = j+3            
                        
         end do
         

         !M( theta^S|_op + delta theta^S )
         CALL Transfer_Point_to_Point( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))     

         ! delta theta^D = M( theta^S|_op + delta theta^S ) - x^D|_op
         j=1
         do i=1,Mesh2_O%nnodes
            LinVec_2(j:j+2) = Mesh2_I%RotationAcc(:,i) - Mesh2_I_op%RotationAcc(:,i)          
            j = j+3
         end do
                  
         print *, 'Rotational RotationAcc: ', n1, TwoNorm( LinVec_2 - LinVec_2_theta )                              
      end do            
                  
      ! ..................
      ! Translational acceleration:
      ! ..................
      
      call meshcopy( Mesh1_O_op, Mesh1_O, MESH_UPDATECOPY, ErrStat, ErrMsg) 
      do n1=1,15
         
         ! perturb x^S:
         call random_number( LinVec_1 )
         LinVec_1 = LinVec_1/real(n1**3)               
         call random_number( LinVec_1_a )
         LinVec_1_a = LinVec_1_a/real(n1**3)               
         call random_number( LinVec_1_b )
         LinVec_1_b = LinVec_1_b/real(n1**3)               
         call random_number( LinVec_1_c )
         LinVec_1_c = LinVec_1_c/real(n1**3)               
         
         LinVec_2_theta = matmul( map_mod1_mod2%dm%ta1, LinVec_1   ) + matmul( map_mod1_mod2%dm%ta2,  LinVec_1_a ) &
                        + matmul( map_mod1_mod2%dm%mi,  LinVec_1_b ) + matmul( map_mod1_mod2%dm%fx_p, LinVec_1_c )! approximate delta theta^D         
                      
         j=1
         do i=1,Mesh1_O%nnodes
            
               ! Mesh1_O%Orientation = theta^S|_op + delta theta^S
            call SmllRotTrans( 'orientation', LinVec_1_theta(j  )+LinVec_1(j  ) &
                                            , LinVec_1_theta(j+1)+LinVec_1(j+1) &
                                            , LinVec_1_theta(j+2)+LinVec_1(j+2) & 
                                            , Mesh1_O%Orientation(:,:,i), ErrStat=ErrStat, ErrMsg=ErrMsg)            

                        
               ! Mesh1_O%RotationVel = theta^S|_op + delta theta^S
            Mesh1_O%RotationVel(:,i) = Mesh1_O_op%RotationVel(:,i) + LinVec_1_a(j:j+2)
            
               ! Mesh1_O%TranslationAcc = theta^S|_op + delta theta^S
            Mesh1_O%TranslationAcc(:,i) = Mesh1_O_op%TranslationAcc(:,i) + LinVec_1_b(j:j+2)
            
               ! Mesh1_O%RotationAcc = theta^S|_op + delta theta^S
            Mesh1_O%RotationAcc(:,i) = Mesh1_O_op%RotationAcc(:,i) + LinVec_1_c(j:j+2)
            
            j = j+3            
            !IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))
                        
         end do
         

         !M( theta^S|_op + delta theta^S )
         CALL Transfer_Point_to_Point( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))     

         ! delta theta^D = M( theta^S|_op + delta theta^S ) - x^D|_op
         j=1
         do i=1,Mesh2_O%nnodes
            LinVec_2(j:j+2) = Mesh2_I%TranslationAcc(:,i) - Mesh2_I_op%TranslationAcc(:,i)          
            j = j+3
         end do
                  
         print *, 'Translational Acceleration: ', n1, TwoNorm( LinVec_2 - LinVec_2_theta )
                              
      end do      
      
      
end if


else      
      
end if ! linearization or mapping test

      ! ..............................................................................................................................   
      ! Destroy data structures:
      ! ..............................................................................................................................   

      CALL MeshDestroy( mesh1_I, ErrStat, ErrMsg );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      CALL MeshDestroy( mesh1_O, ErrStat, ErrMsg );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      CALL MeshDestroy( mesh2_I, ErrStat, ErrMsg );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      CALL MeshDestroy( mesh2_O, ErrStat, ErrMsg );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

      call MeshMapDestroy(Map_Mod1_Mod2, ErrStat, ErrMsg);IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      call MeshMapDestroy(Map_Mod2_Mod1, ErrStat, ErrMsg);IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
                    
      if (allocated(LinVec_1)) deallocate(LinVec_1)
      if (allocated(LinVec_2)) deallocate(LinVec_2)      
      if (allocated(LinVec_1_theta)) deallocate(LinVec_1_theta)
      if (allocated(LinVec_2_theta)) deallocate(LinVec_2_theta)      
      if (allocated(LinVec_1_a)) deallocate(LinVec_1_a)
      if (allocated(LinVec_1_b)) deallocate(LinVec_1_b)
      if (allocated(LinVec_1_c)) deallocate(LinVec_1_c)      
      
   end do
   
   
   
   
contains   
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
      real(Dbki) :: DCM(3,3)
      
      
      Mesh1Type = ELEMENT_LINE2
      Mesh2Type = ELEMENT_LINE2
      
      !.........................
      ! Mesh1 (Output: Motions)
      !.........................
      
      Nnodes = 5            
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
         
         CALL MeshPositionNode ( mesh1_O, j, position, ErrStat, ErrMsg, orient=orientation )     
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
         DCM = GetDCM(Angle,2)
         
            ! place nodes in a line
         z = (j-1.0)/(NNodes-1.0)
         position = (/(sqrt(3.0)/3.0)*z-(sqrt(3.0)/6.0), 0.0_ReKi, z /)
         CALL MeshPositionNode ( mesh2_O, j, position, ErrStat, ErrMsg, Orient=DCM )     
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
   
      Mesh1Type = ELEMENT_Line2
      Mesh2Type = ELEMENT_Line2
                        
      
      !.........................
      ! Mesh1 (Output: Motions)
      !.........................
      
      Nnodes = 11
            
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
         x        = (j-1.0_ReKi) * dx
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
         
         x        = (j-1.0_ReKi) * dx 
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
                    
               
      !.........................
      ! Mesh2 (Output: Loads)
      !.........................
                    
      NNodes = 31 !19   
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
         x           = (j-1.0_ReKi) * dx
         yp          = a*omega*cos(omega*x)
         position    = (/x, a*sin(omega*x), 0.0_ReKi /)                    
         Angle       = atan(yp) 
         
         Orientation(:,1) = (/ COS(Angle), -1.*SIN(Angle), 0.0_ReKi /)
         Orientation(:,2) = (/ SIN(Angle),     COS(Angle), 0.0_ReKi /)
         Orientation(:,3) = (/      0.,        0.0,        1.0 /)
                 
         CALL MeshPositionNode ( mesh2_O, j, position, ErrStat, ErrMsg, &
               Orient= Orientation )    
         
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
   
END PROGRAM

