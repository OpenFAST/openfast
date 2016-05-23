   
subroutine Test_TestMeshMapping()

   USE NWTC_Library
   use TestMeshMapping_Mod
      
   IMPLICIT NONE   

   REAL(ReKi), allocatable :: LinVec_2_tmp(:), LinVec_2_a_tmp(:)
   
   REAL(ReKi), allocatable :: LinVec_1(:), LinVec_1_a(:), LinVec_1_b(:), LinVec_1_c(:)
   REAL(ReKi), allocatable :: LinVec_2(:), LinVec_2_a(:), LinVec_2_b(:), LinVec_2_c(:) ! data for linearization testing
   logical, parameter :: TestLinearization = .true.
   integer :: seed(15) 
   
   INTEGER :: TestNumber 
   CHARACTER(256) :: BinOutputName 
   
   seed = 0
   call random_seed(SIZE=n1)
   CALL RANDOM_SEED (PUT = SEED (1 : n1))
   
   ! CALL NWTC_Init(  )   
   
   DO TestNumber=1,13

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
      
      Mesh1_O%RemapFlag = .false.
      Mesh2_O%RemapFlag = .false.
      
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
      
      call AllocAry( LinVec_1,   Mesh1_O%nnodes*3, 'LinVec_1',   ErrStat, ErrMsg);       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      call AllocAry( LinVec_1_a, Mesh1_O%nnodes*3, 'LinVec_1_a', ErrStat, ErrMsg);       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      call AllocAry( LinVec_1_b, Mesh1_O%nnodes*3, 'LinVec_1_b', ErrStat, ErrMsg);       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      call AllocAry( LinVec_1_c, Mesh1_O%nnodes*3, 'LinVec_1_c', ErrStat, ErrMsg);       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      call AllocAry( LinVec_2,   Mesh2_O%nnodes*3, 'LinVec_2',   ErrStat, ErrMsg);       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      call AllocAry( LinVec_2_a, Mesh2_O%nnodes*3, 'LinVec_2_a', ErrStat, ErrMsg);       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      call AllocAry( LinVec_2_b, Mesh2_O%nnodes*3, 'LinVec_2_b', ErrStat, ErrMsg);       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      call AllocAry( LinVec_2_c, Mesh2_O%nnodes*3, 'LinVec_2_c', ErrStat, ErrMsg);       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      
      call AllocAry( LinVec_2_tmp,   Mesh2_O%nnodes*9, 'LinVec_2_tmp',   ErrStat, ErrMsg);       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      call AllocAry( LinVec_2_a_tmp, Mesh2_O%nnodes*9, 'LinVec_2_a_tmp', ErrStat, ErrMsg);       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

         ! make copies of operating point values
      CALL MeshCopy (mesh1_O, mesh1_O_op, MESH_NEWCOPY, ErrStat , ErrMsg ); IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      CALL MeshCopy (mesh2_O, mesh2_O_op, MESH_NEWCOPY, ErrStat , ErrMsg ); IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))            
      CALL MeshCopy (mesh1_I, mesh1_I_op, MESH_NEWCOPY, ErrStat , ErrMsg ); IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))            
      CALL MeshCopy (mesh2_I, mesh2_I_op, MESH_NEWCOPY, ErrStat , ErrMsg ); IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))            
            
write(99, *) '---------------------------------------------------------------'
write(99, *) '   Test ', TestNumber
write(99, *) '---------------------------------------------------------------'
      
         ! all of these meshes are set to their operating point values

      IF (Mesh1Type == ELEMENT_POINT ) THEN
         IF ( Mesh2Type == ELEMENT_POINT ) THEN 
            write(*,*) 'POINT-to-POINT'
            CALL Linearize_Point_to_Point( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))     
            CALL Linearize_Point_to_Point( Mesh2_O, Mesh1_I, Map_Mod2_Mod1, ErrStat, ErrMsg, Mesh2_I, Mesh1_O );   IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))                                                
         ELSEIF ( Mesh2Type == ELEMENT_LINE2) THEN                                                                
            write(*,*) 'POINT-to-LINE2 MOTIONS and LINE2-to-POINT LOADS' 
            CALL Linearize_Point_to_Line2( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))     
            CALL Linearize_Line2_to_Point( Mesh2_O, Mesh1_I, Map_Mod2_Mod1, ErrStat, ErrMsg, Mesh2_I, Mesh1_O );   IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))         
         END IF                                                                                                   
      ELSEIF ( Mesh1Type == ELEMENT_LINE2 ) THEN                                                                  
         IF ( Mesh2Type == ELEMENT_LINE2 ) THEN                                                                   
            write(*,*) 'LINE2-to-LINE2'
            CALL Linearize_Line2_to_Line2( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))   
            CALL Linearize_Line2_to_Line2( Mesh2_O, Mesh1_I, Map_Mod2_Mod1, ErrStat, ErrMsg, Mesh2_I, Mesh1_O );   IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))        
         ELSEIF ( Mesh2Type == ELEMENT_POINT ) THEN        
            write(*,*) 'LINE2-to-POINT MOTIONS and POINT-to-LINE2 LOADS' 
            CALL Linearize_Line2_to_Point( Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg );                     IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))                  
            CALL Linearize_Point_to_Line2( Mesh2_O, Mesh1_I, Map_Mod2_Mod1, ErrStat, ErrMsg, Mesh2_I, Mesh1_O );   IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))              
         END IF            
      END IF
      
             
      if (allocated(map_mod1_mod2%dm%mi))   call wrmatrix(map_mod1_mod2%dm%mi,  99,Fmt,'Mi')
      if (allocated(map_mod1_mod2%dm%fx_p)) call wrmatrix(map_mod1_mod2%dm%fx_p,99,Fmt,'fx_p')
      if (allocated(map_mod1_mod2%dm%tv))   call wrmatrix(map_mod1_mod2%dm%tv,  99,Fmt,'tv')
      if (allocated(map_mod1_mod2%dm%ta1))  call wrmatrix(map_mod1_mod2%dm%ta1, 99,Fmt,'ta1')
      if (allocated(map_mod1_mod2%dm%ta2))  call wrmatrix(map_mod1_mod2%dm%ta2, 99,Fmt,'ta2')
            
      if (allocated(map_mod2_mod1%dm%li ))  call wrmatrix(map_mod2_mod1%dm%li,  99,Fmt,'li')
      if (allocated(map_mod2_mod1%dm%M_u))  call wrmatrix(map_mod2_mod1%dm%M_u, 99,Fmt,'m_u')
      if (allocated(map_mod2_mod1%dm%M_t))  call wrmatrix(map_mod2_mod1%dm%M_t, 99,Fmt,'m_t')
      if (allocated(map_mod2_mod1%dm%M_f))  call wrmatrix(map_mod2_mod1%dm%M_f, 99,Fmt,'m_f')      
      
      
      write(*,ErrTxtFmt) 'Field','#','Absolute Err','Relative Err','Max perturb:D','Max OP:S','Max delta:S','Max delta:S','Max delta:S','Max delta:S'
      write(*,ErrTxtFmt) '-----','-','------------','------------','-------------','--------','-----------','-----------','-----------','-----------'
      ! ..................
      ! Rotational displacement:
      ! ..................
      
      !call meshcopy( Mesh1_O_op, Mesh1_O, MESH_UPDATECOPY, ErrStat, ErrMsg) 
      do n1=1,15
         
         ! perturb x^S:
         call getRotationPerturb(LinVec_1)
                  
         LinVec_2_a = matmul( map_mod1_mod2%dm%mi, LinVec_1 ) ! approximate delta theta^D         
             
         !call wrmatrix(LinVec_1,  99,Fmt,'LinVec_1')
         !call wrmatrix(LinVec_2_a,  99,Fmt,'LinVec_2_a')
         
         j=1
         do i=1,Mesh1_O%nnodes

               ! delta theta^S
            !Orientation = EulerConstruct(LinVec_1(j:j+2))
            call SmllRotTrans( 'orientation', LinVec_1(j  ) &
                                            , LinVec_1(j+1) &
                                            , LinVec_1(j+2) & 
                                            , Orientation, ErrStat=ErrStat, ErrMsg=ErrMsg)            
            
               ! Mesh1_O%Orientation = theta^S|_op + delta theta^S
            Mesh1_O%Orientation(:,:,i) = matmul( Mesh1_O_op%Orientation(:,:,i), Orientation )
                                    
            j = j+3            
            !IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))
         end do

         !M( x^S|_op + delta x^S )
         call TransferMotionData()

         ! delta x^D = M( x^S|_op + delta x^S ) - x^D|_op
         j=1
         do i=1,Mesh2_I%nnodes
            !theta^D|_op + delta theta^D = theta^D -> 
            !delta theta^D = (theta^D|_op)^(-1) * theta^D
            Orientation = matmul( transpose(Mesh2_I_op%Orientation(:,:,i)), Mesh2_I%Orientation(:,:,i) )
            
            !LinVec_2(j:j+2) = EulerExtract(Orientation)
            LinVec_2(j:j+2) = GetSmllRotAngs(Orientation,ErrStat,ErrMsg)    
            IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))
            
            j = j+3
         end do
         !call wrmatrix(LinVec_2-LinVec_2_a,  99,Fmt,'LinVec_2')
         
         e=TwoNorm( LinVec_2 - LinVec_2_a )               
         !if (.not. equalrealNos(e,0.0_ReKi)) &
            write(*,ErrFmt) 'Rotational Displacement: ', n1, e, e/max(mmin, TwoNorm(LinVec_2)), maxval(abs(LinVec_2)), 0.0, maxval(abs(LinVec_1)) 
            
      end do
      write(*,*) ! blank line
      
      ! ..................
      ! Translational displacement:
      ! ..................
      mm = maxval( abs(Mesh1_O_op%TranslationDisp ) )
      call meshcopy( Mesh1_O_op, Mesh1_O, MESH_UPDATECOPY, ErrStat, ErrMsg) 
      do n1=1,15
         
         ! perturb x^S:
         call getRandomVector(LinVec_1, mm)
         call getRotationPerturb(LinVec_1_a)
                  
                  
         LinVec_2_a = matmul( map_mod1_mod2%dm%mi, LinVec_1 ) + matmul( map_mod1_mod2%dm%fx_p, LinVec_1_a ) ! approximate delta theta^D         
                      
         j=1
         do i=1,Mesh1_O%nnodes

               ! Mesh1_O%TranslationDisp = u^S|_op + delta u^S
            Mesh1_O%TranslationDisp(:,i) = Mesh1_O_op%TranslationDisp(:,i) + LinVec_1(j:j+2)
            
               ! Mesh1_O%Orientation = theta^S|_op + delta theta^S                                    
            call SmllRotTrans( 'orientation', LinVec_1_a(j  ) &
                                            , LinVec_1_a(j+1) &
                                            , LinVec_1_a(j+2) & 
                                            , Orientation, ErrStat=ErrStat, ErrMsg=ErrMsg)            
            Mesh1_O%Orientation(:,:,i) = matmul( Mesh1_O_op%Orientation(:,:,i), Orientation )
            
            
            j = j+3            
            !IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))
                        
         end do
         

         !M( x^S|_op + delta x^S )
         call TransferMotionData()

         ! delta x^D = M( x^S|_op + delta x^S ) - x^D|_op
         j=1
         do i=1,Mesh2_I%nnodes
            LinVec_2(j:j+2) = Mesh2_I%TranslationDisp(:,i) - Mesh2_I_op%TranslationDisp(:,i)          
            j = j+3
         end do
         
         e=TwoNorm( LinVec_2 - LinVec_2_a )               
         !if (.not. equalrealNos(e,0.0_ReKi)) &
            write(*,ErrFmt) 'Translational Displacement: ', n1, e, e/max(mmin, TwoNorm(LinVec_2)),  maxval(abs(LinVec_2)), &
               mm , maxval(abs(LinVec_1)), maxval(abs(LinVec_1_a))  
                              
      end do      
      
      write(*,*) ! blank line
      
      ! ..................
      ! Rotational velocity:
      ! ..................
      mm = maxval( abs(Mesh1_O_op%RotationVel ) )
      call meshcopy( Mesh1_O_op, Mesh1_O, MESH_UPDATECOPY, ErrStat, ErrMsg) 
      do n1=1,15
         
         ! perturb x^S:
         call getRandomVector(LinVec_1, mm)
         
         LinVec_2_a = matmul( map_mod1_mod2%dm%mi, LinVec_1 )  ! approximate delta theta^D         
                      
         j=1
         do i=1,Mesh1_O%nnodes

               ! Mesh1_O%RotationVel = omega^S|_op + delta omega^S
            Mesh1_O%RotationVel(:,i) = Mesh1_O_op%RotationVel(:,i) + LinVec_1(j:j+2)
                                                
            j = j+3            
                        
         end do
         
         !M( x^S|_op + delta x^S )
         call TransferMotionData()

         ! delta x^D = M( x^S|_op + delta x^S ) - x^D|_op
         j=1
         do i=1,Mesh2_I%nnodes
            LinVec_2(j:j+2) = Mesh2_I%RotationVel(:,i) - Mesh2_I_op%RotationVel(:,i)          
            j = j+3
         end do
                  
         e=TwoNorm( LinVec_2 - LinVec_2_a )               
         !if (.not. equalrealNos(e,0.0_ReKi)) &
            write(*,ErrFmt) 'Rotational Velocity: ', n1, e, e/max(mmin, TwoNorm(LinVec_2)), maxval(abs(LinVec_2)), &
               mm , maxval(abs(LinVec_1)) 
      end do            
      
      write(*,*) ! blank line
      
      ! ..................
      ! Translational velocity:
      ! ..................
      mm = maxval( abs(Mesh1_O_op%TranslationVel ) )
      call meshcopy( Mesh1_O_op, Mesh1_O, MESH_UPDATECOPY, ErrStat, ErrMsg) 
      do n1=1,15
         
         ! perturb x^S:
         call getRotationPerturb(LinVec_1)
         call getRandomVector(LinVec_1_a)
         call getRandomVector(LinVec_1_b)
         
         LinVec_2_a = matmul( map_mod1_mod2%dm%tv, LinVec_1 ) + matmul( map_mod1_mod2%dm%mi, LinVec_1_a ) &
                    + matmul( map_mod1_mod2%dm%fx_p, LinVec_1_b )! approximate delta theta^D         
                      
         j=1
         do i=1,Mesh1_O%nnodes
            
               ! Mesh1_O%Orientation = theta^S|_op + delta theta^S
            call SmllRotTrans( 'orientation', LinVec_1(j  ) &
                                            , LinVec_1(j+1) &
                                            , LinVec_1(j+2) & 
                                            , Orientation, ErrStat=ErrStat, ErrMsg=ErrMsg)            
            Mesh1_O%Orientation(:,:,i) = matmul( Mesh1_O_op%Orientation(:,:,i), Orientation )

                        
               ! Mesh1_O%TranslationVel = v^S|_op + delta v^S
            Mesh1_O%TranslationVel(:,i) = Mesh1_O_op%TranslationVel(:,i) + LinVec_1_a(j:j+2)
            
               ! Mesh1_O%RotationVel = omega^S|_op + delta omega^S
            Mesh1_O%RotationVel(:,i) = Mesh1_O_op%RotationVel(:,i) + LinVec_1_b(j:j+2)
            
            j = j+3            
            !IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))
                        
         end do
         

         !M( x^S|_op + delta x^S )
         call TransferMotionData()

         ! delta x^D = M( x^S|_op + delta x^S ) - x^D|_op
         j=1
         do i=1,Mesh2_I%nnodes
            LinVec_2(j:j+2) = Mesh2_I%TranslationVel(:,i) - Mesh2_I_op%TranslationVel(:,i)          
            j = j+3
         end do
                  
         e=TwoNorm( LinVec_2 - LinVec_2_a )               
         !if (.not. equalrealNos(e,0.0_ReKi)) &
            write(*,ErrFmt) 'Translational Velocity: ', n1, e, e/max(mmin, TwoNorm(LinVec_2)), maxval(abs(LinVec_2)), &
               mm, maxval(abs(LinVec_1)), maxval(abs(LinVec_1_a)), maxval(abs(LinVec_1_b)) 
                              
      end do      
      write(*,*) ! blank line
      
      ! ..................
      ! Rotational acceleration:
      ! ..................
      mm = maxval( abs(Mesh1_O_op%RotationAcc ) )
            
      call meshcopy( Mesh1_O_op, Mesh1_O, MESH_UPDATECOPY, ErrStat, ErrMsg) 
      do n1=1,15
         
         ! perturb x^S:
         call getRandomVector(LinVec_1)
         
         LinVec_2_a = matmul( map_mod1_mod2%dm%mi, LinVec_1 )  ! approximate delta theta^D         
                      
         j=1
         do i=1,Mesh1_O%nnodes

               ! Mesh1_O%RotationAcc = alpha^S|_op + delta alpha^S
            Mesh1_O%RotationAcc(:,i) = Mesh1_O_op%RotationAcc(:,i) + LinVec_1(j:j+2)
                                                
            j = j+3            
                        
         end do
         

         !M( x^S|_op + delta x^S )
         call TransferMotionData()

         ! delta x^D = M( x^S|_op + delta x^S ) - x^D|_op
         j=1
         do i=1,Mesh2_I%nnodes
            LinVec_2(j:j+2) = Mesh2_I%RotationAcc(:,i) - Mesh2_I_op%RotationAcc(:,i)          
            j = j+3
         end do
                  
         e=TwoNorm( LinVec_2 - LinVec_2_a )               
         !if (.not. equalrealNos(e,0.0_ReKi)) &
            write(*,ErrFmt) 'Rotational Acceleration: ', n1, e, e/max(mmin, TwoNorm(LinVec_2)), maxval(abs(LinVec_2)), &
              mm , maxval(abs(LinVec_1)) 
      end do            
      write(*,*) ! blank line
      
      ! ..................
      ! Translational acceleration:
      ! ..................

      call meshcopy( Mesh1_O_op, Mesh1_O, MESH_UPDATECOPY, ErrStat, ErrMsg) 
      do n1=1,15
         
         ! perturb x^S:
         call getRotationPerturb(LinVec_1) ! delta theta
         call getRandomVector(LinVec_1_a)  ! delta omega 
         call getRandomVector(LinVec_1_b)  ! delta a
         call getRandomVector(LinVec_1_c)  ! delta alpha
         
         LinVec_2_a = matmul( map_mod1_mod2%dm%ta1, LinVec_1   ) + matmul( map_mod1_mod2%dm%ta2,  LinVec_1_a ) &
                    + matmul( map_mod1_mod2%dm%mi,  LinVec_1_b ) + matmul( map_mod1_mod2%dm%fx_p, LinVec_1_c )! approximate delta theta^D         
             
         
         j=1
         do i=1,Mesh1_O%nnodes
            
               ! Mesh1_O%Orientation = theta^S|_op + delta theta^S
            call SmllRotTrans( 'orientation', LinVec_1(j  ) &
                                            , LinVec_1(j+1) &
                                            , LinVec_1(j+2) & 
                                            , Orientation, ErrStat=ErrStat, ErrMsg=ErrMsg)            
            Mesh1_O%Orientation(:,:,i) = matmul( Mesh1_O_op%Orientation(:,:,i), Orientation )

                        
               ! Mesh1_O%RotationVel = omega^S|_op + delta omega^S
            Mesh1_O%RotationVel(:,i) = Mesh1_O_op%RotationVel(:,i) + LinVec_1_a(j:j+2)
            
               ! Mesh1_O%TranslationAcc = a^S|_op + delta a^S
            Mesh1_O%TranslationAcc(:,i) = Mesh1_O_op%TranslationAcc(:,i) + LinVec_1_b(j:j+2)
            
               ! Mesh1_O%RotationAcc = alpha^S|_op + delta alpha^S
            Mesh1_O%RotationAcc(:,i) = Mesh1_O_op%RotationAcc(:,i) + LinVec_1_c(j:j+2)
            
            j = j+3            
            !IF (ErrStat /= ErrID_None) CALL WrScr("*******"//TRIM(ErrMsg))
                        
         end do
         

         !M( x^S|_op + delta x^S )
         call TransferMotionData()

         ! delta x^D = M( x^S|_op + delta x^S ) - x^D|_op
         j=1
         do i=1,Mesh2_I%nnodes
            LinVec_2(j:j+2) = Mesh2_I%TranslationAcc(:,i) - Mesh2_I_op%TranslationAcc(:,i)          
            j = j+3
         end do
                  
         e=TwoNorm( LinVec_2 - LinVec_2_a )               
         !if (.not. equalrealNos(e,0.0_ReKi)) &
            write(*,ErrFmt) 'Translational Acceleration: ', n1, e, e/max(mmin, TwoNorm(LinVec_2)), maxval(abs(LinVec_2)), &
               mm, maxval(abs(LinVec_1)), maxval(abs(LinVec_1_a)), maxval(abs(LinVec_1_b)), maxval(abs(LinVec_1_c)) 
                              
      end do      
      write(*,*) ! blank line
      
      
      ! ..................
      ! Forces:
      ! ..................

            ! get max Force :
      mf=maxval( abs(Mesh2_O_op%Force ) )
      !mf = sqrt( EPSILON(mf)) * max(1.0_ReKi, mf)
      
      
      call meshcopy( Mesh2_O_op, Mesh2_O, MESH_UPDATECOPY, ErrStat, ErrMsg) 
      do n1=1,15
         
         ! perturb x^S:
         call getRandomVector(LinVec_2, mf)
         
         LinVec_1_a = matmul( map_mod2_mod1%dm%li, LinVec_2 )  ! approximate delta theta^D         
                      
         j=1
         do i=1,Mesh2_O%nnodes

               ! Mesh2_O%Force = theta^S|_op + delta theta^S
            Mesh2_O%Force(:,i) = Mesh2_O_op%Force(:,i) + LinVec_2(j:j+2)
                                                
            j = j+3            
                        
         end do
         

         !M( x^S|_op + delta x^S )
         call TransferLoadData()
         
         ! delta x^D = M( x^S|_op + delta x^S ) - x^D|_op
         j=1
         do i=1,Mesh1_I%nnodes
            LinVec_1(j:j+2) = Mesh1_I%Force(:,i) - Mesh1_I_op%Force(:,i)   
            j = j+3
         end do
                  
         e=TwoNorm( LinVec_1 - LinVec_1_a )               
         !if (.not. equalrealNos(mf+e,mf)) &
            write(*,ErrFmt) 'Forces: ', n1, e, e/max(mmin, maxval(abs(LinVec_1))), maxval(abs(LinVec_1)), mf, maxval(abs(LinVec_2))         
      end do            
      write(*,*) ! blank line
      
      ! ..................
      ! Moments:
      ! ..................
      
            ! get max Moment :
      mf=maxval( abs(Mesh2_O_op%Force ) )     
      mm=maxval( abs(Mesh2_O_op%Moment ) )     
      
      !call meshcopy( Mesh2_I_op, Mesh2_I, MESH_UPDATECOPY, ErrStat, ErrMsg) 
      !call meshcopy( Mesh1_I_op, Mesh1_I, MESH_UPDATECOPY, ErrStat, ErrMsg) 
      
      call meshcopy( Mesh2_O_op, Mesh2_O, MESH_UPDATECOPY, ErrStat, ErrMsg) 
      call meshcopy( Mesh1_O_op, Mesh1_O, MESH_UPDATECOPY, ErrStat, ErrMsg) 
      
      do n1=1,15

         
         ! perturb x^S:
         !......
         call getRandomVector(LinVec_2_b, max(1.0_ReKi,mf)) ! delta f
         call getRandomVector(LinVec_2_c, max(1.0_ReKi,mm)) ! delta m 

            ! for u and theta, we actually perturb x^D and transfer to x^S to get delta x^S
         call getRandomVector(LinVec_1_b) !delta u
         call getRotationPerturb(LinVec_1_c) ! delta theta
         
         
         ! get Mesh1_O (destination, perturbed)
         
         j=1
         do i=1,Mesh1_O%nnodes
            
               ! bjj: we need orientation and translation to be consistent on same meshes as well as on transfered meshes
            
               ! Mesh2_I%TranslationDisp = u^S|_op + delta u^S
            if (allocated(map_mod2_mod1%dm%m_u)) Mesh1_O%TranslationDisp(:,i) = Mesh1_O_op%TranslationDisp(:,i) + LinVec_1_b(j:j+2)
            
            
               ! Mesh2_I%Orientation = theta^S|_op + delta theta^S
            call SmllRotTrans( 'orientation', LinVec_1_c(j  ) &
                                            , LinVec_1_c(j+1) &
                                            , LinVec_1_c(j+2) & 
                                            , Orientation, ErrStat=ErrStat, ErrMsg=ErrMsg)            
            Mesh1_O%Orientation(:,:,i) = matmul( Mesh1_O_op%Orientation(:,:,i), Orientation )
         end do
         
         call TransferMotionData() ! this sets Mesh2_I, the displaced source mesh
                  
            ! now figure out what that means for the x^S perturbations:
         
            ! LinVec_2 for delta u^S
         if (allocated(map_mod2_mod1%dm%m_u)) then
            j=1
            do i=1,Mesh2_I%nnodes
               LinVec_2(j:j+2) = Mesh2_I%TranslationDisp(:,i) - Mesh2_I_op%TranslationDisp(:,i)
               j = j+3
            end do
         else
            LinVec_2 = 0
         end if
         
            ! LinVec_2_a for delta theta^S
         if (allocated(map_mod2_mod1%dm%m_t)) then
            j=1
            do i=1,Mesh2_I%nnodes
               Orientation = matmul( transpose( Mesh2_I_op%Orientation(:,:,i) ), Mesh2_I%Orientation(:,:,i) )
               LinVec_2_a(j:j+2) = GetSmllRotAngs(Orientation,ErrStat,ErrMsg)
               j = j+3
            end do                  
         end if
         !......
         
         LinVec_1_a =   matmul( map_mod2_mod1%dm%m_f, LinVec_2_b ) & ! approximate delta theta^D (after we add the u and t components on the next lines)
                      + matmul( map_mod2_mod1%dm%li,  LinVec_2_c )  
                    
         if (allocated(map_mod2_mod1%dm%m_u)) LinVec_1_a = LinVec_1_a + matmul( map_mod2_mod1%dm%m_u, LinVec_2   ) ! u
         if (allocated(map_mod2_mod1%dm%m_t)) LinVec_1_a = LinVec_1_a + matmul( map_mod2_mod1%dm%m_t, LinVec_2_a ) ! theta
                                   
         j=1
         do i=1,Mesh2_O%nnodes
            
               ! bjj: we need orientation and translation to be consistent on same meshes as well as on transfered meshes;
               !      so we need to actually perturb the destination and transfer to the source to get the source perturbations
            
               ! Mesh2_O%Force = Force^S|_op + delta Force^S
            Mesh2_O%Force(:,i) = Mesh2_O_op%Force(:,i) + LinVec_2_b(j:j+2)

               ! Mesh2_O%Moment = Moment^S|_op + delta Moment^S
            Mesh2_O%Moment(:,i) = Mesh2_O_op%Moment(:,i) + LinVec_2_c(j:j+2)
            
            j = j+3            
                        
         end do
         

                           
         !M( x^S|_op + delta x^S )
         call TransferLoadData()
                  
         ! delta x^D = M( x^S|_op + delta x^S ) - x^D|_op
         j=1
         do i=1,Mesh1_I%nnodes
            LinVec_1(j:j+2) = Mesh1_I%Moment(:,i) - Mesh1_I_op%Moment(:,i)   
            j = j+3
         end do
                  
         e=TwoNorm( LinVec_1 - LinVec_1_a )               
         !if (.not. equalrealNos(e,0.0_ReKi)) &
            write(*,ErrFmt) 'Moments: ', n1, e, e/max(mmin, maxval(abs(LinVec_1))), maxval(abs(LinVec_1)), mm, &
                 maxval(abs(LinVec_2)), maxval(abs(LinVec_2_a)), maxval(abs(LinVec_2_b)), maxval(abs(LinVec_2_c))
      end do   
      
end if ! linearization or mapping test

      ! ..............................................................................................................................   
      ! Destroy data structures:
      ! ..............................................................................................................................   

      CALL MeshDestroy( mesh1_I, ErrStat, ErrMsg );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      CALL MeshDestroy( mesh1_O, ErrStat, ErrMsg );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      CALL MeshDestroy( mesh2_I, ErrStat, ErrMsg );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      CALL MeshDestroy( mesh2_O, ErrStat, ErrMsg );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

      CALL MeshDestroy( mesh1_I_op, ErrStat, ErrMsg );    IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      CALL MeshDestroy( mesh1_O_op, ErrStat, ErrMsg );    IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      CALL MeshDestroy( mesh2_I_op, ErrStat, ErrMsg );    IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      CALL MeshDestroy( mesh2_O_op, ErrStat, ErrMsg );    IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      
      
      call MeshMapDestroy(Map_Mod1_Mod2, ErrStat, ErrMsg);IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      call MeshMapDestroy(Map_Mod2_Mod1, ErrStat, ErrMsg);IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
                    
      if (allocated(LinVec_1)) deallocate(LinVec_1)
      if (allocated(LinVec_1_a)) deallocate(LinVec_1_a)
      if (allocated(LinVec_1_b)) deallocate(LinVec_1_b)
      if (allocated(LinVec_1_c)) deallocate(LinVec_1_c)      
      if (allocated(LinVec_2)) deallocate(LinVec_2)      
      if (allocated(LinVec_2_a)) deallocate(LinVec_2_a)      
      if (allocated(LinVec_2_b)) deallocate(LinVec_2_b)      
      if (allocated(LinVec_2_c)) deallocate(LinVec_2_c)      
      
      if (allocated(LinVec_2_tmp)) deallocate(LinVec_2_tmp)      
      if (allocated(LinVec_2_a_tmp)) deallocate(LinVec_2_a_tmp)      
      
   end do
              
   
   
END subroutine Test_TestMeshMapping

