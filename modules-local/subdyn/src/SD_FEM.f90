!**********************************************************************************************************************************
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
MODULE SD_FEM
  USE NWTC_Library
  USE SubDyn_Types
  
  
  
    CONTAINS
    
SUBROUTINE NodeCon(Init,p, ErrStat, ErrMsg)
  
!This Subroutine maps nodes to elements
! allocate for NodesConnE and NodesConnN                                                                               
  USE NWTC_Library
  USE SubDyn_Types
  USE qsort_c_module
  IMPLICIT NONE

  TYPE(SD_InitInputType),         INTENT( INOUT )  ::Init   
  TYPE(SD_ParameterType),         INTENT( IN    )  ::p  
  INTEGER(IntKi),                 INTENT(   OUT )  :: ErrStat     ! Error status of the operation
  CHARACTER(1024),                INTENT(   OUT )  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   
  !local variable
  INTEGER(IntKi) :: SortA(Init%MaxMemjnt,2)  !To sort nodes and elements
  INTEGER(IntKi) :: I,J,K  !counter
  
  ALLOCATE(Init%NodesConnE(Init%NNode, Init%MaxMemJnt+1), STAT=ErrStat)                !the row index is the number of the real node, i.e. ID, 1st col has number of elements attached to node, and 2nd col has element numbers (up to 10)                                    
  IF ( ErrStat /= 0 )  THEN                                                                                                
      ErrMsg = ' Error allocating NodesConnE matrix'                                                    
      ErrStat = ErrID_Fatal                                                                                                         
      RETURN                                                                                                              
  ENDIF                                                                                                                  
  Init%NodesConnE = 0                                                                                                    
                                                                                                                          
  ALLOCATE(Init%NodesConnN(Init%NNode, Init%MaxMemJnt+2), STAT=ErrStat)                                                    
  IF ( ErrStat /= 0 )  THEN                                                                                                
      ErrMsg = ' Error allocating NodesConnN matrix'                                                    
      ErrStat = ErrID_Fatal                                                                                                         
      RETURN                                                                                                              
  ENDIF                                                                                                                  
  Init%NodesConnN = 0                                                                                                    
                                                                                                                          
!   ! find the node connectivity, nodes/elements that connect to a common node                                             
!   ! initialize the temp array for sorting                                                                             
   SortA(Init%MaxMemjnt,2) = 0                                                                                            
                                                                                                                          
   DO I = 1, Init%NNode                                                                                                   
      Init%NodesConnN(I, 1) = NINT( Init%Nodes(I, 1) )      !This should not be needed, could remove the extra 1st column like for the other array                                                                      
                                                                                                                          
      k = 0                                                                                                               
      DO J = 1, Init%NElem                          !This should be vectorized                                                                      
         IF ( ( Init%Nodes(I, 1)==p%Elems(J, 2)) .OR. (Init%Nodes(I, 1)==p%Elems(J, 3) ) ) THEN   !If i-th nodeID matches 1st node or 2nd of j-th element                                                                   
            k = k + 1                                                                                                     
            Init%NodesConnE(I, k + 1) = p%Elems(J, 1)                                                                  
            Init%NodesConnN(I, k + 1) = p%Elems(J, 3)                                                                  
            IF ( Init%Nodes(I, 1)==p%Elems(J, 3) ) Init%NodesConnN(I, k + 1) = p%Elems(J, 2)     !If nodeID matches 2nd node of element                                                                
         ENDIF                                                                                                            
      ENDDO                                                                                                               
                                                                                                                          
      IF( k>1 )THEN ! sort the nodes ascendingly                                                                          
         SortA(1:k, 1) = Init%NodesConnN(I, 3:(k+2))                                                                      
         CALL QsortC( SortA(1:k, 1:2) )                                                                                   
         Init%NodesConnN(I, 3:(k+2)) = SortA(1:k, 1)                                                                      
      ENDIF                                                                                                               
                                                                                                                          
      Init%NodesConnE(I, 1) = k    !Store how many elements connect i-th node in 2nd column                                                                                       
      Init%NodesConnN(I, 2) = k                                                                                           
   ENDDO                            

END SUBROUTINE NodeCon

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE SD_Discrt(Init,p, ErrStat, ErrMsg)
   USE NWTC_Library
   USE SubDyn_Types
   IMPLICIT NONE

   TYPE(SD_InitInputType),       INTENT(INOUT)  ::Init
   TYPE(SD_ParameterType),       INTENT(INOUT)  ::p
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(1024),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
      ! local variable
   INTEGER                       :: I, J, Node1, Node2, Prop1, Prop2, flg, flg1, flg2
   INTEGER                       :: OldJointIndex(Init%NJoints)
   INTEGER                       :: NPE      ! node per element
   INTEGER(4)                    :: Sttus
   INTEGER                       :: TempNProp
   REAL(ReKi), ALLOCATABLE       :: TempProps(:, :)
   INTEGER, ALLOCATABLE          :: TempMembers(:, :) ,TempReacts(:,:)         
   CHARACTER(1024)               :: OutFile
   CHARACTER(  50)               :: tempStr ! string number of nodes in member
   CHARACTER(1024)               :: OutFmt
   INTEGER                       :: knode, kelem, kprop, nprop
   REAL(ReKi)                    :: x1, y1, z1, x2, y2, z2, dx, dy, dz, dd, dt, d1, d2, t1, t2
   
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   ! number of nodes per element
   IF( ( Init%FEMMod .GE. 0 ) .and. (Init%FEMMod .LE. 3) ) THEN
      NPE = 2 
   ENDIF
   
   ! Calculate total number of nodes according to divisions
   Init%NNode = Init%NJoints + ( Init%NDiv - 1 )*p%NMembers
   ! Total number of element
   Init%NElem = p%NMembers*Init%NDiv
   ! Total number of property sets (temp)
   TempNProp = Init%NElem*NPE
   
   ! Calculate total number of nodes and elements according to element types
   ! for 3-node or 4-node beam elements
   Init%NNode = Init%NNode + (NPE - 2)*Init%NElem
   Init%MembersCol = Init%MembersCol + (NPE - 2) 
   
   ! check the number of interior modes
   IF ( p%Nmodes .GT. 6*(Init%NNode - Init%NInterf - p%NReact) ) THEN
      WRITE(tempStr, *) 6*(Init%NNode - Init%NInterf - p%NReact)
      ErrMsg = ' The NModes must be less than or equal to'//TRIM(tempStr)
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF
   
   
   
   ! Allocate for nodes and elements and membernodes
   ALLOCATE(Init%Nodes(Init%NNode, Init%JointsCol), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      ErrMsg = ' Error allocating Nodes arrays'
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF
   Init%Nodes = 0
   
   ALLOCATE(p%Elems(Init%NElem, Init%MembersCol), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      ErrMsg = ' Error allocating Elems arrays'
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF
   p%Elems = 0

   ! for two-node element only, otherwise the number of nodes in one element is different
   ALLOCATE(Init%MemberNodes(p%NMembers, Init%NDiv+1), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      ErrMsg = ' Error allocating MemberNodes arrays'
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF
   Init%MemberNodes = 0

      ! Allocate Temp members and property sets
   ALLOCATE(TempMembers(p%NMembers, Init%MembersCol), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      ErrMsg = ' Error allocating TempProps arrays'
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF
   TempMembers = 0

   ALLOCATE(TempProps(TempNProp, Init%PropSetsCol), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      ErrMsg = ' Error allocating TempProps arrays'
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF
   TempProps = 0

   ! Initialize Nodes
   DO I = 1,Init%NJoints
      OldJointIndex(I) = Init%Joints(I, 1)
      Init%Nodes(I, 1) = I
      Init%Nodes(I, 2) = Init%Joints(I, 2)
      Init%Nodes(I, 3) = Init%Joints(I, 3)
      Init%Nodes(I, 4) = Init%Joints(I, 4)
   ENDDO
   
   ! Initialize TempMembers and Elems
   DO I = 1, p%NMembers
      TempMembers(I, 1) = I
      TempMembers(I, 4) = Init%Members(I, 4)
      TempMembers(I, 5) = Init%Members(I, 5)
      
      p%Elems(I,     1) = I
      p%Elems(I, NPE+2) = Init%Members(I, 4)
      p%Elems(I, NPE+3) = Init%Members(I, 5)
      
      Node1 = Init%Members(I, 2)
      Node2 = Init%Members(I, 3)
      
      flg1 = 0
      flg2 = 0
      DO J = 1, Init%NJoints
         IF ( Node1 == Init%Joints(J, 1) ) THEN
            TempMembers(I, 2) = J
            p%Elems(I, 2) = J
            flg1 = 1
         ENDIF
         IF ( Node2 == Init%Joints(J, 1) ) THEN
            TempMembers(I, 3) = J
            p%Elems(I, NPE+1) = J
            flg2 = 1
         ENDIF
         
      ENDDO
      
      IF ( (flg1 == 0) .OR. (flg2 == 0) ) THEN
         ErrMsg = ' Member has node not in the node list !'
         ErrStat = ErrID_Fatal
         RETURN
      ENDIF

   ENDDO
   
   ! Initialize Temp property set
   TempProps(1:Init%NPropSets, :) = Init%PropSets(1:Init%NPropSets, :)   
   
   ! Initialize boundary constraint vector
   ! Change the node number
   ALLOCATE(Init%BCs(6*p%NReact, 2), STAT=Sttus)  !!!! RRD: THIS MAY NEED TO CHANGE IF NOT ALL NODES ARE RESTRAINED
   IF ( Sttus /= 0 )  THEN
      ErrMsg = ' Error allocating BCs arrays'
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF
   Init%BCs = 0
      
    !Allocate array that will be p%Reacts renumbered and ordered so that ID does not play a role, just ordinal position number will count -RRD
   ALLOCATE(TempReacts(p%NReact, Init%ReactCol), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      ErrMsg = ' Error allocating TempReacts arrays'
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF
   TempReacts=0 !INitialize -RRD

   DO I = 1, p%NReact
      Node1 = p%Reacts(I, 1);  !NODE ID
      TempReacts(I,2:Init%ReactCol)=p%Reacts(I, 2:Init%ReactCol)  !Assign all the appropriate fixity to the new Reacts array -RRD
      
      flg = 0
      DO J = 1, Init%NJoints
         IF ( Node1 == Init%Joints(J, 1) ) THEN
            Node2 = J
            flg = 1
            TempReacts(I,1)=Node2      !New node ID for p!React  -RRD
            EXIT  !Exit J loop if node found -RRD
         ENDIF
      ENDDO
      
      IF (flg == 0) THEN
         ErrMsg = ' React has node not in the node list !'
         ErrStat = ErrID_Fatal
         RETURN
      ENDIF
      
      
      DO J = 1, 6
         Init%BCs( (I-1)*6+J, 1) = (Node2-1)*6+J;
         Init%BCs( (I-1)*6+J, 2) = p%Reacts(I, J+1);
      ENDDO
      
   ENDDO
   p%Reacts=TempReacts   !UPDATED REACTS
   
   ! Initialize interface constraint vector
   ! Change the node number
   ALLOCATE(Init%IntFc(6*Init%NInterf, 2), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      ErrMsg = ' Error allocating IntFc arrays'
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF
   Init%IntFc = 0
      
   DO I = 1, Init%NInterf
      Node1 = Init%Interf(I, 1);
      flg = 0
      DO J = 1, Init%NJoints
         IF ( Node1 == Init%Joints(J, 1) ) THEN
            Node2 = J
            flg = 1
         ENDIF
      ENDDO
      
      IF (flg == 0) THEN
         ErrMsg = ' Interf has node not in the node list !'
         ErrStat = ErrID_Fatal
         RETURN
      ENDIF
      
      DO J = 1, 6
         Init%IntFc( (I-1)*6+J, 1) = (Node2-1)*6+J;
         Init%IntFc( (I-1)*6+J, 2) = Init%Interf(I, J+1);
      ENDDO
   ENDDO
  
   ! Change numbering in concentrated mass matrix
   DO I = 1, Init%NCMass
      Node1 = NINT( Init%CMass(I, 1) )
      DO J = 1, Init%NJoints
         IF ( Node1 == Init%Joints(J, 1) ) THEN
            Init%CMass(I, 1) = J
         ENDIF
      ENDDO
   ENDDO

! discretize structure according to NDiv 

knode = Init%NJoints
kelem = 0
kprop = Init%NPropSets

IF (Init%NDiv .GT. 1) THEN
   DO I = 1, p%NMembers
      ! create new node
      Node1 = TempMembers(I, 2)
      Node2 = TempMembers(I, 3)
      
      IF ( Node1==Node2 ) THEN
         ErrMsg = ' Same starting and ending node in the member.'
         ErrStat = 4
         RETURN
      ENDIF
    
      
      
      Prop1 = TempMembers(I, 4)
      Prop2 = TempMembers(I, 5)
      
      Init%MemberNodes(I,           1) = Node1
      Init%MemberNodes(I, Init%NDiv+1) = Node2
      
      IF  ( ( .not. EqualRealNos(TempProps(Prop1, 2),TempProps(Prop2, 2) ) ) &
       .OR. ( .not. EqualRealNos(TempProps(Prop1, 3),TempProps(Prop2, 3) ) ) &
       .OR. ( .not. EqualRealNos(TempProps(Prop1, 4),TempProps(Prop2, 4) ) ) )  THEN
      
         ErrMsg = ' Material E,G and rho in a member must be the same'
         ErrStat = 3
         RETURN
      ENDIF

      x1 = Init%Nodes(Node1, 2)
      y1 = Init%Nodes(Node1, 3)
      z1 = Init%Nodes(Node1, 4)

      x2 = Init%Nodes(Node2, 2)
      y2 = Init%Nodes(Node2, 3)
      z2 = Init%Nodes(Node2, 4)
      
      dx = ( x2 - x1 )/Init%NDiv
      dy = ( y2 - y1 )/Init%NDiv
      dz = ( z2 - z1 )/Init%NDiv
      
      d1 = TempProps(Prop1, 5)
      t1 = TempProps(Prop1, 6)

      d2 = TempProps(Prop2, 5)
      t2 = TempProps(Prop2, 6)
      
      dd = ( d2 - d1 )/Init%NDiv
      dt = ( t2 - t1 )/Init%NDiv
      
      
      ! node connect to Node1
      knode = knode + 1
      Init%MemberNodes(I, 2) = knode
      CALL GetNewNode(knode, x1+dx, y1+dy, z1+dz, Init)
      
      IF ( ( .NOT.(EqualRealNos( dd , 0.0_ReKi ) ) ) .OR. &
           ( .NOT.( EqualRealNos( dt , 0.0_ReKi ) ) ) ) THEN   
           ! create a new property set 
           ! k, E, G, rho, d, t, Init
           
           kprop = kprop + 1
           CALL GetNewProp(kprop, TempProps(Prop1, 2), TempProps(Prop1, 3),&
                           TempProps(Prop1, 4), d1+dd, t1+dt, TempProps, TempNProp, Init%PropSetsCol)           
           kelem = kelem + 1
           CALL GetNewElem(kelem, Node1, knode, Prop1, kprop, p)  
           nprop = kprop              
      ELSE
           kelem = kelem + 1
           CALL GetNewElem(kelem, Node1, knode, Prop1, Prop1, p)                
           nprop = Prop1 
      ENDIF
      
      ! interior nodes
      
      DO J = 2, (Init%NDiv-1)
         knode = knode + 1
         Init%MemberNodes(I, J+1) = knode

         CALL GetNewNode(knode, x1 + J*dx, y1 + J*dy, z1 + J*dz, Init)
         
         IF ( ( .NOT.(EqualRealNos( dd , 0.0_ReKi ) ) ) .OR. &
              ( .NOT.( EqualRealNos( dt , 0.0_ReKi ) ) ) ) THEN   
              ! create a new property set 
              ! k, E, G, rho, d, t, Init
              
              kprop = kprop + 1
              CALL GetNewProp(kprop, TempProps(Prop1, 2), TempProps(Prop1, 3),&
                              Init%PropSets(Prop1, 4), d1 + J*dd, t1 + J*dt, &
                              TempProps, TempNProp, Init%PropSetsCol)           
              kelem = kelem + 1
              CALL GetNewElem(kelem, knode-1, knode, nprop, kprop, p)
              nprop = kprop                
         ELSE
              kelem = kelem + 1
              CALL GetNewElem(kelem, knode-1, knode, nprop, nprop, p)                
               
         ENDIF
      ENDDO
      
      ! the element connect to Node2
      kelem = kelem + 1
      CALL GetNewElem(kelem, knode, Node2, nprop, Prop2, p)                

   ENDDO ! loop over all members

ELSE ! NDiv = 1

   Init%MemberNodes(1:p%NMembers, 1:2) = p%Elems(1:Init%NElem, 2:3)   

ENDIF ! if NDiv is greater than 1

! set the props in Init
ALLOCATE(Init%Props(kprop, Init%PropSetsCol), STAT=Sttus)
IF ( Sttus /= 0 )  THEN
   ErrMsg = ' Error allocating TempProps arrays'
   ErrStat = ErrID_Fatal
   RETURN
ENDIF
Init%NProp = kprop
Init%Props(1:kprop, 1:Init%PropSetsCol) = TempProps

! deallocate temp matrices
IF (ALLOCATED(TempProps)) DEALLOCATE(TempProps)
IF (ALLOCATED(TempMembers)) DEALLOCATE(TempMembers)
IF (ALLOCATED(TempReacts)) DEALLOCATE(TempReacts)

END SUBROUTINE SD_Discrt
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE GetNewNode(k, x, y, z, Init)
   USE NWTC_Library
   USE SubDyn_Types
   IMPLICIT NONE

   TYPE(SD_InitInputType)   ::Init
   
   ! local variables
   INTEGER          :: k
   REAL(ReKi)       :: x, y, z
   
   Init%Nodes(k, 1) = k
   Init%Nodes(k, 2) = x
   Init%Nodes(k, 3) = y
   Init%Nodes(k, 4) = z


END SUBROUTINE GetNewNode
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE GetNewElem(k, n1, n2, p1, p2, p)
   USE NWTC_Library
   USE SubDyn_Types
   IMPLICIT NONE

   INTEGER,                INTENT(IN   )   :: k
   INTEGER,                INTENT(IN   )   :: n1
   INTEGER,                INTENT(IN   )   :: n2
   INTEGER,                INTENT(IN   )   :: p1
   INTEGER,                INTENT(IN   )   :: p2
   TYPE(SD_ParameterType), INTENT(INOUT)   :: p
  
   
   ! Local Variables
   
   p%Elems(k, 1) = k
   p%Elems(k, 2) = n1
   p%Elems(k, 3) = n2
   p%Elems(k, 4) = p1
   p%Elems(k, 5) = p2

END SUBROUTINE GetNewElem
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE GetNewProp(k, E, G, rho, d, t, TempProps, NTempProps, PropCol)
   USE NWTC_Library
   USE SubDyn_Types
   IMPLICIT NONE

   TYPE(SD_InitInputType)   ::Init
   
   ! local variables
   INTEGER                  :: k, NTempProps, PropCol
   REAL(ReKi)               :: E, G, rho, d, t
   REAL(ReKi)               :: TempProps(NTempProps, PropCol)
   
   TempProps(k, 1) = k
   TempProps(k, 2) = E
   TempProps(k, 3) = G
   TempProps(k, 4) = rho
   TempProps(k, 5) = d
   TempProps(k, 6) = t

END SUBROUTINE GetNewProp
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!SUBROUTINE InitIAJA(Init)                                                                                                 
!   ! for 2-node element only                                                                                              
!   USE NWTC_Library                                                                                                       
!   USE SubDyn_Types                                                                                                       
!   USE qsort_c_module                                                                                                     
!   IMPLICIT NONE                                                                                                          
!                                                                                                                          
!   TYPE(SD_InitInputType)   ::Init                                                                                        
!                                                                                                                          
!   ! local variables                                                                                                      
!   INTEGER                      :: I, J, k, UnDbg, Sttus, ERRSTAT, r, s                                                   
!   INTEGER                      :: NNZ, TDOF                                                                              
!   INTEGER, ALLOCATABLE         :: IA(:), JA(:)                                                                           
!   CHARACTER(1024)              :: OutFile, tempStr, OutFmt                                                               
!   INTEGER, ALLOCATABLE         :: Col_Arr(:)                                                                             
!   INTEGER                      :: ND                                                                                     
!   INTEGER                      :: SortA(Init%MaxMemjnt,2)                                                                
!                                                                                                                          
!   ! allocate for NodesConnE and NodesConnN                                                                               
!   ALLOCATE(Init%NodesConnE(Init%NNode, Init%MaxMemJnt+2), STAT=Sttus)                                                    
!   IF ( Sttus /= 0 )  THEN                                                                                                
!      ErrMsg = ' Error allocating NodesConnE matrix'                                                    
!      ErrStat = 1                                                                                                         
!      RETURN                                                                                                              
!   ENDIF                                                                                                                  
!   Init%NodesConnE = 0                                                                                                    
!                                                                                                                          
!   ALLOCATE(Init%NodesConnN(Init%NNode, Init%MaxMemJnt+2), STAT=Sttus)                                                    
!   IF ( Sttus /= 0 )  THEN                                                                                                
!      ErrMsg = ' Error allocating NodesConnE matrix'                                                    
!      ErrStat = 1                                                                                                         
!      RETURN                                                                                                              
!   ENDIF                                                                                                                  
!   Init%NodesConnN = 0                                                                                                    
!                                                                                                                          
!   ! find the node connectivity, nodes/elements that connect to a common node                                             
!      ! initialize the temp array for sorting                                                                             
!   SortA(Init%MaxMemjnt,2) = 0                                                                                            
!                                                                                                                          
!   DO I = 1, Init%NNode                                                                                                   
!      Init%NodesConnE(I, 1) = Init%Nodes(I, 1)                                                                            
!      Init%NodesConnN(I, 1) = Init%Nodes(I, 1)                                                                            
!                                                                                                                          
!      k = 0                                                                                                               
!      DO J = 1, Init%NElem                                                                                                
!         IF ( Init%Nodes(I, 1)==p%Elems(J, 2) ) THEN                                                                   
!            k = k + 1                                                                                                     
!            Init%NodesConnE(I, k + 2) = p%Elems(J, 1)                                                                  
!            Init%NodesConnN(I, k + 2) = p%Elems(J, 3)                                                                  
!         ENDIF                                                                                                            
!         IF ( Init%Nodes(I, 1)==p%Elems(J, 3) ) THEN                                                                   
!            k = k + 1                                                                                                     
!            Init%NodesConnE(I, k + 2) = p%Elems(J, 1)                                                                  
!            Init%NodesConnN(I, k + 2) = p%Elems(J, 2)                                                                  
!         ENDIF                                                                                                            
!      ENDDO                                                                                                               
!                                                                                                                          
!      IF( k>1 )THEN ! sort the nodes ascendingly                                                                          
!         SortA(1:k, 1) = Init%NodesConnN(I, 3:(k+2))                                                                      
!         CALL QsortC( SortA(1:k, 1:2) )                                                                                   
!         Init%NodesConnN(I, 3:(k+2)) = SortA(1:k, 1)                                                                      
!      ENDIF                                                                                                               
!                                                                                                                          
!      Init%NodesConnE(I, 2) = k                                                                                           
!      Init%NodesConnN(I, 2) = k                                                                                           
!   ENDDO                                                                                                                  
!                                                                                                                          
!                                                                                                                          
!                                                                                                                          
!! allocate the column array - column numbers that have nonzero component                                                  
!   ALLOCATE(Col_Arr( 6*(Init%MaxMemJnt+1)), STAT=Sttus)                                                                   
!   IF ( Sttus /= 0 )  THEN                                                                                                
!      ErrMsg = ' Error allocating Col_Arr'                                                              
!      ErrStat = 1                                                                                                         
!      RETURN                                                                                                              
!   ENDIF                                                                                                                  
!   Col_Arr = 0                                                                                                            
!                                                                                                                          
!! allocate the row array IA                                                                                               
!   TDOF = 6*Init%NNode                                                                                                    
!   Init%TDOF = TDOF                                                                                                       
!                                                                                                                          
!   ALLOCATE(IA( TDOF + 1 ), STAT=Sttus)                                                                                   
!   IF ( Sttus /= 0 )  THEN                                                                                                
!      ErrMsg = ' Error allocating IA'                                                                   
!      ErrStat = 1                                                                                                         
!      RETURN                                                                                                              
!   ENDIF                                                                                                                  
!   IA = 0                                                                                                                 
!                                                                                                                          
!! allocate the row array JA                                                                                               
!   ALLOCATE(JA( (TDOF*TDOF+TDOF)/2 ), STAT=Sttus)                                                                         
!   IF ( Sttus /= 0 )  THEN                                                                                                
!      ErrMsg = ' Error allocating IA'                                                                   
!      ErrStat = 1                                                                                                         
!      RETURN                                                                                                              
!   ENDIF                                                                                                                  
!   JA = 0                                                                                                                 
!                                                                                                                          
!                                                                                                                          
!   ! for each node, find the nodes that connect to this node                                                              
!   ! find the target columns in the global K and M for these nodes                                                        
!   NNZ = 0                                                                                                                
!   DO I = 1, Init%NNode                                                                                                   
!                                                                                                                          
!      DO J = 1, 6 ! the common node: first 6 columns of row I                                                             
!         Col_Arr(J) = (I-1)*6 + J                                                                                         
!      ENDDO                                                                                                               
!                                                                                                                          
!      k = 1                                                                                                               
!      DO J = 1, Init%NodesConnN(I, 2)                                                                                     
!         nd = Init%NodesConnN(I, J+2) ! nodes that connect to the common node                                             
!         IF ( nd > I ) THEN ! only count the node number that is greater than the common node                             
!            k = k + 1                                                                                                     
!            DO r = 1, 6                                                                                                   
!               Col_Arr( (k-1)*6 + r ) = (nd-1)*6 + r                                                                      
!            ENDDO ! r                                                                                                     
!         ENDIF                                                                                                            
!                                                                                                                          
!      ENDDO ! J                                                                                                           
!                                                                                                                          
!!     write(*, *) ' Col_Arr '                                                                                             
!!     write(*, *) 'I = ', I                                                                                               
!!     WRITE(tempStr, '(I10)') k*6                                                                                         
!!     OutFmt = ('('//trim(tempStr)//'(I5))')                                                                              
!!     WRITE(*, trim(OutFmt)) (Col_Arr(1:k*6))                                                                             
!                                                                                                                          
!                                                                                                                          
!     ! total number of columns is k*6                                                                                     
!     ! only count the diagonal and the upper triangle                                                                     
!     DO r = 1, 6                                                                                                          
!        IA( (I-1)*6 + r ) = NNZ + 1                                                                                       
!        DO s = r, k*6                                                                                                     
!            NNZ = NNZ + 1                                                                                                 
!            JA(NNZ) = Col_Arr(s)                                                                                          
!        ENDDO ! s                                                                                                         
!                                                                                                                          
!     ENDDO ! r                                                                                                            
!                                                                                                                          
!   ENDDO ! I                                                                                                              
!   IA( TDOF + 1 ) = NNZ+1                                                                                                 
!                                                                                                                          
!                                                                                                                          
!! allocate the row array IA                                                                                               
!   ALLOCATE(Init%IA( TDOF +1 ), STAT=Sttus)                                                                               
!   IF ( Sttus /= 0 )  THEN                                                                                                
!      ErrMsg = ' Error allocating Init%IA'                                                              
!      ErrStat = 1                                                                                                         
!      RETURN                                                                                                              
!   ENDIF                                                                                                                  
!   Init%IA = IA                                                                                                           
!                                                                                                                          
!! allocate the row array JA                                                                                               
!   ALLOCATE(Init%JA(NNZ), STAT=Sttus)                                                                                     
!   IF ( Sttus /= 0 )  THEN                                                                                                
!      ErrMsg = ' Error allocating Init%JA'                                                              
!      ErrStat = 1                                                                                                         
!      RETURN                                                                                                              
!   ENDIF                                                                                                                  
!   Init%JA(1:NNZ) = JA(1:NNZ)                                                                                             
!   Init%NNZ = NNZ                                                                                                         
!                                                                                                                          
!! deallocate temp matrices                                                                                                
!IF (ALLOCATED(IA)) DEALLOCATE(IA)                                                                                         
!IF (ALLOCATED(JA)) DEALLOCATE(JA)                                                                                         
!IF (ALLOCATED(Col_Arr)) DEALLOCATE(Col_Arr)                                                                               
!                                                                                                                          
!! test the qsortc subroutine                                                                                              
!!TestA(:,1) = (/1,3,4,2,6/)                                                                                               
!!TestA(:,2) = (/1,3,4,2,6/)                                                                                               
!                                                                                                                          
!!CALL QsortC(TestA(1:5, 1:1))                                                                                             
!                                                                                                                          
!!write(*, '(2(I5))') ((TestA(i, j), j = 1, 2), i = 1, 5)                                                                  
!                                                                                                                          
!!--------------------------------------                                                                                   
!! write node connectivity data and IA, JA to a txt file                                                                   
!CALL GetNewUnit( UnDbg )                                                                                                  
!                                                                                                                          
!OutFile = (trim(Init%RootName)//'_Connectivity.txt' )                                                                     
!CALL OpenFOutFile ( UnDbg, OutFile , ErrStat )                                                                            
!                                                                                                                          
!IF ( ErrStat /= 0 ) THEN                                                                                                  
!   CLOSE( UnDbg )                                                                                                         
!   RETURN                                                                                                                 
!END IF                                                                                                                    
!                                                                                                                          
!                                                                                                                          
!WRITE(tempStr, '(I10)') Init%MaxMemJnt + 2                                                                                
!OutFmt = ('('//trim(tempStr)//'(I6))')                                                                                    
!WRITE(UnDbg, '(A)') 'Elements connect to a common node'                                                                   
!WRITE(UnDbg, trim(OutFmt)) ((Init%NodesConnE(i, j), j = 1, Init%MaxMemJnt+2), i = 1, Init%NNode)                          
!WRITE(UnDbg, '(A)') 'Nodes connect to a common node'                                                                      
!WRITE(UnDbg, trim(OutFmt)) ((Init%NodesConnN(i, j), j = 1, Init%MaxMemJnt+2), i = 1, Init%NNode)                          
!                                                                                                                          
!WRITE(UnDbg, *) 'TDOF = ', TDOF                                                                                           
!WRITE(UnDbg, '(I6, I6)') ( (i, Init%IA(i)), i = 1, TDOF+1)                                                                
!Write(UnDbg, *) 'NNZ = ', NNZ                                                                                             
!write(Undbg, '(I6, I6)') ( (i, Init%JA(i)), i = 1, NNZ)                                                                   
!                                                                                                                          
!CLOSE(UnDbg)                                                                                                              
!                                                                                                                          
!END SUBROUTINE InitIAJA                                                                                                   

!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE AssembleKM(Init,p, ErrStat, ErrMsg)

   USE NWTC_Library
   USE SubDyn_Types
   IMPLICIT NONE

   TYPE(SD_InitInputType),       INTENT(INOUT)  ::Init
   TYPE(SD_ParameterType),       INTENT(INOUT)  ::p
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(1024),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   
   !TYPE(SD_InitInputType),               INTENT(  OUT)     :: ElemProps(Init%NElem)
   
   INTEGER                  :: I, J, K, Jn, Kn
   
   INTEGER                  :: NNE        ! number of nodes in one element
   INTEGER                  :: N1, N2     ! starting node and ending node in the element
   INTEGER                  :: P1, P2     ! property set numbers for starting and ending nodes

   REAL(ReKi)               :: D1, D2, t1, t2, E, G, rho ! properties of a section
   REAL(ReKi)               :: x1, y1, z1, x2, y2, z2    ! coordinates of the nodes
   REAL(ReKi)               :: DirCos(3, 3)              ! direction cosine matrices
   REAL(ReKi)               :: L                         ! length of the element
   REAL(ReKi)               :: r1, r2, t, Iyy, Jzz, Ixx, A, kappa
   LOGICAL                  :: shear
   REAL(ReKi), ALLOCATABLE  :: Ke(:,:), Me(:, :), FGe(:) ! element stiffness and mass matrices gravity force vector
   INTEGER, ALLOCATABLE     :: nn(:)                     ! node number in element 
   INTEGER                  :: tgt_row_sys(6), row_in_elem(6), tgt_col_sys(6), col_in_elem(6)
   
   INTEGER                  :: ei, ej, ti, tj, r, s, beg_jA, end_jA, ii, jj
   CHARACTER(1024)          :: tempstr1, tempstr2, outfile
   INTEGER                  :: UnDbg, Sttus 
   
   

   
   
   
      ! for current application
   IF ( (Init%FEMMod .LE. 3) .and. (Init%FEMMod .GE. 0)) THEN
      NNE = 2   
   ENDIF                              
   
   ! total degree of freedoms of the system 
   Init%TDOF = 6*Init%NNode
   
   
   ! allocate element stiffness matrix
   ALLOCATE( Ke(NNE*6, NNE*6), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      ErrMsg = ' Error allocating element stiffness matrix Ke'
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF

   ! allocate element mass matrix
   ALLOCATE( Me(NNE*6, NNE*6), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      ErrMsg = ' Error allocating element mass matrix Me'
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF
   
      ! allocate element gravity force vector
   ALLOCATE( FGe(NNE*6), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      ErrMsg = ' Error allocating element gravity force vector FGe'
      ErrStat = 1
      RETURN
   ENDIF
   
   
   ! allocate system stiffness matrix
   ALLOCATE( Init%K(Init%TDOF, Init%TDOF), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      ErrMsg = ' Error allocating system stiffness matrix K'
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF
   Init%K = 0.0_ReKi

   ! allocate system mass matrix
   ALLOCATE( Init%M(Init%TDOF, Init%TDOF), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      ErrMsg = ' Error allocating element mass matrix M '
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF
   Init%M = 0.0_ReKi
   
      ! allocate system gravity force vector
   ALLOCATE( Init%FG(Init%TDOF), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      ErrMsg = ' Error allocating element gravity force vector FG'
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF
   Init%FG = 0.0_ReKi

   
   ! allocate node number in element array
   ALLOCATE( nn(NNE), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      ErrMsg = ' Error allocating nn'
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF
   
      ! loop over all elements
   DO I = 1, Init%NElem
   
      DO J = 1, NNE
         NN(J) = p%Elems(I, J + 1)
      ENDDO
   
      N1 = p%Elems(I,       2)
      N2 = p%Elems(I, NNE + 1)
      
      P1 = p%Elems(I, NNE + 2)
      P2 = p%Elems(I, NNE + 3)
      
      
      E   = Init%Props(P1, 2)
      G   = Init%Props(P1, 3)
      rho = Init%Props(P1, 4)
      D1  = Init%Props(P1, 5)
      t1  = Init%Props(P1, 6)
      D2  = Init%Props(P2, 5)
      t2  = Init%Props(P2, 6)
      
      x1  = Init%Nodes(N1, 2)
      y1  = Init%Nodes(N1, 3)
      z1  = Init%Nodes(N1, 4)
      
      x2  = Init%Nodes(N2, 2)
      y2  = Init%Nodes(N2, 3)
      z2  = Init%Nodes(N2, 4)

      CALL GetDirCos(X1, Y1, Z1, X2, Y2, Z2, DirCos, L, ErrStat, ErrMsg)
      IF ( ErrStat /= 0 ) RETURN
      CALL SetConstants()
      
      ! 1: uniform Euler-Bernouli
      ! 3: uniform Timoshenko
      IF ( (Init%FEMMod == 1).OR.(Init%FEMMod == 3)) THEN ! uniform element 
         r1 = 0.25*(D1 + D2)
         t  = 0.5*(t1+t2)
         
         IF ( EqualRealNos(t, 0.0_ReKi) ) THEN
            r2 = 0
         ELSE
            r2 = r1 - t
         ENDIF
         
         A = Pi_D*(r1*r1-r2*r2)
         Ixx = 0.25*Pi_D*(r1**4-r2**4)
         Iyy = Ixx
         Jzz = 2.0*Ixx
         
         IF( Init%FEMMod == 1 ) THEN ! uniform Euler-Bernoulli
            Shear = .false.
            kappa = 0
         ELSEIF( Init%FEMMod == 3 ) THEN ! uniform Timoshenko
            Shear = .true.
            kappa = 0.53
         ENDIF
         
         p%ElemProps(i)%Area = A
         p%ElemProps(i)%Length = L
         p%ElemProps(i)%Ixx = Ixx
         p%ElemProps(i)%Iyy = Iyy
         p%ElemProps(i)%Jzz = Jzz
         p%ElemProps(i)%Shear = Shear
         p%ElemProps(i)%kappa = kappa
         p%ElemProps(i)%YoungE = E
         p%ElemProps(i)%ShearG = G
         p%ElemProps(i)%Rho = rho
         p%ElemProps(i)%DirCos = DirCos
         
         
         CALL ElemK(A, L, Ixx, Iyy, Jzz, Shear, kappa, E, G, DirCos, Ke)
         CALL ElemM(A, L, Ixx, Iyy, Jzz, rho, DirCos, Me)
         CALL ElemG(A, L, rho, DirCos, FGe, Init%g)
         
      ELSEIF  (Init%FEMMod == 2) THEN ! tapered Euler-Bernoulli
      
         ErrStat = ErrID_Fatal
         ErrMsg = 'FEMMod = 2 is not implemented!'
         RETURN
         
      ELSEIF  (Init%FEMMod == 4) THEN ! tapered Timoshenko
         ErrStat = ErrID_Fatal
         ErrMsg = 'FEMMod = 4 is not implemented!'
         RETURN
         
      ELSE
         ErrStat = ErrID_Fatal
         ErrMsg = ' FEMMod is not valid! Please choose from 1, 2, 3, and 4. '
         RETURN
         
         
      ENDIF  
      
!      !~~~~~~~~~~ assemble to system K and M in compressed row format                                                                                            
!      !           for two node element                                                                                                                           
!      DO J = 1, NNE ! NNE = 2, 2 nodes in one element                                                                                                           
!                                                                                                                                                                
!         DO ja = 1, 6                                                                                                                                           
!            tgt_row_sys(ja) = ( nn(J) - 1 )*6 + ja                                                                                                              
!            row_in_elem(ja) = ( J - 1 )*6 + ja                                                                                                                  
!         ENDDO !ja                                                                                                                                              
!                                                                                                                                                                
!         DO K = 1, NNE                                                                                                                                          
!                                                                                                                                                                
!            DO ja = 1, 6                                                                                                                                        
!               tgt_col_sys(ja) = ( nn(K) - 1 )*6 + ja                                                                                                           
!               col_in_elem(ja) = ( K - 1 )*6 + ja                                                                                                               
!            ENDDO !ja                                                                                                                                           
!                                                                                                                                                                
!            DO ii = 1, 6                                                                                                                                        
!               ei = row_in_elem(ii)                                                                                                                             
!               ti = tgt_row_sys(ii)                                                                                                                             
!                                                                                                                                                                
!               DO jj = 1, 6                                                                                                                                     
!                  ej = col_in_elem(jj)                                                                                                                          
!                  tj = tgt_col_sys(jj)                                                                                                                          
!                                                                                                                                                                
!                  IF(ti .LE. tj) THEN ! store the upper triangle and the diagonal                                                                               
!                     beg_jA = Init%IA(ti)                                                                                                                       
!                     end_jA = Init%IA(ti+1) - 1                                                                                                                 
!                                                                                                                                                                
!                     s = 0                                                                                                                                      
!                     DO r = beg_jA, end_jA                                                                                                                      
!                        IF ( Init%JA(r) == tj ) THEN                                                                                                            
!                           s = r                                                                                                                                
!                        ENDIF                                                                                                                                   
!                     ENDDO ! r                                                                                                                                  
!                                                                                                                                                                
!                     IF ( s == 0) THEN                                                                                                                          
!                        write(tempstr1, *)  ti                                                                                                                  
!                        write(tempstr2, *)  tj                                                                                                                  
!                        ErrMsg = ('A( '//trim(tempstr1)//','//trim(tempstr2)//') not found in AssembleKM')                                    
!                        ErrStat = 1                                                                                                                             
!                        RETURN                                                                                                                                  
!                     ENDIF                                                                                                                                      
!                                                                                                                                                                
!                     Init%K(s) = Init%K(s) + Ke(ei, ej)                                                                                                         
!                     Init%M(s) = Init%M(s) + Me(ei, ej)                                                                                                         
!                                                                                                                                                                
!                                                                                                                                                                
!                  ENDIF                                                                                                                                         
!                                                                                                                                                                
!                                                                                                                                                                
!               ENDDO ! jj                                                                                                                                       
!            ENDDO !ii                                                                                                                                           
!                                                                                                                                                                
!         ENDDO ! K                                                                                                                                              
!                                                                                                                                                                
!      ENDDO ! J                                                                                                                                                 
!                                                                                                                                                                
!                                                                                                                                                                

      
      ! assemble element matrices to global matrices
         
      DO J = 1, NNE
         jn = nn(j)
         
         Init%FG( (jn*6-5):(jn*6) ) = Init%FG( (jn*6-5):(jn*6) ) &
                                    + FGe( (J*6-5):(J*6) )
         
         DO K = 1, NNE
            kn = nn(k)
            
            Init%K( (jn*6-5):(jn*6), (kn*6-5):(kn*6) ) = Init%K( (jn*6-5):(jn*6), (kn*6-5):(kn*6) ) &
                                                  + Ke( (J*6-5):(J*6), (K*6-5):(K*6) )
                  
            Init%M( (jn*6-5):(jn*6), (kn*6-5):(kn*6) ) = Init%M( (jn*6-5):(jn*6), (kn*6-5):(kn*6) ) &
                                                  + Me( (J*6-5):(J*6), (K*6-5):(K*6) )
                     
                     
         ENDDO !K
                     
      ENDDO !J
                     
                     
   ENDDO ! I end loop over elements
               
               
!   ! add concentrated mass (compressed row format)
!   DO I = 1, Init%NCMass
!      DO J = 1, 3
!          r = ( NINT( Init%CMass(I, 1) ) - 1 )*6 + J
!          Init%M(Init%IA(r)) = Init%M(Init%IA(r)) + Init%CMass(I, 2)
!      ENDDO
!      DO J = 4, 6
!          r = ( NINT( Init%CMass(I, 1) ) - 1 )*6 + J
!          Init%M(Init%IA(r)) = Init%M(Init%IA(r)) + Init%CMass(I, J-1)
!      ENDDO
!
!   ENDDO ! I concentrated mass
         
         
      
      ! add concentrated mass 
   DO I = 1, Init%NCMass
      DO J = 1, 3
          r = ( NINT(Init%CMass(I, 1)) - 1 )*6 + J
          Init%M(r, r) = Init%M(r, r) + Init%CMass(I, 2)
          
      ENDDO
      DO J = 4, 6
          r = ( NINT(Init%CMass(I, 1)) - 1 )*6 + J
          Init%M(r, r) = Init%M(r, r) + Init%CMass(I, J-1)
      ENDDO

   ENDDO ! I concentrated mass
 
      ! add concentrated mass induced gravity force
   DO I = 1, Init%NCMass
      
      r = ( NINT(Init%CMass(I, 1)) - 1 )*6 + 3
      Init%FG(r) = Init%FG(r) - Init%CMass(I, 2)*Init%g !TODO Changed this sign for concentrated load because now g is positive. GJH 5/6/13

   ENDDO ! I concentrated mass induced gravity
   
   

! deallocate temp matrices
IF (ALLOCATED(Ke)) DEALLOCATE(Ke)
IF (ALLOCATED(Me)) DEALLOCATE(Me)
IF (ALLOCATED(nn)) DEALLOCATE(nn)
   


END SUBROUTINE AssembleKM
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------

SUBROUTINE GetDirCos(X1, Y1, Z1, X2, Y2, Z2, DirCos, xyz, ErrStat, ErrMsg)
   !This should be from local to global -RRD
   !Convention used: keep x (local) in global X-Z plane in the general x positive direction  THIS IS THE OLD WAY (huimin's)
   USE NWTC_Library
   IMPLICIT NONE

   REAL(ReKi) , INTENT(IN)         :: x1, y1, z1, x2, y2, z2
   REAL(ReKi) , INTENT(OUT)        :: DirCos(3, 3)
   
   REAL(ReKi)         :: xz,  xyz, Dx,Dy,Dz, Delta,Dxy
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(1024),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   
   
   xz = sqrt( (x1-x2)**2 + (z1-z2)**2 )
   Dxy = sqrt( (x2-x1)**2 + (y2-y1)**2 )
   xyz = sqrt( (x1-x2)**2 + (z1-z2)**2 + (y1-y2)**2 )
   
   IF ( EqualRealNos(xyz, 0.0_ReKi) ) THEN
      ErrMsg = ' Same starting and ending location in the element.'
      ErrStat = 4
      RETURN
   ENDIF
   
   DirCos = 0
   
   IF ( EqualRealNos(xz, 0.0_ReKi) ) THEN
      IF( y2 < y1) THEN  !x is kept along global x
         DirCos(1, 1) = 1.0
         DirCos(2, 3) = -1.0
         DirCos(3, 2) = 1.0
      ELSE
         DirCos(1, 1) = 1.0
         DirCos(2, 3) = 1.0
         DirCos(3, 2) = -1.0
      ENDIF
 
   ELSE
 
      DirCos(1, 1) =  (z2-z1)/xz
      DirCos(1, 2) = -(x1-x2)*(y1-y2)/(xz*xyz)
      DirCos(1, 3) = (x2-x1)/xyz

      DirCos(2, 2) = xz/xyz
      DirCos(2, 3) = (y2-y1)/xyz

      DirCos(3, 1) = -(x2-x1)/xz
      DirCos(3, 2) = -(y1-y2)*(z1-z2)/(xz*xyz)
      DirCos(3, 3) = (z2-z1)/xyz
   ENDIF
      !RRD new DirCos
      DirCos=0
      Dy=y2-y1
      Dx=x2-x1
      Dz=z2-z1
      Delta=xyz
      IF ( EqualRealNos(Dxy, 0.0_ReKi) ) THEN !RRD
       IF( Dz < 0) THEN  !x is kept along global x
         DirCos(1, 1) =  1.0
         DirCos(2, 3) = -1.0
         DirCos(3, 3) = -1.0
      ELSE
         DirCos(1, 1) = 1.0
         DirCos(2, 2) = 1.0
         DirCos(3, 3) = 1.0
     
      ENDIF 
     ELSE
      DirCos(1, 1) =  Dy/Dxy
      DirCos(1, 2) = +Dx*Dz/(Delta*Dxy)
      DirCos(1, 3) =  Dx/Delta
      
      DirCos(2, 1) = -Dx/Dxy
      DirCos(2, 2) = +Dz*Dy/(Delta*Dxy)
      DirCos(2, 3) =  Dy/Delta
     
      DirCos(3, 2) = -Dxy/Delta
      DirCos(3, 3) = +Dz/Delta
     ENDIF
     !RRD end of new DirCos

END SUBROUTINE GetDirCos
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------

SUBROUTINE ElemK(A, L, Ixx, Iyy, Jzz, Shear, kappa, E, G, DirCos, K)
   ! element stiffness matrix for classical beam elements
   ! shear is true  -- non-tapered Timoshenko beam 
   ! shear is false -- non-tapered Euler-Bernoulli beam 
   USE NWTC_Library
   IMPLICIT NONE

   REAL(ReKi), INTENT( IN)               :: A, L, Ixx, Iyy, Jzz, E, G, kappa
   REAL(ReKi), INTENT( IN)               :: DirCos(3,3)
   LOGICAL, INTENT( IN)                  :: Shear
   
   REAL(ReKi), INTENT(OUT)             :: K(12, 12)  !RRD:  Ke and Me  need to be modified if convention of dircos is not followed?
         
   REAL(ReKi)                            :: Ax, Ay, Kx, Ky
   REAL(ReKi)                            :: DC(12, 12)
   
   Ax = kappa*A
   Ay = kappa*A
   
   K = 0
   
   IF (Shear) THEN
      Kx = 12.0*E*Iyy / (G*Ax*L*L)
      Ky = 12.0*E*Ixx / (G*Ay*L*L)
   ELSE
      Kx = 0.0
      Ky = 0.0
   ENDIF
      
   K( 9,  9) = E*A/L
   K( 7,  7) = 12.0*E*Iyy/( L*L*L*(1.0 + Kx) )
   K( 8,  8) = 12.0*E*Ixx/( L*L*L*(1.0 + Ky) )
   K(12, 12) = G*Jzz/L
   K(10, 10) = (4.0 + Ky)*E*Ixx / ( L*(1.0+Ky) )  
   K(11, 11) = (4.0 + Kx)*E*Iyy / ( L*(1.0+Kx) )
   K( 2,  4) = -6.*E*Ixx / ( L*L*(1.0+Ky) )
   K( 1,  5) =  6.*E*Iyy / ( L*L*(1.0+Kx) )
   K( 4, 10) = (2.0-Ky)*E*Ixx / ( L*(1.0+Ky) )
   K( 5, 11) = (2.0-Kx)*E*Iyy / ( L*(1.0+Kx) )
   
   K( 3,  3)  = K(9,9)
   K( 1,  1)  = K(7,7)
   K( 2,  2)  = K(8,8)
   K( 6,  6)  = K(12,12)
   K( 4,  4)  = K(10,10)
   K(5,5)  = K(11,11)
   K(4,2)  = K(2,4)
   K(5,1)  = K(1,5)
   K(10,4) = K(4,10)
   K(11,5) = K(5,11)
   K(12,6)= -K(6,6)
   K(10,2)=  K(4,2)
   K(11,1)=  K(5,1)
   K(9,3) = -K(3,3)
   K(7,1) = -K(1,1)
   K(8,2) = -K(2,2)
   K(6, 12) = -K(6,6)
   K(2, 10) =  K(4,2)
   K(1, 11) =  K(5,1)
   K(3, 9)  = -K(3,3)
   K(1, 7)  = -K(1,1)
   K(2, 8)  = -K(2,2)
   K(11,7) = -K(5,1)
   K(10,8) = -K(4,2)
   K(7,11) = -K(5,1)
   K(8,10) = -K(4,2)
   K(7,5) = -K(5,1)
   K(5,7) = -K(5,1)
   K(8,4) = -K(4,2)
   K(4,8) = -K(4,2)
   
   DC = 0
   DC( 1: 3,  1: 3) = DirCos
   DC( 4: 6,  4: 6) = DirCos
   DC( 7: 9,  7: 9) = DirCos
   DC(10:12, 10:12) = DirCos
   
   K = MATMUL( MATMUL(DC, K), TRANSPOSE(DC) )
   
   !write(*, *) K - TRANSPOSE(K)

END SUBROUTINE ElemK
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------

SUBROUTINE ElemM(A, L, Ixx, Iyy, Jzz, rho, DirCos, M)
   ! element mass matrix for classical beam elements

   USE NWTC_Library
   IMPLICIT NONE

   REAL(ReKi), INTENT( IN)               :: A, L, Ixx, Iyy, Jzz, rho
   REAL(ReKi), INTENT( IN)               :: DirCos(3,3)
   
   REAL(ReKi)             :: M(12, 12)
         
   REAL(ReKi)                            :: t, rx, ry, po
   REAL(ReKi)                            :: DC(12, 12)
   
   t = rho*A*L;
   rx = rho*Ixx;
   ry = rho*Iyy;
   po = rho*Jzz*L;   

   M = 0
   
      
   M( 9,  9) = t/3.0
   M( 7,  7) = 13.0*t/35.0 + 6.0*ry/(5.0*L)
   M( 8,  8) = 13.0*t/35.0 + 6.0*rx/(5.0*L)
   M(12, 12) = po/3.0
   M(10, 10) = t*L*L/105.0 + 2.0*L*rx/15.0
   M(11, 11) = t*L*L/105.0 + 2.0*L*ry/15.0
   M( 2,  4) = -11.0*t*L/210.0 - rx/10.0
   M( 1,  5) =  11.0*t*L/210.0 + ry/10.0
   M( 3,  9) = t/6.0
   M( 5,  7) =  13.*t*L/420. - ry/10.
   M( 4,  8) = -13.*t*L/420. + rx/10. 
   M( 6, 12) = po/6.
   M( 2, 10) =  13.*t*L/420. - rx/10. 
   M( 1, 11) = -13.*t*L/420. + ry/10.
   M( 8, 10) =  11.*t*L/210. + rx/10.
   M( 7, 11) = -11.*t*L/210. - ry/10. 
   M( 1,  7) =  9.*t/70. - 6.*ry/(5.*L)
   M( 2,  8) =  9.*t/70. - 6.*rx/(5.*L)
   M( 4, 10) = -L*L*t/140. - rx*L/30. 
   M( 5, 11) = -L*L*t/140. - ry*L/30.
   
   M( 3,  3) = M( 9,  9)
   M( 1,  1) = M( 7,  7)
   M( 2,  2) = M( 8,  8)
   M( 6,  6) = M(12, 12)
   M( 4,  4) = M(10, 10)
   M( 5,  5) = M(11, 11)
   M( 4,  2) = M( 2,  4)
   M( 5,  1) = M( 1,  5)
   M( 9,  3) = M( 3,  9)
   M( 7,  5) = M( 5,  7)
   M( 8,  4) = M( 4,  8)
   M(12,  6) = M( 6, 12)
   M(10,  2) = M( 2, 10)
   M(11,  1) = M( 1, 11)
   M(10,  8) = M( 8, 10)
   M(11,  7) = M( 7, 11)
   M( 7,  1) = M( 1,  7)
   M( 8,  2) = M( 2,  8)
   M(10,  4) = M( 4, 10)
   M(11,  5) = M( 5, 11)
   
   DC = 0
   DC( 1: 3,  1: 3) = DirCos
   DC( 4: 6,  4: 6) = DirCos
   DC( 7: 9,  7: 9) = DirCos
   DC(10:12, 10:12) = DirCos
   
   M = MATMUL( MATMUL(DC, M), TRANSPOSE(DC) )

END SUBROUTINE ElemM


!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------

SUBROUTINE ApplyConstr(Init,p)

   USE NWTC_Library
  
   USE SubDyn_Types
   
   IMPLICIT NONE

   TYPE(SD_InitInputType), INTENT(INOUT)   :: Init
   TYPE(SD_ParameterType), INTENT(IN)   :: p
   
   INTEGER                  :: I, J, k
   INTEGER                  :: bgn_j, end_j, row_n
   
   DO I = 1, p%NReact*6
      row_n = Init%BCs(I, 1)
      
      IF (Init%BCs(I, 2) == 1) THEN
         Init%K(row_n, :) = 0
         Init%K(:, row_n) = 0
         Init%K(row_n, row_n) = 1
      
         Init%M(row_n, :) = 0
         Init%M(:, row_n) = 0
         Init%M(row_n, row_n) = 0 !0.00001          !what is this???? I changed this to 0.  RRD 7/31
      ENDIF
      
   ENDDO ! I

      

END SUBROUTINE ApplyConstr
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE ElemG(A, L, rho, DirCos, F, g)
         
   REAL(ReKi), INTENT( OUT)           :: F(12)
   REAL(ReKi), INTENT( IN )           :: A
   REAL(ReKi), INTENT( IN )           :: g
   REAL(ReKi), INTENT( IN )           :: L
   REAL(ReKi), INTENT( IN )           :: rho
   REAL(ReKi), INTENT( IN )           :: DirCos(3, 3)

   REAL(ReKi)                         :: TempCoeff

   F = 0
   F(3) = -0.5*L*rho*A*g ! TODO: Check this, I changed the sign because the force should be negative.  GJH May 5
   F(9) = F(3)

   TempCoeff = 1.0/12.0*g*L*L*rho*A  !RRD : I am changing this to >0 sign 6/10/13
      
   !F(4) = TempCoeff*( DirCos(1, 3)*DirCos(2, 1) - DirCos(1, 1)*DirCos(2, 3) ) !These do not work if convnetion on z2>z1, x2>x1, y2>y1 are not followed as I have discovered 7/23
   !F(5) = TempCoeff*( DirCos(1, 3)*DirCos(2, 2) - DirCos(1, 2)*DirCos(2, 3) ) 
   
   !RRD attempt at new dircos which keeps x in the X-Y plane
         F(4) = -TempCoeff * SQRT(1-DirCos(3,3)**2) * DirCos(1,1)
         F(5) = -TempCoeff * SQRT(1-DirCos(3,3)**2) * DirCos(2,1)
    !RRD ends
   F(10) = -F(4)
   F(11) = -F(5)
   !F(12) is 0 for g along z alone

   
END SUBROUTINE ElemG

!------------------------------------------------------------------------------------------------------
SUBROUTINE LumpForces(Area1,Area2,crat,L,rho, g, DirCos, F)
         !This rountine calculates the lumped gravity forces at the nodes given the element geometry
         !It assumes a linear variation of the dimensions from node 1 to node 2, thus the area may be quadratically varying if crat<>1
   REAL(ReKi), INTENT( OUT)           :: F(12)
   REAL(ReKi), INTENT( IN )           :: Area1,Area2,crat !X-sectional areas at node 1 and node 2, t2/t1 thickness ratio
   REAL(ReKi), INTENT( IN )           :: g !gravity
   REAL(ReKi), INTENT( IN )           :: L !Length of element
   REAL(ReKi), INTENT( IN )           :: rho !density
   REAL(ReKi), INTENT( IN )           :: DirCos(3, 3)

    !LOCALS
   REAL(ReKi)                         :: TempCoeff,a0,a1,a2  !coefficients of the gravity quadratically distributed force

   
   !Calculate quadratic polynomial coefficients
   a0=A1
   a2=( (Area1+A2) - (Area1*crat+Area2/crat) )/L**2.  !*x**2
   a1= (Area2-Area1)/L -a2*L                       !*x                   
   
   !Now calculate the Lumped Forces
   F = 0
   F(3) = -(a0*L/2. +a1*L**2/6. +a2*L**3/12. )*rho*g  !Forces along z (must be negative on earth)
   F(9) = -(a0*L/2. +a1*L**2/3. +a2*L**3/4.  )*rho*g  !Forces along z (must be negative on earth)

   !Now calculate the Lumped Moments
   !HERE TO BE COMPLETED FOR THE BELOW
   TempCoeff = 1.0/12.0*g*L*L*rho*Area2  !RRD : I am changing this to >0 sign 6/10/13
      
   !F(4) = TempCoeff*( DirCos(1, 3)*DirCos(2, 1) - DirCos(1, 1)*DirCos(2, 3) ) !These do not work if convnetion on z2>z1, x2>x1, y2>y1 are not followed as I have discovered 7/23
   !F(5) = TempCoeff*( DirCos(1, 3)*DirCos(2, 2) - DirCos(1, 2)*DirCos(2, 3) ) 
   
   !RRD attempt at new dircos which keeps x in the X-Y plane
         F(4) = -TempCoeff * SQRT(1-DirCos(3,3)**2) * DirCos(1,1)
         F(5) = -TempCoeff * SQRT(1-DirCos(3,3)**2) * DirCos(2,1)
    !RRD ends
   F(10) = -F(4)
   F(11) = -F(5)
   !F(12) is 0 for g along z alone

   
END SUBROUTINE LumpForces

END MODULE SD_FEM
