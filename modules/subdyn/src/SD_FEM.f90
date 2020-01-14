!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!    This file is part of SubDyn.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!**********************************************************************************************************************************
MODULE SD_FEM
  USE NWTC_Library
  USE SubDyn_Types
  IMPLICIT NONE

  INTEGER,          PARAMETER  :: LAKi            = R8Ki                  ! Define the kind to be used for LAPACK routines for getting eigenvalues/vectors. Apparently there is a problem with SGGEV's eigenvectors
 
  INTEGER(IntKi),   PARAMETER  :: MaxMemjnt       = 10                    ! Maximum number of members at one joint
  INTEGER(IntKi),   PARAMETER  :: MaxOutChs       = 2000                  ! Max number of Output Channels to be read in
  INTEGER(IntKi),   PARAMETER  :: TPdofL          = 6                     ! 6 degrees of freedom (length of u subarray [UTP])
   
  ! values of these parameters are ordered by their place in SubDyn input file:
  INTEGER(IntKi),   PARAMETER  :: JointsCol       = 4                     ! Number of columns in Joints (JointID, JointXss, JointYss, JointZss)
  INTEGER(IntKi),   PARAMETER  :: ReactCol        = 7                     ! Number of columns in reaction dof array (JointID,RctTDxss,RctTDYss,RctTDZss,RctRDXss,RctRDYss,RctRDZss)
  INTEGER(IntKi),   PARAMETER  :: InterfCol       = 7                     ! Number of columns in interf matrix (JointID,ItfTDxss,ItfTDYss,ItfTDZss,ItfRDXss,ItfRDYss,ItfRDZss)
  INTEGER(IntKi),   PARAMETER  :: MaxNodesPerElem = 2                     ! Maximum number of nodes per element (currently 2)
  INTEGER(IntKi),   PARAMETER  :: MembersCol      = MaxNodesPerElem + 3   ! Number of columns in Members (MemberID,MJointID1,MJointID2,MPropSetID1,MPropSetID2,COSMID) 
  INTEGER(IntKi),   PARAMETER  :: PropSetsCol     = 6                     ! Number of columns in PropSets  (PropSetID,YoungE,ShearG,MatDens,XsecD,XsecT)  !bjj: this really doesn't need to store k, does it? or is this supposed to be an ID, in which case we shouldn't be storing k (except new property sets), we should be storing IDs
  INTEGER(IntKi),   PARAMETER  :: XPropSetsCol    = 10                    ! Number of columns in XPropSets (PropSetID,YoungE,ShearG,MatDens,XsecA,XsecAsx,XsecAsy,XsecJxx,XsecJyy,XsecJ0)
  INTEGER(IntKi),   PARAMETER  :: COSMsCol        = 10                    ! Number of columns in (cosine matrices) COSMs (COSMID,COSM11,COSM12,COSM13,COSM21,COSM22,COSM23,COSM31,COSM32,COSM33)
  INTEGER(IntKi),   PARAMETER  :: CMassCol        = 5                     ! Number of columns in Concentrated Mass (CMJointID,JMass,JMXX,JMYY,JMZZ)
  ! Indices in Members table
  INTEGER(IntKi),   PARAMETER  :: iMProp= 4 ! Index in Members table where the PropSet1 and 2 are stored

  INTEGER(IntKi),   PARAMETER  :: SDMaxInpCols    = MAX(JointsCol,ReactCol,InterfCol,MembersCol,PropSetsCol,XPropSetsCol,COSMsCol,CMassCol)


  INTERFACE FINDLOCI ! In the future, use FINDLOC from intrinsic
     MODULE PROCEDURE FINDLOCI_ReKi
     MODULE PROCEDURE FINDLOCI_IntKi
  END INTERFACE

CONTAINS
    
!> Maps nodes to elements 
!! allocate NodesConnE and NodesConnN                                                                               
SUBROUTINE NodeCon(Init,p, ErrStat, ErrMsg)
  USE qsort_c_module ,only: QsortC
  TYPE(SD_InitType),              INTENT( INOUT ) :: Init
  TYPE(SD_ParameterType),         INTENT( IN    ) :: p
  INTEGER(IntKi),                 INTENT(   OUT ) :: ErrStat     ! Error status of the operation
  CHARACTER(*),                   INTENT(   OUT ) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
  ! Local variables
  INTEGER(IntKi) :: SortA(MaxMemJnt,1)  !To sort nodes and elements
  INTEGER(IntKi) :: I,J,K  !counter
  
  ! The row index is the number of the real node, i.e. ID, 1st col has number of elements attached to node, and 2nd col has element numbers (up to 10)                                    
  CALL AllocAry(Init%NodesConnE, Init%NNode, MaxMemJnt+1,'NodesConnE', ErrStat, ErrMsg); if (ErrStat/=0) return;
  CALL AllocAry(Init%NodesConnN, Init%NNode, MaxMemJnt+2,'NodesConnN', ErrStat, ErrMsg); if (ErrStat/=0) return;
  Init%NodesConnE = 0                                                                                                    
  Init%NodesConnN = 0                                                                                                    
                                                                                                                          
   ! find the node connectivity, nodes/elements that connect to a common node                                             
   DO I = 1, Init%NNode                                                                                                   
      Init%NodesConnN(I, 1) = NINT( Init%Nodes(I, 1) )      !This should not be needed, could remove the extra 1st column like for the other array                                                                      
      k = 0                                                                                                               
      DO J = 1, Init%NElem                          !This should be vectorized                                                                      
         IF ( ( NINT(Init%Nodes(I, 1))==p%Elems(J, 2)) .OR. (NINT(Init%Nodes(I, 1))==p%Elems(J, 3) ) ) THEN   !If i-th nodeID matches 1st node or 2nd of j-th element                                                                   
            k = k + 1                                                                                                     
            if (k > MaxMemJnt+1) then 
               CALL SetErrStat(ErrID_Fatal, 'Maximum number of members reached on node'//trim(Num2LStr(NINT(Init%Nodes(I,1)))), ErrStat, ErrMsg, 'NodeCon');
            endif
            Init%NodesConnE(I, k + 1) = p%Elems(J, 1)                                                                  
            Init%NodesConnN(I, k + 1) = p%Elems(J, 3)                                                                  
            IF ( NINT(Init%Nodes(I, 1))==p%Elems(J, 3) ) Init%NodesConnN(I, k + 1) = p%Elems(J, 2)     !If nodeID matches 2nd node of element                                                                
         ENDIF                                                                                                            
      ENDDO                                                                                                               
                                                                                                                          
      IF( k>1 )THEN ! sort the nodes ascendingly                                                                          
         SortA(1:k, 1) = Init%NodesConnN(I, 3:(k+2))  
         CALL QsortC( SortA(1:k, 1:1) )                                                                                   
         Init%NodesConnN(I, 3:(k+2)) = SortA(1:k, 1)                                                                      
      ENDIF                                                                                                               
                                                                                                                          
      Init%NodesConnE(I, 1) = k    !Store how many elements connect i-th node in 2nd column                                                                                       
      Init%NodesConnN(I, 2) = k                                                                                           
   ENDDO                            

END SUBROUTINE NodeCon
!----------------------------------------------------------------------------
!>
! - Removes the notion of "ID" and use Index instead
! - Creates Nodes (use indices instead of ID), similar to Joints array
! - Creates Elems (use indices instead of ID)  similar to Members array
! - Updates Reacts (use indices instead of ID)
! - Updates Interf (use indices instead of ID)
SUBROUTINE SD_ReIndex_CreateNodesAndElems(Init,p, ErrStat, ErrMsg)
   TYPE(SD_InitType),            INTENT(INOUT)  ::Init
   TYPE(SD_ParameterType),       INTENT(INOUT)  ::p
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variable
   INTEGER                       :: I, n, iMem, iNode, JointID   
   CHARACTER(1024)               :: sType !< String for element type
   INTEGER(IntKi)                :: ErrStat2
   CHARACTER(1024)               :: ErrMsg2
   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL AllocAry(p%Elems,         Init%NElem,    MembersCol, 'p%Elems',         ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry(Init%Nodes,      Init%NNode,    JointsCol,  'Init%Nodes',      ErrStat2, ErrMsg2); if(Failed()) return

   ! --- Initialize Nodes
   Init%Nodes = 0   
   do I = 1,Init%NJoints
      Init%Nodes(I, 1) = I                ! JointID replaced by index I
      Init%Nodes(I, 2) = Init%Joints(I, 2)
      Init%Nodes(I, 3) = Init%Joints(I, 3)
      Init%Nodes(I, 4) = Init%Joints(I, 4)
   enddo

   ! --- Re-Initialize Reactions, pointing to index instead of JointID
   do I = 1, p%NReact
      JointID=p%Reacts(I,1)
      p%Reacts(I,1) = FINDLOCI(Init%Joints(:,1), JointID ) ! Replace JointID with Index
      if (p%Reacts(I,1)<=0) then
         CALL Fatal('Reaction joint table: line '//TRIM(Num2LStr(I))//' refers to JointID '//trim(Num2LStr(JointID))//' which is not in the joint list!')
         return
      endif
   enddo

   ! --- Re-Initialize interface joints, pointing to index instead of JointID
   Init%IntFc = 0
   do I = 1, Init%NInterf
      JointID=Init%Interf(I,1)
      Init%Interf(I,1) = FINDLOCI(Init%Joints(:,1), JointID )
      if (Init%Interf(I,1)<=0) then
         CALL Fatal('Interface joint table: line '//TRIM(Num2LStr(I))//' refers to JointID '//trim(Num2LStr(JointID))//' which is not in the joint list!')
         return
      endif
   enddo

   ! Change numbering in concentrated mass matrix
   do I = 1, Init%NCMass
      JointID = Init%CMass(I,1)
      Init%CMass(I,1) = FINDLOCI(Init%Joints(:,1), JointID )
      if (Init%Interf(I,1)<=0) then
         CALL Fatal('Concentrated mass table: line '//TRIM(Num2LStr(I))//' refers to JointID '//trim(Num2LStr(JointID))//' which is not in the joint list!')
         return
      endif
   enddo


   ! --- Initialize Elems, starting with each member as an element (we'll take NDiv into account later)
   p%Elems = 0
   ! --- Replacing "MemberID"  "JointID", and "PropSetID" by simple index in this tables
   DO iMem = 1, p%NMembers
      ! Column 1  : member index (instead of MemberID)
      p%Elems(iMem,     1)  = iMem
      ! Column 2-3: Joint index (instead of JointIDs)
      p%Elems(iMem,     1)  = iMem  ! NOTE: element/member number (not MemberID)
      do iNode=2,3
         p%Elems(iMem,iNode) = FINDLOCI(Init%Joints(:,1), Init%Members(iMem, iNode) ) 
         if (p%Elems(iMem,iNode)<=0) then
            CALL Fatal(' MemberID '//TRIM(Num2LStr(Init%Members(iMem,1)))//' has JointID'//TRIM(Num2LStr(iNode-1))//' = '// TRIM(Num2LStr(Init%Members(iMem, iNode)))//' which is not in the joint list!')
            return
         endif
      enddo
      ! Column 4-5: PropIndex 1-2 (instead of PropSetID1&2)
      ! NOTE: this index has different meaning depending on the member type !
      DO n=iMProp,iMProp+1

         sType='Member x-section property'
         p%Elems(iMem,n) = FINDLOCI(Init%PropSets(:,1), Init%Members(iMem, n) ) 

         if (p%Elems(iMem,n)<=0) then
            CALL Fatal('For MemberID '//TRIM(Num2LStr(Init%Members(iMem,1)))//'the PropSetID'//TRIM(Num2LStr(n-3))//' is not in the'//trim(sType)//' table!')
         endif
      END DO !n, loop through property ids         
   END DO !iMem, loop through members
    
   ! TODO in theory, we shouldn't need these anymore
   ! deallocate(Init%Members)
   ! deallocate(Init%Joints)
CONTAINS
   LOGICAL FUNCTION Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_ReIndex_CreateNodesAndElems') 
      Failed =  ErrStat >= AbortErrLev
   END FUNCTION Failed
   SUBROUTINE Fatal(ErrMsg_in)
      CHARACTER(len=*), intent(in) :: ErrMsg_in
      CALL SetErrStat(ErrID_Fatal, ErrMsg_in, ErrStat, ErrMsg, 'SD_ReIndex_CreateNodesAndElems');
   END SUBROUTINE Fatal
END SUBROUTINE SD_ReIndex_CreateNodesAndElems

!----------------------------------------------------------------------------
SUBROUTINE SD_Discrt(Init,p, ErrStat, ErrMsg)
   TYPE(SD_InitType),            INTENT(INOUT)  ::Init
   TYPE(SD_ParameterType),       INTENT(INOUT)  ::p
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variable
   INTEGER                       :: I, J, n, Node1, Node2, Prop1, Prop2   
   INTEGER                       :: NNE      ! number of nodes per element
   INTEGER                       :: MaxNProp
   REAL(ReKi), ALLOCATABLE       :: TempProps(:, :)
   INTEGER, ALLOCATABLE          :: TempMembers(:, :)
   INTEGER                       :: knode, kelem, kprop, nprop
   REAL(ReKi)                    :: x1, y1, z1, x2, y2, z2, dx, dy, dz, dd, dt, d1, d2, t1, t2
   LOGICAL                       :: found, CreateNewProp
   INTEGER(IntKi)                :: ErrStat2
   CHARACTER(ErrMsgLen)          :: ErrMsg2
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! number of nodes per element
   IF( ( Init%FEMMod >= 0 ) .and. (Init%FEMMod <= 3) ) THEN
      NNE = 2 
   ELSE
      CALL Fatal('FEMMod '//TRIM(Num2LStr(Init%FEMMod))//' not implemented.')
      RETURN
   ENDIF
   
   ! Total number of element   
   Init%NElem = p%NMembers*Init%NDiv  
   ! Total number of nodes - Depends on division and number of nodes per element
   Init%NNode = Init%NJoints + ( Init%NDiv - 1 )*p%NMembers 
   Init%NNode = Init%NNode + (NNE - 2)*Init%NElem 
   
   ! check the number of interior modes
   IF ( p%Nmodes > 6*(Init%NNode - Init%NInterf - p%NReact) ) THEN
      CALL Fatal(' NModes must be less than or equal to '//TRIM(Num2LStr( 6*(Init%NNode - Init%NInterf - p%NReact) )))
      RETURN
   ENDIF
   
   CALL AllocAry(Init%MemberNodes,p%NMembers,    Init%NDiv+1,'Init%MemberNodes',ErrStat2, ErrMsg2); if(Failed()) return ! for two-node element only, otherwise the number of nodes in one element is different
   CALL AllocAry(Init%IntFc,      6*Init%NInterf,2,          'Init%IntFc',      ErrStat2, ErrMsg2); if(Failed()) return

   ! --- Reindexing JointsID and MembersID into Nodes and Elems arrays
   ! NOTE: need NNode and NElem 
   CALL SD_ReIndex_CreateNodesAndElems(Init,p, ErrStat2, ErrMsg2);  if(Failed()) return
   
   ! --- Initialize boundary constraint vector - TODO: assumes order of DOF, NOTE: Needs Reindexing first
   CALL AllocAry(Init%BCs, 6*p%NReact, 2, 'Init%BCs', ErrStat2, ErrMsg2); if(Failed()) return
   CALL InitConstr(Init, p)
      
   ! --- Initialize interface constraint vector - TODO: assumes order of DOF, NOTE: Needs Reindexing first
   DO I = 1, Init%NInterf
      DO J = 1, 6
         Init%IntFc( (I-1)*6+J, 1) = (Init%Interf(I,1)-1)*6+J; ! TODO assumes order of DOF, and needs Interf1 reindexed
         Init%IntFc( (I-1)*6+J, 2) = Init%Interf(I, J+1);
      ENDDO
   ENDDO

    Init%MemberNodes = 0
    ! --- Setting up MemberNodes (And Elems, Props, Nodes if divisions)
    if (Init%NDiv==1) then
       ! NDiv = 1
       Init%MemberNodes(1:p%NMembers, 1:2) = p%Elems(1:Init%NElem, 2:3) 
       Init%NProp = Init%NPropSets

    else if (Init%NDiv > 1) then

       ! Discretize structure according to NDiv 
       ! - Elems is fully reinitialized, connectivity needs to be done again using SetNewElem
       ! - Nodes are not  reinitialized, but appended to NNodes
       ! 

       ! Initialize Temp arrays that will contain user inputs + input from the subdivided members
       !  We don't know how many properties will be needed, so allocated to size MaxNProp
       MaxNProp   = Init%NPropSets + Init%NElem*NNE ! Maximum possible number of property sets (temp): This is property set per element node, for all elements (bjj, added Init%NPropSets to account for possibility of entering many unused prop sets)
       CALL AllocAry(TempMembers, p%NMembers,    MembersCol , 'TempMembers', ErrStat2, ErrMsg2); if(Failed()) return
       CALL AllocAry(TempProps,  MaxNProp,      PropSetsCol,'TempProps',  ErrStat2, ErrMsg2); if(Failed()) return
       TempProps = -9999.
       TempMembers                    = p%Elems(1:p%NMembers,:)
       TempProps(1:Init%NPropSets, :) = Init%PropSets   

       kelem = 0
       knode = Init%NJoints
       kprop = Init%NPropSets
       DO I = 1, p%NMembers !the first p%NMembers rows of p%Elems contain the element information
          ! Member data
          Node1 = TempMembers(I, 2)
          Node2 = TempMembers(I, 3)
          Prop1 = TempMembers(I, iMProp  )
          Prop2 = TempMembers(I, iMProp+1)

          IF ( Node1==Node2 ) THEN
             CALL Fatal(' Same starting and ending node in the member.')
             RETURN
          ENDIF
  
          
          Init%MemberNodes(I,           1) = Node1
          Init%MemberNodes(I, Init%NDiv+1) = Node2
          
          IF  ( ( .not. EqualRealNos(TempProps(Prop1, 2),TempProps(Prop2, 2) ) ) &
           .OR. ( .not. EqualRealNos(TempProps(Prop1, 3),TempProps(Prop2, 3) ) ) &
           .OR. ( .not. EqualRealNos(TempProps(Prop1, 4),TempProps(Prop2, 4) ) ) )  THEN
          
             CALL Fatal(' Material E,G and rho in a member must be the same')
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
          
             ! If both dd and dt are 0, no interpolation is needed, and we can use the same property set for new nodes/elements. otherwise we'll have to create new properties for each new node
          CreateNewProp = .NOT. ( EqualRealNos( dd , 0.0_ReKi ) .AND.  EqualRealNos( dt , 0.0_ReKi ) )  
          
          ! node connect to Node1
          knode = knode + 1
          Init%MemberNodes(I, 2) = knode
          CALL SetNewNode(knode, x1+dx, y1+dy, z1+dz, Init)
              
          IF ( CreateNewProp ) THEN   
               ! create a new property set 
               ! k, E, G, rho, d, t, Init
               kprop = kprop + 1
               CALL SetNewProp(kprop, TempProps(Prop1, 2), TempProps(Prop1, 3), TempProps(Prop1, 4), d1+dd, t1+dt, TempProps)           
               kelem = kelem + 1
               CALL SetNewElem(kelem, Node1, knode, Prop1, kprop, p)  
               nprop = kprop              
          ELSE
               kelem = kelem + 1
               CALL SetNewElem(kelem, Node1, knode, Prop1, Prop1, p)                
               nprop = Prop1 
          ENDIF
          
          ! interior nodes
          DO J = 2, (Init%NDiv-1)
             knode = knode + 1
             Init%MemberNodes(I, J+1) = knode

             CALL SetNewNode(knode, x1 + J*dx, y1 + J*dy, z1 + J*dz, Init)
             
             IF ( CreateNewProp ) THEN   
                  ! create a new property set 
                  ! k, E, G, rho, d, t, Init               
                  kprop = kprop + 1
                  CALL SetNewProp(kprop, TempProps(Prop1, 2), TempProps(Prop1, 3),&
                                  Init%PropSets(Prop1, 4), d1 + J*dd, t1 + J*dt, &
                                  TempProps)           
                  kelem = kelem + 1
                  CALL SetNewElem(kelem, knode-1, knode, nprop, kprop, p)
                  nprop = kprop                
             ELSE
                  kelem = kelem + 1
                  CALL SetNewElem(kelem, knode-1, knode, nprop, nprop, p)         
                   
             ENDIF
          ENDDO
          
          ! the element connect to Node2
          kelem = kelem + 1
          CALL SetNewElem(kelem, knode, Node2, nprop, Prop2, p)                

       ENDDO ! loop over all members
       !
       Init%NProp = kprop

    ENDIF ! if NDiv is greater than 1

    ! set the props in Init
    CALL AllocAry(Init%Props, Init%NProp, PropSetsCol, 'Init%PropsBeams', ErrStat2, ErrMsg2); if(Failed()) return

    if (Init%NDiv==1) then
       Init%Props(1:Init%NProp, 1:PropSetsCol) = Init%PropSets(1:Init%NProp, 1:PropSetsCol)
    else if (Init%NDiv>1) then
       Init%Props(1:Init%NProp, 1:PropSetsCol) = TempProps(1:Init%NProp, 1:PropSetsCol)
    endif
    !Init%Props(1:kprop, 1:Init%PropSetsCol) = TempProps
    Init%Props = TempProps(1:Init%NProp, :)  !!RRD fixed it on 1/23/14 to account for NDIV=1

    CALL CleanUp_Discrt()

CONTAINS
   LOGICAL FUNCTION Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_Discrt') 
      Failed =  ErrStat >= AbortErrLev
      if (Failed) CALL CleanUp_Discrt()
   END FUNCTION Failed

   SUBROUTINE Fatal(ErrMsg_in)
      CHARACTER(len=*), intent(in) :: ErrMsg_in
      CALL SetErrStat(ErrID_Fatal, ErrMsg_in, ErrStat, ErrMsg, 'SD_Discrt');
      CALL CleanUp_Discrt()
   END SUBROUTINE Fatal

   SUBROUTINE CleanUp_Discrt()
      ! deallocate temp matrices
      IF (ALLOCATED(TempProps))   DEALLOCATE(TempProps)
      IF (ALLOCATED(TempMembers)) DEALLOCATE(TempMembers)
   END SUBROUTINE CleanUp_Discrt

END SUBROUTINE SD_Discrt


!> Returns index of val in Array (val is an integer!)
! NOTE: in the future use intrinsinc function findloc
FUNCTION FINDLOCI_ReKi(Array, Val) result(i)
   real(ReKi)    , dimension(:), intent(in) :: Array !< Array to search in
   integer(IntKi), intent(in)               :: val   !< Val
   integer(IntKi)                           :: i     !< Index of joint in joint table
   logical :: found
   i = 1
   do while ( i <= size(Array) )
      if ( Val == NINT(Array(i)) ) THEN
         return ! Exit when found
      else
         i = i + 1
      endif
   enddo
   i=-1
END FUNCTION
!> Returns index of val in Array (val is an integer!)
! NOTE: in the future use intrinsinc function findloc
FUNCTION FINDLOCI_IntKi(Array, Val) result(i)
   integer(IntKi), dimension(:), intent(in) :: Array !< Array to search in
   integer(IntKi), intent(in)               :: val   !< Val
   integer(IntKi)                           :: i     !< Index of joint in joint table
   logical :: found
   i = 1
   do while ( i <= size(Array) )
      if ( Val == Array(i) ) THEN
         return ! Exit when found
      else
         i = i + 1
      endif
   enddo
   i=-1
END FUNCTION


!------------------------------------------------------------------------------------------------------
!> Set properties of node k
SUBROUTINE SetNewNode(k, x, y, z, Init)
   TYPE(SD_InitType),      INTENT(INOUT) :: Init
   INTEGER,                INTENT(IN)    :: k
   REAL(ReKi),             INTENT(IN)    :: x, y, z
   
   Init%Nodes(k, 1) = k
   Init%Nodes(k, 2) = x
   Init%Nodes(k, 3) = y
   Init%Nodes(k, 4) = z

END SUBROUTINE SetNewNode

!------------------------------------------------------------------------------------------------------
!> Set properties of element k
SUBROUTINE SetNewElem(k, n1, n2, p1, p2, p)
   INTEGER,                INTENT(IN   )   :: k
   INTEGER,                INTENT(IN   )   :: n1
   INTEGER,                INTENT(IN   )   :: n2
   INTEGER,                INTENT(IN   )   :: p1
   INTEGER,                INTENT(IN   )   :: p2
   TYPE(SD_ParameterType), INTENT(INOUT)   :: p
   
   p%Elems(k, 1) = k
   p%Elems(k, 2) = n1
   p%Elems(k, 3) = n2
   p%Elems(k, iMProp  ) = p1
   p%Elems(k, iMProp+1) = p2

END SUBROUTINE SetNewElem

!------------------------------------------------------------------------------------------------------
!> Set material properties of element k
SUBROUTINE SetNewProp(k, E, G, rho, d, t, TempProps)
   INTEGER   , INTENT(IN)   :: k
   REAL(ReKi), INTENT(IN)   :: E, G, rho, d, t
   REAL(ReKi), INTENT(INOUT):: TempProps(:, :)
   
   TempProps(k, 1) = k
   TempProps(k, 2) = E
   TempProps(k, 3) = G
   TempProps(k, 4) = rho
   TempProps(k, 5) = d
   TempProps(k, 6) = t

END SUBROUTINE SetNewProp

!------------------------------------------------------------------------------------------------------
!> Assemble stiffness and mass matrix, and gravity force vector
SUBROUTINE AssembleKM(Init,p, ErrStat, ErrMsg)
   TYPE(SD_InitType),            INTENT(INOUT) :: Init
   TYPE(SD_ParameterType),       INTENT(INOUT) :: p
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! Local variables
   INTEGER                  :: I, J, K, Jn, Kn
   INTEGER, PARAMETER       :: NNE=2      ! number of nodes in one element, fixed to 2
   INTEGER                  :: N1, N2     ! starting node and ending node in the element
   INTEGER                  :: P1, P2     ! property set numbers for starting and ending nodes
   REAL(ReKi)               :: D1, D2, t1, t2, E, G, rho ! properties of a section
   REAL(ReKi)               :: x1, y1, z1, x2, y2, z2    ! coordinates of the nodes
   REAL(ReKi)               :: DirCos(3, 3)              ! direction cosine matrices
   REAL(ReKi)               :: L                         ! length of the element
   REAL(ReKi)               :: r1, r2, t, Iyy, Jzz, Ixx, A, kappa, nu, ratioSq, D_inner, D_outer
   LOGICAL                  :: shear
   REAL(ReKi), ALLOCATABLE  :: Ke(:,:), Me(:, :), FGe(:) ! element stiffness and mass matrices gravity force vector
   INTEGER, DIMENSION(NNE)  :: nn                        ! node number in element 
   INTEGER                  :: r
   INTEGER(IntKi)           :: ErrStat2
   CHARACTER(ErrMsgLen)     :: ErrMsg2

   
   ! for current application
   if    (Init%FEMMod == 2) THEN ! tapered Euler-Bernoulli
       CALL Fatal ('FEMMod = 2 is not implemented.')
       return
   elseif (Init%FEMMod == 4) THEN ! tapered Timoshenko
       CALL Fatal ('FEMMod = 2 is not implemented.')
       return
   elseif ((Init%FEMMod == 1) .or. (Init%FEMMod == 3)) THEN !
      ! 1: uniform Euler-Bernouli,  3: uniform Timoshenko
   else
       CALL Fatal('FEMMod is not valid. Please choose from 1, 2, 3, and 4. ')
       return
   endif
   
   ! total degrees of freedom of the system 
   Init%TDOF = 6*Init%NNode
   
   ALLOCATE( p%ElemProps(Init%NElem), STAT=ErrStat2)
   IF (ErrStat2 /= 0) THEN
       CALL Fatal('Error allocating p%ElemProps')
       return
   ENDIF

   CALL AllocAry( Ke,     NNE*6,         NNE*6 , 'Ke',      ErrStat2, ErrMsg2); if(Failed()) return; ! element stiffness matrix
   CALL AllocAry( Me,     NNE*6,         NNE*6 , 'Me',      ErrStat2, ErrMsg2); if(Failed()) return; ! element mass matrix 
   CALL AllocAry( FGe,    NNE*6,                 'FGe',     ErrStat2, ErrMsg2); if(Failed()) return; ! element gravity force vector 
   CALL AllocAry( Init%K, Init%TDOF, Init%TDOF , 'Init%K',  ErrStat2, ErrMsg2); if(Failed()) return; ! system stiffness matrix 
   CALL AllocAry( Init%m, Init%TDOF, Init%TDOF , 'Init%M',  ErrStat2, ErrMsg2); if(Failed()) return; ! system mass matrix 
   CALL AllocAry( Init%FG,Init%TDOF,             'Init%FG', ErrStat2, ErrMsg2); if(Failed()) return; ! system gravity force vector 
   Init%K  = 0.0_ReKi
   Init%M  = 0.0_ReKi
   Init%FG = 0.0_ReKi

   
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

      CALL GetDirCos(X1, Y1, Z1, X2, Y2, Z2, DirCos, L, ErrStat2, ErrMsg2); if(Failed()) return
         
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
       ! kappa = 0.53            
         ! equation 13 (Steinboeck et al) in SubDyn Theory Manual 
         nu = E / (2.0_ReKi*G) - 1.0_ReKi
         D_outer = 2.0_ReKi * r1  ! average (outer) diameter
         D_inner = D_outer - 2*t  ! remove 2x thickness to get inner diameter
         ratioSq = ( D_inner / D_outer)**2
         kappa =   ( 6.0 * (1.0 + nu) **2 * (1.0 + ratioSq)**2 ) &
                 / ( ( 1.0 + ratioSq )**2 * ( 7.0 + 14.0*nu + 8.0*nu**2 ) + 4.0 * ratioSq * ( 5.0 + 10.0*nu + 4.0 *nu**2 ) )
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
      
      ! assemble element matrices to global matrices
      DO J = 1, NNE
         jn = nn(j)
         Init%FG( (jn*6-5):(jn*6) ) = Init%FG( (jn*6-5):(jn*6) )  + FGe( (J*6-5):(J*6) )
         DO K = 1, NNE
            kn = nn(k)
            Init%K( (jn*6-5):(jn*6), (kn*6-5):(kn*6) ) = Init%K( (jn*6-5):(jn*6), (kn*6-5):(kn*6) ) + Ke( (J*6-5):(J*6), (K*6-5):(K*6) )
            Init%M( (jn*6-5):(jn*6), (kn*6-5):(kn*6) ) = Init%M( (jn*6-5):(jn*6), (kn*6-5):(kn*6) ) + Me( (J*6-5):(J*6), (K*6-5):(K*6) )
         ENDDO !K
      ENDDO !J
   ENDDO ! I end loop over elements
      
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
   ENDDO ! Loop on concentrated mass

   ! add concentrated mass induced gravity force
   DO I = 1, Init%NCMass
       r = ( NINT(Init%CMass(I, 1)) - 1 )*6 + 3
       Init%FG(r) = Init%FG(r) - Init%CMass(I, 2)*Init%g 
   ENDDO ! I concentrated mass induced gravity
   
   CALL CleanUp_AssembleKM()
   
CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AssembleKM') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call Cleanup_AssembleKM()
   END FUNCTION Failed
   
   SUBROUTINE Fatal(ErrMsg_in)
      character(len=*), intent(in) :: ErrMsg_in
      CALL SetErrStat(ErrID_Fatal, ErrMsg_in, ErrStat, ErrMsg, 'AssembleKM');
      CALL CleanUp_AssembleKM()
   END SUBROUTINE Fatal

   SUBROUTINE CleanUp_AssembleKM()
      IF(ALLOCATED(Ke )) DEALLOCATE(Ke )
      IF(ALLOCATED(Me )) DEALLOCATE(Me )
      IF(ALLOCATED(FGe)) DEALLOCATE(FGe)
   END SUBROUTINE CleanUp_AssembleKM
   
END SUBROUTINE AssembleKM

!------------------------------------------------------------------------------------------------------
!> Computes directional cosine matrix DirCos
!! rrd: This should be from local to global 
!! bjj: note that this is the transpose of what is normally considered the Direction Cosine Matrix  
!!      in the FAST framework. It seems to be used consistantly in the code (i.e., the transpose 
!!      of this matrix is used later).
SUBROUTINE GetDirCos(X1, Y1, Z1, X2, Y2, Z2, DirCos, L, ErrStat, ErrMsg)
   REAL(ReKi) ,      INTENT(IN   )  :: x1, y1, z1, x2, y2, z2  ! (x,y,z) positions of two nodes making up an element
   REAL(ReKi) ,      INTENT(  OUT)  :: DirCos(3, 3)            ! calculated direction cosine matrix
   REAL(ReKi) ,      INTENT(  OUT)  :: L                       ! length of element
   INTEGER(IntKi),   INTENT(  OUT)  :: ErrStat                 ! Error status of the operation
   CHARACTER(*),     INTENT(  OUT)  :: ErrMsg                  ! Error message if ErrStat /= ErrID_None
   REAL(ReKi)                       ::  Dx,Dy,Dz, Dxy          ! distances between nodes
   ErrMsg  = ""
   ErrStat = ErrID_None
   
   Dy=y2-y1
   Dx=x2-x1
   Dz=z2-z1
   Dxy = sqrt( Dx**2 + Dy**2 )
   L   = sqrt( Dx**2 + Dy**2 + Dz**2)
   
   IF ( EqualRealNos(L, 0.0_ReKi) ) THEN
      ErrMsg = ' Same starting and ending location in the element.'
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF
   
   IF ( EqualRealNos(Dxy, 0.0_ReKi) ) THEN 
      DirCos=0.0_ReKi    ! whole matrix set to 0
      IF ( Dz < 0) THEN  !x is kept along global x
         DirCos(1, 1) =  1.0_ReKi
         DirCos(2, 2) = -1.0_ReKi
         DirCos(3, 3) = -1.0_ReKi
      ELSE
         DirCos(1, 1) = 1.0_ReKi
         DirCos(2, 2) = 1.0_ReKi
         DirCos(3, 3) = 1.0_ReKi
      ENDIF 
   ELSE
      DirCos(1, 1) =  Dy/Dxy
      DirCos(1, 2) = +Dx*Dz/(L*Dxy)
      DirCos(1, 3) =  Dx/L
      
      DirCos(2, 1) = -Dx/Dxy
      DirCos(2, 2) = +Dz*Dy/(L*Dxy)
      DirCos(2, 3) =  Dy/L
     
      DirCos(3, 1) = 0.0_ReKi
      DirCos(3, 2) = -Dxy/L
      DirCos(3, 3) = +Dz/L
   ENDIF

END SUBROUTINE GetDirCos

!------------------------------------------------------------------------------------------------------
!> Element stiffness matrix for classical beam elements
!! shear is true  -- non-tapered Timoshenko beam 
!! shear is false -- non-tapered Euler-Bernoulli beam 
SUBROUTINE ElemK(A, L, Ixx, Iyy, Jzz, Shear, kappa, E, G, DirCos, K)
   REAL(ReKi), INTENT( IN) :: A, L, Ixx, Iyy, Jzz, E, G, kappa
   REAL(ReKi), INTENT( IN) :: DirCos(3,3)
   LOGICAL   , INTENT( IN) :: Shear
   REAL(ReKi), INTENT(OUT) :: K(12, 12) 
   ! Local variables
   REAL(ReKi)                            :: Ax, Ay, Kx, Ky
   REAL(ReKi)                            :: DC(12, 12)
   
   Ax = kappa*A
   Ay = kappa*A
   
   K(1:12,1:12) = 0
   
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
   
   K = MATMUL( MATMUL(DC, K), TRANSPOSE(DC) ) ! TODO: change me if DirCos convention is  transposed
   
END SUBROUTINE ElemK

!------------------------------------------------------------------------------------------------------
!> Element mass matrix for classical beam elements
SUBROUTINE ElemM(A, L, Ixx, Iyy, Jzz, rho, DirCos, M)
   REAL(ReKi), INTENT( IN) :: A, L, Ixx, Iyy, Jzz, rho
   REAL(ReKi), INTENT( IN) :: DirCos(3,3)
   REAL(ReKi), INTENT(OUT) :: M(12, 12)

   REAL(ReKi) :: t, rx, ry, po
   REAL(ReKi) :: DC(12, 12)
   
   t = rho*A*L;
   rx = rho*Ixx;
   ry = rho*Iyy;
   po = rho*Jzz*L;   

   M(1:12,1:12) = 0
      
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
   
   M = MATMUL( MATMUL(DC, M), TRANSPOSE(DC) ) ! TODO change me if direction cosine is transposed

END SUBROUTINE ElemM

!------------------------------------------------------------------------------------------------------
!> Sets a list of DOF indices and the value these DOF should have
!! NOTE: need p%Reacts to have an updated first column that uses indices and not JointID
SUBROUTINE InitConstr(Init, p)
   TYPE(SD_InitType     ),INTENT(INOUT) :: Init
   TYPE(SD_ParameterType),INTENT(IN   ) :: p
   !
   INTEGER(IntKi) :: I,J

   Init%BCs = 0
   DO I = 1, p%NReact
      DO J = 1, 6
         Init%BCs( (I-1)*6+J, 1) = (p%Reacts(I,1)-1)*6+J; ! DOF Index, looping through Joints in index order
         Init%BCs( (I-1)*6+J, 2) = p%Reacts(I, J+1);
      ENDDO
   ENDDO
END SUBROUTINE InitConstr

!> Apply constraint (Boundary conditions) on Mass and Stiffness matrices
SUBROUTINE ApplyConstr(Init,p)
   TYPE(SD_InitType     ),INTENT(INOUT):: Init
   TYPE(SD_ParameterType),INTENT(IN   ):: p
   
   INTEGER :: I !, J, k
   INTEGER :: row_n !bgn_j, end_j,
   
   DO I = 1, p%NReact*6
      row_n = Init%BCs(I, 1)
      IF (Init%BCs(I, 2) == 1) THEN
         Init%K(row_n,:    )= 0
         Init%K(:    ,row_n)= 0
         Init%K(row_n,row_n)= 1

         Init%M(row_n,:    )= 0
         Init%M(:    ,row_n)= 0
         Init%M(row_n,row_n)= 0
      ENDIF
   ENDDO ! I, loop on reaction nodes
END SUBROUTINE ApplyConstr

!------------------------------------------------------------------------------------------------------
!> calculates the lumped forces and moments due to gravity on a given element:
!! the element has two nodes, with the loads for both elements stored in array F. Indexing of F is:
!!    Fx_n1=1,Fy_n1=2,Fz_n1=3,Mx_n1= 4,My_n1= 5,Mz_n1= 6,
!!    Fx_n2=7,Fy_n2=8,Fz_n2=9,Mx_n2=10,My_n2=11,Mz_n2=12
SUBROUTINE ElemG(A, L, rho, DirCos, F, g)
   REAL(ReKi), INTENT( IN )           :: A            !< area
   REAL(ReKi), INTENT( IN )           :: L            !< element length
   REAL(ReKi), INTENT( IN )           :: rho          !< density
   REAL(ReKi), INTENT( IN )           :: DirCos(3, 3) !< direction cosine matrix (for determining distance between nodes 1 and 2)
   REAL(ReKi), INTENT( IN )           :: g            !< gravity
   REAL(ReKi), INTENT( OUT)           :: F(12)        !< returned loads. positions 1-6 are the loads for node 1; 7-12 are loads for node 2.
   REAL(ReKi) :: TempCoeff
   REAL(ReKi) :: w            ! weight per unit length
   
   F = 0             ! initialize whole array to zero, then set the non-zero portions
   w = rho*A*g       ! weight per unit length
   
   ! lumped forces on both nodes (z component only):
   F(3) = -0.5*L*w 
   F(9) = F(3)
          
   ! lumped moments on node 1 (x and y components only):
   ! bjj: note that RRD wants factor of 1/12 because of boundary conditions. Our MeshMapping routines use factor of 1/6 (assuming generic/different boundary  
   !      conditions), so we may have some inconsistent behavior. JMJ suggests using line2 elements for SubDyn's input/output meshes to improve the situation.
   TempCoeff = L*L*w/12.0_ReKi ! let's not calculate this twice  
   F(4) = -TempCoeff * DirCos(2,3) ! = -L*w*Dy/12.   !bjj: DirCos(2,3) = Dy/L
   F(5) =  TempCoeff * DirCos(1,3) ! =  L*w*Dx/12.   !bjj: DirCos(1,3) = Dx/L

      ! lumped moments on node 2: (note the opposite sign of node 1 moment)
   F(10) = -F(4)
   F(11) = -F(5)
   !F(12) is 0 for g along z alone
   
END SUBROUTINE ElemG
!------------------------------------------------------------------------------------------------------
!> Calculates the lumped gravity forces at the nodes given the element geometry
!! It assumes a linear variation of the dimensions from node 1 to node 2, thus the area may be quadratically varying if crat<>1
!! bjj: note this routine is a work in progress, intended for future version of SubDyn. Compare with ElemG.
SUBROUTINE LumpForces(Area1,Area2,crat,L,rho, g, DirCos, F)
   REAL(ReKi), INTENT( IN ) :: Area1,Area2,crat !< X-sectional areas at node 1 and node 2, t2/t1 thickness ratio
   REAL(ReKi), INTENT( IN ) :: g                !< gravity
   REAL(ReKi), INTENT( IN ) :: L                !< Length of element
   REAL(ReKi), INTENT( IN ) :: rho              !< density
   REAL(ReKi), INTENT( IN ) :: DirCos(3, 3)     !< Direction cosine matrix
   REAL(ReKi), INTENT( OUT) :: F(12)            !< Lumped forces
   !LOCALS
   REAL(ReKi)                         :: TempCoeff,a0,a1,a2  !coefficients of the gravity quadratically distributed force
   print*,'Error: the function lumpforces is not ready to use'
   STOP

   !Calculate quadratic polynomial coefficients
   a0 = -99999 ! TODO: this is wrong
   a2 = ( (Area1+A2) - (Area1*crat+Area2/crat) )/L**2. ! *x**2
   a1 = (Area2-Area1)/L -a2*L                          ! *x
   
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
   F(4) = -TempCoeff * SQRT(1-DirCos(3,3)**2) * DirCos(1,1) !bjj: compare with ElemG() and verify this lumping is consistent
   F(5) = -TempCoeff * SQRT(1-DirCos(3,3)**2) * DirCos(2,1) !bjj: compare with ElemG() and verify this lumping is consistent
   !RRD ends
   F(10) = -F(4)
   F(11) = -F(5)
   !F(12) is 0 for g along z alone
END SUBROUTINE LumpForces

END MODULE SD_FEM
